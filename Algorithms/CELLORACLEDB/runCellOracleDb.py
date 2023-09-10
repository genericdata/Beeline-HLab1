# refernece https://github.com/asistradition/inferelator_run_scripts/blob/master/gsj_2020/co_experimental.py
import gc
import copy
from optparse import OptionParser
import os
import sys
import pandas as pd
import warnings

from distutils.util import strtobool
from inferelator import utils
from inferelator import workflow
from inferelator import crossvalidation_workflow
from inferelator.benchmarking.celloracle import CellOracleWorkflow ,CellOracleRegression
from inferelator.distributed.inferelator_mp import MPControl

import numpy as np

import scanpy as sc
import celloracle as co

class CellOracleDbWorkflow(CellOracleWorkflow):

    _base_grn_file = None # str

    def set_base_grn_file(self,base_grn_file):
        self._base_grn_file = base_grn_file
    
    def startup_finish(self):
        """
        Skip inferelator preprocessing and do celloracle preprocessing; 
        and skip log transform

        As per https://github.com/morris-lab/CellOracle/issues/58
        """

        self.align_priors_and_expression()

        self.data.convert_to_float()

        adata = self.data._adata

        if "paga" not in adata.uns:
            utils.Debug.vprint("Normalizing data {sh}".format(sh=adata.shape))

            sc.pp.filter_genes(adata, min_counts=1)
            sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

            adata.raw = adata
            adata.layers["raw_count"] = adata.raw.X.copy()

            utils.Debug.vprint("Scaling data")

            #sc.pp.log1p(adata)
            sc.pp.scale(adata)

            utils.Debug.vprint("PCA Preprocessing")

            sc.tl.pca(adata, svd_solver='arpack')

            utils.Debug.vprint("Diffmap Preprocessing")

            sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
            sc.tl.diffmap(adata)
            sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

            utils.Debug.vprint("Clustering Preprocessing")

            sc.tl.louvain(adata, resolution=0.8)
            sc.tl.paga(adata, groups='louvain')
            sc.pl.paga(adata)
            sc.tl.draw_graph(adata, init_pos='paga', random_state=123)

            # Restore counts
            adata.X = adata.layers["raw_count"].copy()

        else:
            # Assume all the preprocessing is done and just move along

            utils.Debug.vprint("Using saved preprocessing for CellOracle")

class CellOracleDbRegression(CellOracleRegression):

    oracle_imputation = True

    def run_regression(self):

        utils.Debug.vprint("Creating Oracle Object")

        # Set up oracle object
        oracle = co.Oracle()
        oracle.import_anndata_as_raw_count(adata=self.data._adata,
                                           cluster_column_name="louvain",
                                           embedding_name="X_pca")

        # Apparently PCA is not transferred from the adata object
        oracle.perform_PCA(100)

        # Add base GRN
        base_grn = pd.read_parquet(self._base_grn_file)
        oracle.import_TF_data(TF_info_matrix=base_grn)
        
        utils.Debug.vprint("Imputation Preprocessing")

        if self.oracle_imputation:

            # Heuristics from Celloracle documentation
            n_comps = np.where(
                np.diff(
                    np.diff(
                        np.cumsum(oracle.pca.explained_variance_ratio_)
                    ) > 0.002
                )
            )[0][0]

            k = int(0.025 * oracle.adata.shape[0])

            # Make sure n_comps is between 10 and 50
            # It likes to go to 0 for noise controls
            n_comps = max(min(n_comps, 50), 10)

            # Make sure k is at least 25 too I guess
            k = max(k, 25)

            oracle.knn_imputation(
                n_pca_dims=n_comps,
                k=k, balanced=True,
                b_sight=k*8,
                b_maxl=k*4,
                n_jobs=4
            )

        # Pretend to do imputation
        else:
            oracle.adata.layers["imputed_count"] = oracle.adata.layers["normalized_count"].copy()

        utils.Debug.vprint("CellOracle GRN inference")

        # Call GRN inference
        links = oracle.get_links(cluster_name_for_GRN_unit="louvain", alpha=10,
                                 verbose_level=0, test_mode=False)

        # Deepcopy the result dict that we want
        result = copy.deepcopy(links.links_dict)

        # Try to clean up some of these circular references
        del links
        del oracle
        del self.data._adata
        del self.data

        # Call an explicit GC cycle
        gc.collect()

        return self.reprocess_co_output_to_inferelator_results(result)

            

def setup_workflow(in_dir, expr_file, tf_file, out_dir):
    
    worker = workflow.inferelator_workflow(regression=CellOracleDbRegression,workflow=CellOracleDbWorkflow)
    worker.set_file_paths(input_dir=in_dir, output_dir=out_dir, expression_matrix_file=expr_file,
                          tf_names_file=tf_file)
    worker.set_file_properties(expression_matrix_columns_are_genes=False)
    worker.set_file_loading_arguments('expression_matrix_file', sep="\t")
    worker.set_output_file_names(nonzero_coefficient_file_name=None, pdf_curve_file_name=None,
                                 curve_data_file_name=None)
    
    return worker

def set_up_cv_seeds(wkf,random_seed_list):
    cv = crossvalidation_workflow.CrossValidationManager(wkf)
    cv.add_gridsearch_parameter('random_seed',random_seed_list)
    return cv

def parseArgs(args):
    parser = OptionParser()

    parser.add_option('', '--in_dir', type='str',
                      help='Path to input directory')
    parser.add_option('', '--out_dir', type='str',
                      help='Path to output directory')
    parser.add_option('', '--expr_file',type='str',
                      help='File name of expression matrix file')
    parser.add_option('', '--gs_file',type='str',
                      help='File name of gold standard adjacency matrix. Do not specify this option to use no gold standard')
    parser.add_option('', '--regulator_file',type='str',
                      help='File name of regulator names')
    parser.add_option('', '--base_grn_file',type='str',
                      help='File name of base GRN file in parquest format'),
    parser.add_option('', '--split_gs_for_cv',type='str',default='True',
                      help='split_gold_standard_for_crossvalidation')
    parser.add_option('', '--cv_split_ratio',type='str',default='None',
                      help='CV split ratio. Ignored when --split_gs_for_cv is set to False.')
    parser.add_option('', '--random_seed_list',type=str,default='100',
                      help='comma-delimieted list of numbers as random seeds')    


    (opts, args) = parser.parse_args(args)

    return opts, args
            

def main(args):
    opts, args = parseArgs(args)

    MPControl.set_multiprocess_engine("local")
    MPControl.connect()
    utils.Debug.set_verbose_level(2)
    warnings.simplefilter("ignore")
    
    worker = setup_workflow(opts.in_dir,
                            opts.expr_file, opts.regulator_file, opts.out_dir)

    if (opts.gs_file is None):
        worker.set_network_data_flags(use_no_gold_standard=True)
        worker.set_file_paths(gold_standard_file=None)
    else:
        worker.set_file_paths(gold_standard_file=opts.gs_file)
    # Do not use priors; use the base GRN
    worker.set_network_data_flags(use_no_prior=True)
    worker.set_file_paths(priors_file=None)
    worker.set_base_grn_file(base_grn_file=opts.base_grn_file)

    worker.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=bool(strtobool(opts.split_gs_for_cv)),
                                          cv_split_ratio=None if opts.cv_split_ratio=='None' else float(opts.cv_split_ratio))
    
    random_seed_list = [int(seed) for seed in opts.random_seed_list.split(',')]
    cv_wrap = set_up_cv_seeds(worker,random_seed_list)
    cv_wrap.run()
    del cv_wrap
    
    #worker.run()
    del worker

                        
if __name__ == "__main__":
    main(sys.argv)


# refernece https://github.com/asistradition/inferelator_run_scripts/blob/668b889ced738287e045359357d2b367c3d0641d/gsj_2020/scenic_experimental.py
from optparse import OptionParser
import os
import sys
import pandas as pd
import warnings

from distutils.util import strtobool
from re import sub
from inferelator import utils
from inferelator import workflow
from inferelator import crossvalidation_workflow
from inferelator.benchmarking.scenic import SCENICWorkflow, SCENICRegression, ADJ_METHODS
from inferelator.distributed.inferelator_mp import MPControl

from arboreto.algo import grnboost2, genie3
import scanpy as sc

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df
from pyarrow.feather import write_feather


class SCENICWorkflow2(SCENICWorkflow):

    _filter_cells_gene_minimum = None #int; only keep cells with at least this many genes expressed
    _filter_genes_cell_minimum = None #int; only keep genes with at least this many cells expressed

    def set_filter_minimum(self,
                           filter_cells_gene_minimum=None,
                           filter_genes_cell_minimum=None):
        self._filter_cells_gene_minimum = filter_cells_gene_minimum
        self._filter_genes_cell_minimum = filter_genes_cell_minimum
    
    def startup_finish(self):

        self.align_priors_and_expression()

        tf_names = self.tf_names if self.tf_names is not None else self.priors_data.columns
        self.tf_names = [t for t in tf_names if t in self.data.gene_names]

        utils.Debug.vprint("Generating SCENIC prior files", level=0)

        self._feather_rank_file = self.create_feather_file_from_prior()
        self._motif_link_table_file = self.create_motif_table_from_prior()

        utils.Debug.vprint("feather_rank_file:",self._feather_rank_file, level=1)
        utils.Debug.vprint("motif_link_table_file:",self._motif_link_table_file, level=1)
        
        utils.Debug.vprint("Preprocessing data")
        utils.Debug.vprint(self.data.to_df(),level=1)
        if (self._filter_cells_gene_minimum != None):
            sc.pp.filter_cells(self.data._adata, min_genes=self._filter_cells_gene_minimum)
        if (self._filter_genes_cell_minimum != None):
            sc.pp.filter_genes(self.data._adata, min_cells=self._filter_genes_cell_minimum)

        self.data.convert_to_float()

        #sc.pp.normalize_per_cell(self.data._adata, counts_per_cell_after=1e4)
        #sc.pp.log1p(self.data._adata)
        sc.pp.scale(self.data._adata, max_value=10)
        utils.Debug.vprint(self.data.to_df(),level=0)


class SCENICRegression2(SCENICRegression):

    adjacency_method = "grnboost2"
    do_scenic = True

    _rank_threshold = None # int
    _nes_threshold = None # float

    def set_prune_thresholds(self,rank_threshold,nes_threshold):
        self._rank_threshold = rank_threshold
        self._nes_threshold = nes_threshold

    def run_regression(self):
        
        data_df = self.data.to_df()

        utils.Debug.vprint("Calculating {m} adjacencies".format(m=self.adjacency_method), level=0)

        # Get adjacencies
        adj_method = ADJ_METHODS[self.adjacency_method]

        if MPControl.is_dask():
            client_or_address = MPControl.client.client
            MPControl.client.check_cluster_state()
        else:
            client_or_address = 'local'

        adjacencies = adj_method(data_df, tf_names=self.tf_names, verbose=True, client_or_address=client_or_address,
                                 seed=self.random_seed)
        utils.Debug.vprint("Result {m} adjacencies".format(m=self.adjacency_method), adjacencies,level=0)

        if self.do_scenic:

            # Convert adjacencies to modules
            modules = list(modules_from_adjacencies(adjacencies, data_df))
            utils.Debug.vprint('modules:',modules[0].__repr__(), level=1) # ok
            
            # Load feather (rank) databases
            dbs = [RankingDatabase(fname = self._feather_rank_file, name = "RANKING_PRIOR")]
            utils.Debug.vprint('dbs:',dbs[0].genes)  # ok

            utils.Debug.vprint("Pruning adjacencies with SCENIC", level=0)

            # Prune to df
            if (self._rank_threshold is None):
                rank_threshold = dbs[0].total_genes - 1
            else:
                rank_threshold = self._rank_threshold
                
            df = prune2df(dbs, modules, self._motif_link_table_file, client_or_address=client_or_address,
                          rank_threshold=rank_threshold, nes_threshold=self._nes_threshold)
            utils.Debug.vprint('prune2df:',df, level=1) # ok
            utils.Debug.vprint('prune2df shape:',df.shape, level=1) # ok
            utils.Debug.vprint('priors_data:',self.priors_data,level=1)

            return self.reprocess_scenic_output_to_inferelator_results2(df, self.priors_data)

        else:

            return self.reprocess_adj_to_inferelator_results(adjacencies)


    @staticmethod
    def reprocess_scenic_output_to_inferelator_results2(scenic_df, prior_data):

        # if there's nothing in the scenic output make an empty dataframe of 0s
        if scenic_df.shape[0] == 0:
            mat = pd.DataFrame(0.0, index=prior_data.index, columns=prior_data.columns)

        else:
            scenic_df = scenic_df.copy()
            scenic_df.index = scenic_df.index.droplevel(1)
            scenic_df.columns = scenic_df.columns.droplevel(0)

            mat = [pd.DataFrame(data).set_index(0).rename({1: tf}, axis=1)
                    for tf, data in scenic_df['TargetGenes'].iteritems() if len(data)>0]

            if (len(mat)>0):
                mat = pd.concat(mat, axis=0).fillna(0)
                mat = mat.groupby(mat.index).agg('max')
                mat = mat.reindex(prior_data.columns, axis=1).reindex(prior_data.index, axis=0).fillna(0)
            else:
                mat = pd.DataFrame(0.0, index=prior_data.index, columns=prior_data.columns)
            print(mat)
    
        return [mat], [mat.copy()]

        
def setup_workflow(in_dir, expr_file, tf_file, out_dir):
    
    worker = workflow.inferelator_workflow(regression=SCENICRegression2, workflow=SCENICWorkflow2)
    worker.set_file_paths(input_dir=in_dir, output_dir=out_dir, expression_matrix_file=expr_file,
                          tf_names_file=tf_file)
    worker.set_file_properties(expression_matrix_columns_are_genes=False)
    worker.set_file_loading_arguments('expression_matrix_file', sep="\t")
    worker.set_output_file_names(nonzero_coefficient_file_name=None, pdf_curve_file_name=None,
                                 curve_data_file_name=None)
    
    return worker

def set_up_cv_seeds(wkf,random_seed_list):
    cv = crossvalidation_workflow.CrossValidationManager(wkf)
    cv.add_gridsearch_parameter('random_seed', random_seed_list)
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
    parser.add_option('', '--priors_file',type='str',
                      help='File name of prior interaction adjacency matrix. Do not specify this to use no priors')
    parser.add_option('', '--split_gs_for_cv',type='str',default='True',
                      help='split_gold_standard_for_crossvalidation')
    parser.add_option('', '--cv_split_ratio',type='str',default='None',
                      help='CV split ratio. Ignored when --split_gs_for_cv is set to False.')
    parser.add_option('', '--filter_cells_gene_minimum',type=int,default=None,
                      help='only keep cells with at least this many genes expressed')
    parser.add_option('', '--filter_genes_cell_minimum',type=int,default=None,
                      help='only keep genes with at least this many cells expressed')
    parser.add_option('', '--random_seed_list',type=str,default='100',
                      help='comma-delimieted list of numbers as random seeds')    
    parser.add_option('', '--rank_threshold',type=int,default=None,
                      help='rank_threshold for SCENIC')
    parser.add_option('', '--nes_threshold',type=float,default=None,
                      help='nes_threshold for SCENIC')

    
    (opts, args) = parser.parse_args(args)

    return opts, args
            

def main(args):
    opts, args = parseArgs(args)

    MPControl.set_multiprocess_engine("dask-local")
    MPControl.connect()
    utils.Debug.set_verbose_level(2)
    warnings.simplefilter("ignore")
    
    worker = setup_workflow(opts.in_dir,
                            opts.expr_file, opts.regulator_file, opts.out_dir)
    worker.dask_temp_path = os.environ['TMPDIR']

    if (opts.gs_file is None):
        worker.set_network_data_flags(use_no_gold_standard=True)
        worker.set_file_paths(gold_standard_file=None)
    else:
        worker.set_file_paths(gold_standard_file=opts.gs_file)
    if (opts.priors_file is None):
        worker.set_network_data_flags(use_no_prior=True)
        worker.set_file_paths(priors_file=None)
    else:
        worker.set_file_paths(priors_file=opts.priors_file)

    worker.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=bool(strtobool(opts.split_gs_for_cv)),
                                          cv_split_ratio=None if opts.cv_split_ratio=='None' else float(opts.cv_split_ratio))
    worker.set_prune_thresholds(rank_threshold=opts.rank_threshold,
                                nes_threshold=opts.nes_threshold)

    random_seed_list = [int(seed) for seed in opts.random_seed_list.split(',')]
    cv_wrap = set_up_cv_seeds(worker,random_seed_list)
    cv_wrap.run()
    del cv_wrap
    
    #worker.run()
    del worker

                        
if __name__ == "__main__":
    main(sys.argv)


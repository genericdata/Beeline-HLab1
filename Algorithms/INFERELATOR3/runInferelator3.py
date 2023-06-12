# refernece https://github.com/asistradition/inferelator_run_scripts/blob/1a4ab21c3e6fe3df1a28566786da90fbd468d418/gsj_2020/beeline_bbsr.py
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
from inferelator.postprocessing import MetricHandler
from inferelator.distributed.inferelator_mp import MPControl

def setup_workflow(method, wkf, in_dir, expr_file, tf_file, out_dir):
    
    worker = workflow.inferelator_workflow(regression=method, workflow=wkf)
    worker.set_file_paths(input_dir=in_dir, output_dir=out_dir, expression_matrix_file=expr_file,
                          tf_names_file=tf_file)
    worker.set_file_properties(expression_matrix_columns_are_genes=False)
    worker.set_file_loading_arguments('expression_matrix_file', sep="\t")
    worker.set_output_file_names(nonzero_coefficient_file_name=None, pdf_curve_file_name=None,
                                 curve_data_file_name=None)
    
    return worker

def parseArgs(args):
    parser = OptionParser()

    parser.add_option('', '--method', type='str',
                      help='bbsr etc')
    parser.add_option('', '--workflow', type='str',
                      help='tfa etc')
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
    parser.add_option('', '--random_seed',type='int',default=100,
                      help='random seed')
    parser.add_option('', '--num_bootstraps',type='int',default=5,
                      help='number of bootstraps')
    
    (opts, args) = parser.parse_args(args)

    return opts, args
            

def main(args):
    opts, args = parseArgs(args)

    MPControl.set_multiprocess_engine("local")
    MPControl.connect()
    utils.Debug.set_verbose_level(2)
    warnings.simplefilter("ignore")
    
    worker = setup_workflow(opts.method, opts.workflow, opts.in_dir,
                            opts.expr_file, opts.regulator_file, opts.out_dir)
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
    worker.set_run_parameters(num_bootstraps=opts.num_bootstraps,random_seed=opts.random_seed)
    worker.get_data()
    worker.run()

    del worker

                        
if __name__ == "__main__":
    main(sys.argv)



import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil

from inferelator import utils
from inferelator import workflow
from inferelator import crossvalidation_workflow
from inferelator.postprocessing import MetricHandler
from inferelator.distributed.inferelator_mp import MPControl

def generateInputs(RunnerObj):
    algName = RunnerObj.name
    runDir = RunnerObj.params.get('run_dir','')
    if not RunnerObj.inputDir.joinpath(algName,runDir).exists():
        print("Input folder for "+algName+" "+runDir+" does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath(algName,runDir).mkdir(exist_ok = False)

    print('Expression data: '+str(RunnerObj.inputDir.joinpath(RunnerObj.exprData)))
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                 header = 0, index_col = 0)
    if not RunnerObj.inputDir.joinpath(algName,runDir,"ExpressionData.csv").exists():
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath(algName,runDir,"ExpressionData.tsv"),
                             sep = '\t', header  = True, index = True)

    # reformat the true edges file to adj matrix
    print('True edges: '+str(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges)))
    if not RunnerObj.inputDir.joinpath(algName,runDir,"gold_standard.tsv").exists():
        gold_standard = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),sep=',')
        gold_standard['Type'] = gold_standard.get('Type',1)# some refNetworks.csv file do not have the type column
        gold_standard = gold_standard.pivot(index="Gene2", columns="Gene1")
        gold_standard[~pd.isna(gold_standard)] = 1.
        gold_standard = gold_standard.fillna(0).astype(int)
        gold_standard.columns = gold_standard.columns.droplevel(0)
        gold_standard = gold_standard.reindex(ExpressionData.index, axis=1).reindex(ExpressionData.index, axis=0).fillna(0).astype(int)
        gold_standard.to_csv(RunnerObj.inputDir.joinpath(algName,runDir,'gold_standard.tsv'),sep='\t')

    # reformat priors to adj matrix
    print('True edges prior: '+
          (RunnerObj.inputDir.joinpath(algName,runDir,RunnerObj.params['trueEdgesPrior']) if 'trueEdgesPrior' in RunnerObj.params else ''))
    if ('trueEdgesPrior' in RunnerObj.params) and (not RunnerObj.inputDir.joinpath(algName,runDir,"priors.tsv").exists()):
        priors = pd.read_csv(RunnerObj.inputDir.joinpath(algName,runDir,RunnerObj.params['trueEdgesPrior']),sep=',')
        priors['Type'] = priors.get('Type',1)
        priors = priors.pivot(index="Gene2", columns="Gene1")
        priors[~pd.isna(priors)] = 1.
        priors = priors.fillna(0).astype(int)
        priors.columns = priors.columns.droplevel(0)
        priors = priors.reindex(ExpressionData.index, axis=1).reindex(ExpressionData.index, axis=0).fillna(0).astype(int)
        priors.to_csv(RunnerObj.inputDir.joinpath(algName,runDir,'priors.tsv'),sep='\t')

def setup_workflow(in_dir, gs_file, tf_file, out_dir):
        worker = workflow.inferelator_workflow(regression=METHOD, workflow="tfa")
        worker.set_file_paths(input_dir=in_dir, output_dir=out_dir, expression_matrix_file=EXPR_FILE,
                              gold_standard_file=gs_file, tf_names_file=tf_file)
        worker.set_file_properties(expression_matrix_columns_are_genes=False)
        worker.set_file_loading_arguments('expression_matrix_file', sep=",")
        worker.set_run_parameters(num_bootstraps=5)
        worker.set_output_file_names(nonzero_coefficient_file_name=None, pdf_curve_file_name=None,
                                     curve_data_file_name=None)

        return worker
        
def run(RunnerObj):
    
    worker = setup_workflow(sub_pp, sub_pp_gs_reformatted, sub_pp_tfs, out_dir)
    worker.append_to_path('output_dir', "prior")
    
    worker.set_file_paths(priors_file=sub_pp_gs_reformatted)
    worker.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=True, cv_split_ratio=0.5)
    worker.set_run_parameters(random_seed=k)
    worker.get_data()
                                          
def parseOutput(RunnerObj):
    deepdrimRunner.parseOutputForPredictFromBestModel(RunnerObj,'DEEPDRIM7',RunnerObj.inputDir,'auprc')
                                          

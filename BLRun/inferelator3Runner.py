import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil

add_bind = ''
def generateInputs(RunnerObj):
    algName = RunnerObj.name
    runDir = RunnerObj.params.get('run_dir','')
    if not RunnerObj.inputDir.joinpath(algName,runDir).exists():
        print("Input folder for "+algName+" "+runDir+" does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath(algName,runDir).mkdir(exist_ok = False)

    print('Expression data: ',str(RunnerObj.inputDir.joinpath(RunnerObj.exprData)))
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                 header = 0, index_col = 0)
    if not RunnerObj.inputDir.joinpath(algName,runDir,"ExpressionData.tsv").exists():
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath(algName,runDir,"ExpressionData.tsv"),
                             sep = '\t', header  = True, index = True)

    # reformat the true edges file to adj matrix
    print('Gold standard edges: ',RunnerObj.params['useGS'] if RunnerObj.params['useGS']=='False' else str(RunnerObj.inputDir.joinpath(algName,runDir,RunnerObj.params['useGS'])))
    if (RunnerObj.params['useGS']!='False'):
        gold_standard = pd.read_csv(RunnerObj.inputDir.joinpath(algName,runDir,RunnerObj.params['useGS']),sep='\t')
        gold_standard['Type'] = gold_standard.get('Type',1)# some refNetworks.csv file do not have the type column
        gold_standard = gold_standard.drop_duplicates(subset=['Gene1','Gene2','Type']).pivot(index="Gene2", columns="Gene1", values="Type")
        gold_standard[~pd.isna(gold_standard)] = 1.
        gold_standard = gold_standard.fillna(0).astype(int)
        #gold_standard.columns = gold_standard.columns.droplevel(0)
        gold_standard = gold_standard.reindex(ExpressionData.index, axis=1).reindex(ExpressionData.index, axis=0).fillna(0).astype(int)
        gold_standard.to_csv(RunnerObj.inputDir.joinpath(algName,runDir,'gold_standard.tsv'),sep='\t')

    # reformat priors to adj matrix
    print('Prior edges: ',RunnerObj.params['usePrior'] if RunnerObj.params['usePrior']=='False' else str(RunnerObj.inputDir.joinpath(algName,runDir,RunnerObj.params['usePrior'])))
    if (RunnerObj.params['usePrior']!='False'):
        priors = pd.read_csv(RunnerObj.inputDir.joinpath(algName,runDir,RunnerObj.params['usePrior']),sep='\t')
        priors['Type'] = priors.get('Type',1)
        priors = priors.drop_duplicates(subset=['Gene1','Gene2','Type']).pivot(index="Gene2", columns="Gene1", values="Type")
        priors[~pd.isna(priors)] = 1.
        priors = priors.fillna(0).astype(int)
        #priors.columns = priors.columns.droplevel(0)
        priors = priors.reindex(ExpressionData.index, axis=1).reindex(ExpressionData.index, axis=0).fillna(0).astype(int)
        priors.to_csv(RunnerObj.inputDir.joinpath(algName,runDir,'priors.tsv'),sep='\t')

def run(RunnerObj):
    '''
    Function to run INFERELATOR3 algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
    runDir = RunnerObj.params.get('run_dir','')
    algName = RunnerObj.name
    inputDirB = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),algName,runDir)

    # outputs/L0/hESC/DEEPDRIM
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"+('' if runDir=='' else runDir+"/")
    outDirB = Path('/ext3/data').joinpath(outDir)
    os.makedirs(outDir, exist_ok = True)
    
    cmdToRun = ' '.join([
        'singularity exec --no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate INFERELATOR3; ',
        'cd',str(Path('/ext3/data')),';'
        'command time -v -o '+str(outDirB.joinpath('time.txt')),
        'python '+str(Path('/ext3/data').joinpath('Algorithms/INFERELATOR3/runInferelator3.py')),
        '--method ',RunnerObj.params['method'],
        '--workflow ',RunnerObj.params['workflow'],
        '--in_dir ',str(inputDirB),
        '--out_dir ',str(outDirB),
        '--expr_file ',str(inputDirB.joinpath('ExpressionData.tsv')),
        '--gs_file '+str(inputDirB.joinpath('gold_standard.tsv')) if RunnerObj.params['useGS']!='False' else '',
        '--regulator_file ',str(inputDirB.joinpath(RunnerObj.params['regulators'])),
        '--priors_file '+str(inputDirB.joinpath('priors.tsv')) if RunnerObj.params['usePrior']!='False' else '',
        '--split_gs_for_cv '+RunnerObj.params['splitGSForCV'],
        '--cv_split_ratio '+RunnerObj.params['CVSplitRatio'],
        '--random_seed '+RunnerObj.params['randomSeed'],
        '--num_bootstraps '+RunnerObj.params['numBootstraps'],
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    
def parseOutput(RunnerObj):
    '''
    Function to parse outputs from Inferelator3. 

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    runDir = RunnerObj.params.get('run_dir','')
    algName = RunnerObj.name
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"+('' if runDir=='' else runDir+"/")

    # Read output and reformat
    outFile = Path(outDir).joinpath('network.tsv.gz')
    if outFile.exists():
        netDF = pd.read_csv(outFile,compression='gzip',sep='\t',header=0)    
        outDF = netDF[['regulator','target','combined_confidences']].copy()
        outDF = outDF.sort_values(by=['combined_confidences','target','regulator'], ascending=[False,True,True])
        outDF.to_csv(Path(outDir).joinpath('rankedEdges.csv'),
                     sep='\t',index=False,header=['Gene1','Gene2','EdgeWeight'])
                                          

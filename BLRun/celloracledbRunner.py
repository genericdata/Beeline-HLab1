import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from BLRun import scenicdbRunner,inferelator3Runner

add_bind = ''

BASE_GRN_FNAMES = {'promoter_base_GRN':{
    "mm9_gimmemotifsv5_fpr1": "resources/celloracle_data/promoter_base_GRN/mm9_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
    "mm9_gimmemotifsv5_fpr2": "resources/celloracle_data/promoter_base_GRN/mm9_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
    "mm10_gimmemotifsv5_fpr1": "resources/celloracle_data/promoter_base_GRN/mm10_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
    "mm10_gimmemotifsv5_fpr2": "resources/celloracle_data/promoter_base_GRN/mm10_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
    "hg19_gimmemotifsv5_fpr1": "resources/celloracle_data/promoter_base_GRN/hg19_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
    "hg19_gimmemotifsv5_fpr2": "resources/celloracle_data/promoter_base_GRN/hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
    "hg38_gimmemotifsv5_fpr1": "resources/celloracle_data/promoter_base_GRN/hg38_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
    "hg38_gimmemotifsv5_fpr2": "resources/celloracle_data/promoter_base_GRN/hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"
}
}

def generateInputs(RunnerObj):

    algName = RunnerObj.name
    datasetName = RunnerObj.datasetName
    baseGRNType = RunnerObj.dataset_params['baseGRNType']
    baseGRNFName = BASE_GRN_FNAMES[baseGRNType][RunnerObj.dataset_params['baseGRNName'][datasetName]]
    
    runDir = RunnerObj.params.get('run_dir','')

    if not RunnerObj.inputDir.joinpath(algName,runDir).exists():
        print("Input folder for "+algName+" "+runDir+" does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath(algName,runDir).mkdir(exist_ok = False)

    print('Regulator ID data: ',RunnerObj.inputDir.joinpath(algName,RunnerObj.params['regulators']))
    print('mapRegulatorIDs:',RunnerObj.params['mapRegulatorIds'])
    exprDataIdMapOut,exprDataOut,tfNamesMapOut,tfNamesOut =\
        scenicdbRunner.processExpressionData(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                             RunnerObj.inputDir.joinpath(RunnerObj.exprDataIdMap),
                                             RunnerObj.exprDataIdMapType,
                                             [baseGRNFName+'.genes.'+RunnerObj.dataset_params['baseGRNIdMapType'][datasetName]],
                                             RunnerObj.dataset_params['baseGRNIdMapType'][datasetName],
                                             RunnerObj.inputDir.joinpath(algName,RunnerObj.params['regulators']),
                                             RunnerObj.params['mapRegulatorIds'])
    exprDataIdMapOut.to_csv(RunnerObj.inputDir.joinpath(algName,"ExpressionData-IDMap.tsv"),
                            sep='\t',header=True,index=False)
    exprDataOut.to_csv(RunnerObj.inputDir.joinpath(algName,"ExpressionData.tsv"),
                       sep = '\t', header  = True, index = True,index_label=False)
    if (RunnerObj.params['mapRegulatorIds']):
        tfNamesMapOut.to_csv(RunnerObj.inputDir.joinpath(algName,"TFs-IDMap.tsv"),
                             sep='\t',header=True,index=False)
        tfNamesOut.to_csv(RunnerObj.inputDir.joinpath(algName,"TFs.txt"),
                          sep = '\t', header=False, index =False)
    else:
        shutil.copyfile(RunnerObj.inputDir.joinpath(algName,RunnerObj.params['regulators']),
                        RunnerObj.inputDir.joinpath(algName,'TFs.txt'))

        
def run(RunnerObj):
    '''
    Function to run CellOracle algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
    algName = RunnerObj.name
    datasetName = RunnerObj.datasetName
    baseGRNType = RunnerObj.dataset_params['baseGRNType']
    baseGRNFName = BASE_GRN_FNAMES[baseGRNType][RunnerObj.dataset_params['baseGRNName'][datasetName]]
    
    runDir = RunnerObj.params.get('run_dir','')
    inputDirB = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),algName,runDir)

    # outputs/L0/hESC/DEEPDRIM
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"+('' if runDir=='' else runDir+"/")
    outDirB = Path('/ext3/data').joinpath(outDir)
    os.makedirs(outDir, exist_ok = True)

    cmdToRun = ' '.join([
        'singularity exec',
        '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate CELLORACLE; ',
        'cd',str(Path('/ext3/data')),';'
        'NUMBA_CACHE_DIR='+os.environ['TMPDIR'],
        'command time -v -o '+str(outDirB.joinpath('time.txt')),
        'python '+str(Path('/ext3/data').joinpath('Algorithms/CELLORACLEDB/runCellOracleDb.py')),
        '--in_dir ',str(inputDirB),
        '--out_dir ',str(outDirB),
        '--expr_file ',str(inputDirB.joinpath('ExpressionData.tsv')),
        '--gs_file '+str(inputDirB.joinpath('gold_standard.tsv')) if RunnerObj.params['useGS']!='False' else '',
        '--regulator_file ',str(inputDirB.joinpath('TFs.txt')),
        '--base_grn_file '+baseGRNFName,
        '--split_gs_for_cv '+RunnerObj.params['splitGSForCV'],
        '--cv_split_ratio '+RunnerObj.params['CVSplitRatio'],
        '--random_seed_list '+RunnerObj.params['randomSeedList'],
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from CellOracle. 

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    algName = RunnerObj.name

    random_seed_list = [int(seed) for seed in RunnerObj.params['randomSeedList'].split(',')]
    for seed in random_seed_list:
        runDir = 'random_seed_'+str(seed)
        outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"+runDir+"/"
        RunnerObj.params['run_dir'] = runDir
        inferelator3Runner.parseOutput(RunnerObj)

        outFile = Path(outDir).joinpath('rankedEdges.csv')
        if outFile.exists():
            shutil.copy(Path(outDir).joinpath('rankedEdges.csv'),
                        Path(outDir).joinpath('rankedEdges_Db.csv'))
            scenicdbRunner.processRankedEdgesDb(Path(outDir).joinpath('rankedEdges_Db.csv'),
                                                pd.read_csv(RunnerObj.inputDir.joinpath(algName,'ExpressionData-IDMap.tsv'),sep='\t',header=0),
                                                Path(outDir).joinpath('rankedEdges.csv'))


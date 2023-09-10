import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from BLRun import inferelator3Runner

add_bind = ''
def generateInputs(RunnerObj):
    inferelator3Runner.generateInputs(RunnerObj)
        
def run(RunnerObj):
    '''
    Function to run SCENIC algorithm

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
        '/bin/bash -c \"source /ext3/env.sh; conda activate SCENIC; ',
        'cd',str(Path('/ext3/data')),';'
        'NUMBA_CACHE_DIR='+os.environ['TMPDIR'],
        'command time -v -o '+str(outDirB.joinpath('time.txt')),
        'python '+str(Path('/ext3/data').joinpath('Algorithms/SCENIC/runScenic.py')),
        '--in_dir ',str(inputDirB),
        '--out_dir ',str(outDirB),
        '--expr_file ',str(inputDirB.joinpath('ExpressionData.tsv')),
        '--gs_file '+str(inputDirB.joinpath('gold_standard.tsv')) if RunnerObj.params['useGS']!='False' else '',
        '--regulator_file ',str(inputDirB.joinpath(RunnerObj.params['regulators'])),
        '--priors_file '+str(inputDirB.joinpath('priors.tsv')) if RunnerObj.params['usePrior']!='False' else '',
        '--split_gs_for_cv '+RunnerObj.params['splitGSForCV'],
        '--cv_split_ratio '+RunnerObj.params['CVSplitRatio'],
        '--random_seed_list '+RunnerObj.params['randomSeedList'],
        '--rank_threshold '+RunnerObj.params['rankThreshold'] if RunnerObj.params['rankThreshold']!='None' else '',
        '--nes_threshold '+RunnerObj.params['NESThreshold'] if RunnerObj.params['NESThreshold']!='None' else '',
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    
def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCENIC. 

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    random_seed_list = [int(seed) for seed in RunnerObj.params['randomSeedList'].split(',')]

    for seed in random_seed_list:
        RunnerObj.params['run_dir'] = 'random_seed_'+str(seed)
        inferelator3Runner.parseOutput(RunnerObj)

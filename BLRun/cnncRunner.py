import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for DeepDRIM.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("CNNC").exists():
        print("Input folder for CNNC does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("CNNC").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("CNNC/training_pairsCNNC.txt").exists():
        # input expression data, genes x cells
        print('Expression data: '+str(RunnerObj.inputDir.joinpath(RunnerObj.exprData)))

        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        # CNNC needs DataFrame of cells x genes in h5 file with key 'rpkm'
        ExpressionDataT = ExpressionData.T
        ExpressionDataT.to_hdf(RunnerObj.inputDir.joinpath('CNNC/ExpressionDataT.h5'),
                               key='rpkm',mode='w')
        
        geneNameMap = pd.concat([ExpressionData.index.to_series().str.lower(),ExpressionData.index.to_series()],axis=1)
        geneNameMap.to_csv(RunnerObj.inputDir.joinpath("CNNC/CNNC_geneName_map.txt"),
                           sep='\t',header=False,index=False)
        
        # Read file for trueEdges and create positive pairs file (no header)
        print('True edges: '+str(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges)))
        trueEdgesDF = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                sep = ',', header = 0, index_col = None)
        trueEdgesDF.to_csv(RunnerObj.inputDir.joinpath("CNNC/PositivePairs.txt"),
                                 sep=',',header=False,index=False)

        # relative to current path, current path is mounted as /ext3/data
        #inputDir = str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + "/CNNC"
        inputDir = RunnerObj.inputDir.relative_to(Path.cwd()).joinpath('CNNC')

        # Create training pairs, training pairs TF divide, gene name map files
        # training_pairsDEEPDRIM.txt, training_pairsDEEPDRIM.txtTF_divide_pos.txt, DEEPDRIM_geneName_map.txt
        # the script outputs to cwd, so enter into the CNNC directory to have the output files there
        cmdToRun = ' '.join([
            'singularity exec --no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/',
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate CNNC; ',
            'cd',str(Path('/ext3/data').joinpath(inputDir)),';',
            'python /ext3/DeepDRIM/generate_training_realdata.py ',
            '-expr_file '+str(RunnerObj.inputDir.joinpath(RunnerObj.exprData)),
            '-pos_pair_file '+'PositivePairs.txt',
            '-label DEEPDRIM'
            '\"'])
        print(cmdToRun)
        os.system(cmdToRun)

        trainPairsDF = pd.read_csv(RunnerObj.inputDir.joinpath('CNNC/training_pairsDEEPDRIM.txt'),
                                   sep='\t',header=None,names=['gene1','gene2','label'])
        label_new = trainPairsDF.loc[:,'label'].copy()
        for row_i in range(len(trainPairsDF)):
            gene1,gene2,label = trainPairsDF.iloc[row_i]
            if (label==2): # this pair is a reverse of a positive pair, check if this pair itself is positive
                pair_rev = trainPairsDF[(trainPairsDF['gene1']==gene2) & (trainPairsDF['gene2']==gene1)].reset_index(drop=True)
                if (pair_rev.empty): # reverse of this pair does not exist
                    label_new[row_i] = 0
                else: # reverse of this pair exists
                    pair_rev_label = pair_rev.at[0,'label']
                    if (pair_rev_label==1): # reverse of this pair is positive, set this pair to 0
                        label_new[row_i] = 0
                    elif (pair_rev_label==2):# reverse of this pair is a reverse of a positive pair, so this pair is positive
                        label_new[row_i] = 1
                    else:
                        label_new[row_i] = 0
        trainPairsDF_new = pd.concat([trainPairsDF['gene1'],trainPairsDF['gene2'],label_new],axis=1)
        trainPairsDF_new.to_csv(RunnerObj.inputDir.joinpath("CNNC/training_pairsCNNC.txt"),
                           sep='\t',header=False,index=False)

def run(RunnerObj):
    '''
    Function to run CNNC algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/CNNC'
    inputDir = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),'CNNC')
    # make output dirs if they do not exist:
    # outputs/L0/hESC/CNNC
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/CNNC/"
    os.makedirs(outDir, exist_ok = True)
    
    with open(RunnerObj.inputDir.joinpath('CNNC/training_pairsDEEPDRIM.txtTF_divide_pos.txt'), 'r') as fp:
        numBatches = len(fp.readlines()) -1

    print('Generating NEPDF: '+str(outDir+'NEPDF_data'))
    if Path(outDir+'NEPDF_data').exists():
        print("NEPDF_data folder for CNNC exists, removing NEPDF_data folder...")
        shutil.rmtree(outDir+'NEPDF_data')
    
    cmdToRun = ' '.join([
        'singularity exec --no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/',
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM;'
        'cd',str(Path('/ext3/data').joinpath(outDir)),';',
        'command time -v -o time_generate.txt',
        'python /ext3/CNNC/get_xy_label_data_cnn_combine_from_database.py',
        'None',
        str(RunnerObj.inputDir.joinpath('CNNC/CNNC_geneName_map.txt')),
        str(RunnerObj.inputDir.joinpath('CNNC/training_pairsDEEPDRIM.txt')),
        str(RunnerObj.inputDir.joinpath('CNNC/training_pairsDEEPDRIM.txtTF_divide_pos.txt')),
        'None',
        str(RunnerObj.inputDir.joinpath('CNNC/ExpressionDataT.h5')),
        '1',
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    print('Training on entire dataset')
    cmdToRun = ' '.join([
        'singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
        '--no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/',
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM;'
        'cd',str(Path('/ext3/data').joinpath(outDir)),';',
        'command time -v -o time_train.txt',
        'python /ext3/CNNC/train_new_model/train_with_labels_wholedatax.py',
        #'python /ext3/CNNC/train_new_model/train_with_labels_three_foldx.py',
        str(numBatches),
        'NEPDF_data',
        '3',
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    print('Predict on entire dataset')
    cmdToRun = ' '.join([
        'singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
        '--no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/',
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM;'
        'cd',str(Path('/ext3/data').joinpath(outDir)),';',
        'command time -v -o time_predict.txt',
        'python /ext3/CNNC/predict_no_y.py',
        str(numBatches),
        'NEPDF_data',
        '3',
        'xwhole_saved_models_T_32-32-64-64-128-128-512_e200/keras_cnn_trained_model_shallow.h5',
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)
    
def parseOutput(RunnerObj):
    '''
    Function to parse outputs from CNNC.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/CNNC/"

    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)
    
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.iterrows():
        outFile.write('\t'.join([row['TF'],row['target'],str(row['importance'])])+'\n')
    outFile.close()
    

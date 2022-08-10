import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from sklearn.model_selection import GroupKFold

add_bind = ''
def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for DeepDRIM.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("DEEPDRIM").exists():
        print("Input folder for DEEPDRIM does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("DEEPDRIM").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("DEEPDRIM/toPredict_pairsDEEPDRIM.txtTF_divide_pos.txt").exists():
        # input expression data
        print('Expression data: '+str(RunnerObj.inputDir.joinpath(RunnerObj.exprData)))
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        newExpressionData = ExpressionData.copy()
        
        # output expression .csv file (has header)
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)

        geneNameMap = pd.concat([ExpressionData.index.to_series().str.lower(),ExpressionData.index.to_series()],axis=1)
        geneNameMap.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/GeneName_map.txt"),
                           sep='\t',header=False,index=False)
        
        # Read file for trueEdges and create positive pairs file (no header)
        print('True edges: '+str(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges)))
        trueEdgesDF = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                sep = ',', header = 0, index_col = None)
        trueEdgesDF.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/PositivePairs.txt"),
                                 sep=',',header=False,index=False)
        
        
        # relative to current path; current path is mounted as /ext3/data
        #inputDir = str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + "/DEEPDRIM"
        inputDir = RunnerObj.inputDir.relative_to(Path.cwd()).joinpath('DEEPDRIM')

        # Create training pairs, training pairs TF divide, gene name map files
        # the function outputs to cwd, so enter into the DEEPDRIM directory to have the output files there
        cmdToRun = ' '.join([
            'singularity exec --no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
            'cd',str(Path('/ext3/data').joinpath(inputDir)),';',
            'python /ext3/DeepDRIM/generate_training_realdata.py ',
            '-expr_file '+'ExpressionData.csv',
            '-pos_pair_file '+'PositivePairs.txt',
            '-label DEEPDRIM'
            '\"'])
        print(cmdToRun)
        os.system(cmdToRun)

        # Create the cross-validation fold divide file
        # n lines - dividing positions for n-1 factors
        tfDividePos = np.genfromtxt(RunnerObj.inputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtTF_divide_pos.txt'),
                                    delimiter=',',dtype='int')
        tfGroups = np.repeat(np.arange(len(tfDividePos)-1),np.diff(tfDividePos))
        group_kfold = GroupKFold(n_splits=3)
        with open(RunnerObj.inputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide.txt'),'w') as ofh:
            split_it = group_kfold.split(X=np.zeros(len(tfGroups)),
                                         y=np.zeros(len(tfGroups)),groups=tfGroups)
            for train,test in split_it:
                print(np.unique(tfGroups[test]))
                ofh.write(','.join(map(str,np.unique(tfGroups[test]))))# comma seperated list of TF indexes
                ofh.write('\n')

        # Create the pairs to predict file; TF->all genes
        trainingPairs = pd.read_csv(RunnerObj.inputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txt'),
                                    sep='\t',header=None,index_col=None)
        trainingEdgesFrom = trainingPairs.iloc[tfDividePos[:-1],0]
        allGenes = ExpressionData.index.to_series().str.lower()
        pairsToPredict = pd.DataFrame(np.array(np.meshgrid(trainingEdgesFrom, allGenes, indexing='xy')).T.reshape(-1,2),
                                      columns=['Gene1','Gene2'])
        pairsToPredict['value'] = 3
        pairsToPredict.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/toPredict_pairsDEEPDRIM.txt"),
                                      sep='\t',header=False,index=False)

        pairsToPredict['rnum'] = pairsToPredict.index
        predictDividePos = [min(pairsToPredict.loc[pairsToPredict.Gene1==tf,'rnum']) for tf in trainingEdgesFrom]
        predictDividePos.append(pairsToPredict.shape[0])
        np.savetxt(RunnerObj.inputDir.joinpath("DEEPDRIM/toPredict_pairsDEEPDRIM.txtTF_divide_pos.txt"),predictDividePos,fmt='%d')
        
def run(RunnerObj):
    '''
    Function to run DEEPDRIM algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
    inputDir = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),'DEEPDRIM')

    # outputs/L0/hESC/DEEPDRIM
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPDRIM/"
    os.makedirs(outDir, exist_ok = True)
    
    with open(RunnerObj.inputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtTF_divide_pos.txt'), 'r') as fp:
        numBatches = len(fp.readlines()) -1

    loadFromH5 = False
    TFOrderRandom = False
    if Path(outDir).joinpath('representation_train').exists():
        print("Representation folder for DEEPDRIM training exists, removing representation training folder...")
        shutil.rmtree(Path(outDir).joinpath('representation_train'))
    cmdToRun = ' '.join([
        'singularity exec --no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
        'cd',str(Path('/ext3/data').joinpath(outDir)),';'
        'command time -v -o time_generate_train.txt',
        'python /ext3/DeepDRIM/generate_input_realdata.py ',
        '-out_dir representation_train',
        '-expr_file ',str(inputDir.joinpath('ExpressionData.csv')),
        '-pairs_for_predict_file ',str(inputDir.joinpath('training_pairsDEEPDRIM.txt')),
        '-geneName_map_file ',str(inputDir.joinpath('DEEPDRIM_geneName_map.txt')),
        '-flag_load_from_h5 '+str(loadFromH5),
        '-flag_load_split_batch_pos True',
        '-TF_divide_pos_file '+str(inputDir.joinpath('training_pairsDEEPDRIM.txtTF_divide_pos.txt')),
        '-TF_num '+str(numBatches), '-TF_order_random '+str(TFOrderRandom),
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    if Path(outDir).joinpath('representation_predict').exists():
        print("Representation folder for DEEPDRIM predict exists, removing representation predict folder...")
        shutil.rmtree(Path(outDir).joinpath('representation_predict'))
    cmdToRun = ' '.join([
        'singularity exec --no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
        'cd',str(Path('/ext3/data').joinpath(outDir)),';'
        'command time -v -o time_generate_predict.txt',
        'python /ext3/DeepDRIM/generate_input_realdata.py ',
        '-out_dir representation_predict',
        '-expr_file ',str(inputDir.joinpath('ExpressionData.csv')),
        '-pairs_for_predict_file ',str(inputDir.joinpath('toPredict_pairsDEEPDRIM.txt')),
        '-geneName_map_file ',str(inputDir.joinpath('DEEPDRIM_geneName_map.txt')),
        '-flag_load_from_h5 '+str(loadFromH5),
        '-flag_load_split_batch_pos True',
        '-TF_divide_pos_file '+str(inputDir.joinpath('toPredict_pairsDEEPDRIM.txtTF_divide_pos.txt')),
        '-TF_num '+str(numBatches), '-TF_order_random '+str(TFOrderRandom),
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)
    
    cmdToRun = ' '.join([
        'singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
        '--no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
        'cd',str(Path('/ext3/data').joinpath(outDir)),';',
        'command time -v -o time_train.txt',
        'python /ext3/DeepDRIM/DeepDRIM.py',
        '-num_batches '+str(numBatches),
        '-data_path representation_train/version11/',
        '-output_dir .',
        '-cross_validation_fold_divide_file '+str(inputDir.joinpath('training_pairsDEEPDRIM.txtCV_fold_divide.txt')),
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    with open(RunnerObj.inputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide.txt'),'r') as ifh:
        cross_fold = []
        for line in ifh:
            line = line.strip()
            indel_list = [int(i) for i in line.split(',')]
            cross_fold.append(indel_list)

    for indel_list in cross_fold:
        indel_str = '-'.join([str(i) for i in indel_list])
        model_file = os.path.join('test_'+indel_str+'_saved_models200',
                                  'keras_cnn_trained_model_DeepDRIM.h5')
        cmdToRun = ' '.join([
            'singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
            '--no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
            'cd',str(Path('/ext3/data').joinpath(outDir)),';',
            'command time -v -o time_predict.txt',
            'python /ext3/DeepDRIM/DeepDRIM.py',
            '-to_predict True',
            '-num_batches '+str(numBatches),
            '-data_path representation_predict/version11/',
            '-output_dir '+'test_'+indel_str+'_saved_models200_predict/',
            '-weight_path '+model_file,
        '\"'])
        print(cmdToRun)
        os.system(cmdToRun)

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from DEEPDRIM.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPDRIM/"
    #if not Path(outDir+'outFile.txt').exists():
    #    print(outDir+'outFile.txt'+'does not exist, skipping...')
    #    return
    
    geneNameMap = pd.read_csv(RunnerObj.inputDir.joinpath('DEEPDRIM/GeneName_map.txt'),
                              sep='\t',header=None,index_col=0,names=['name_to'])
    with open(RunnerObj.inputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide.txt'),'r') as ifh:
        cross_fold = []
        for line in ifh:
            line = line.strip()
            indel_list = [int(i) for i in line.split(',')]
            cross_fold.append(indel_list)
            
    outDF_list = []
    for indel_list in cross_fold:
        indel_str = '-'.join([str(i) for i in indel_list])
        predict_dir = 'test_'+indel_str+'_saved_models200_predict'
        for indel in indel_list:
            #print(Path(outDir).joinpath(predict_dir))
            y_predict = pd.read_csv(Path(outDir).joinpath(predict_dir,str(indel)+'end_y_predict.csv'),
                                    sep=',',header=None,index_col=None)
            z_test = pd.read_csv(Path(outDir).joinpath(predict_dir,str(indel)+'end_z_test.csv'),
                                 sep=',',header=0,names=['pair_str'],index_col=0)
            z_test[['TF','target']] = z_test['pair_str'].str.split(',', expand=True)

            outDF_indel=pd.DataFrame({'TF':geneNameMap.loc[z_test['TF'],'name_to'].values,
                                      'target':geneNameMap.loc[z_test['target'],'name_to'].values})
            outDF_indel['importance'] = y_predict
            outDF_list.append(outDF_indel)

    outDF = pd.concat(outDF_list)
    outDF = outDF.sort_values(by=['importance'],ascending=False)
    outDF.to_csv(Path(outDir).joinpath('rankedEdges.csv'),
                 sep='\t',index=False,header=['Gene1','Gene2','EdgeWeight'])

    

import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from sklearn.model_selection import KFold
from sklearn import metrics

add_bind = ''
def generateInputsForPredict(RunnerObj,algName,trainingPairsInputDir):
    if not RunnerObj.inputDir.joinpath(algName).exists():
        print("Input folder for "+algName+" does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath(algName).mkdir(exist_ok = False)

    if not RunnerObj.inputDir.joinpath(algName,"GeneName_map.txt").exists():
        print('Expression data: '+str(RunnerObj.inputDir.joinpath(algName,"ExpressionData.csv")))
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(algName,"ExpressionData.csv"),
                                     header = 0, index_col = 0)
        geneNameMap = pd.concat([ExpressionData.index.to_series().str.lower(),ExpressionData.index.to_series()],axis=1)
        geneNameMap.to_csv(RunnerObj.inputDir.joinpath(algName,"GeneName_map.txt"),
                           sep='\t',header=False,index=False)  
    
def generateInputs(RunnerObj):
    generateInputsForPredict(RunnerObj,'DEEPDRIM',RunnerObj.fullInputDir)
                                          
def generateInputs_old(RunnerObj):
    '''
    Function to generate desired inputs for DeepDRIM.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    if not RunnerObj.inputDir.joinpath("DEEPDRIM").exists():
        print("Input folder for DEEPDRIM does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("DEEPDRIM").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("DEEPDRIM/predict_pairsDEEPDRIM.txtTF_divide_pos.txt").exists():
        # input expression data
        print('Expression data: '+str(RunnerObj.inputDir.joinpath(RunnerObj.exprData)))
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        # output expression .csv file (has header)
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)

        geneNameMap = pd.concat([ExpressionData.index.to_series().str.lower(),ExpressionData.index.to_series()],axis=1)
        geneNameMap.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/GeneName_map.txt"),
                           sep='\t',header=False,index=False)  
        # relative to current path; current path is mounted as /ext3/data
        #inputDir = str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + "/DEEPDRIM"
        #inputDir = RunnerObj.inputDir.relative_to(Path.cwd()).joinpath('DEEPDRIM')

        # Create the pairs to predict file; TF->all genes in the filtered expression data;
        # 3-fold CV split by TF, 1/3 train, 1/3 validation and predict on the last 1/3, so only predictions
        #   for the remaining 1/3 TFs will be used to calculate metrics
        trainingPairs = pd.read_csv(RunnerObj.fullInputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txt'),
                                    sep='\t',header=None,index_col=None)
        trainingDividePos = np.genfromtxt(RunnerObj.fullInputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtTF_divide_pos.txt'),
                                    delimiter=',',dtype='int')
        trainingEdgesFrom = trainingPairs.iloc[trainingDividePos[:-1],0]
        exprGenes = ExpressionData.index.to_series().str.lower()
        predictEdgesFrom = np.intersect1d(trainingEdgesFrom.values,exprGenes.values)
        predictPairs = pd.DataFrame(np.array(np.meshgrid(predictEdgesFrom,exprGenes, indexing='xy')).T.reshape(-1,2),
                                      columns=['Gene1','Gene2'])
        predictPairs['value'] = 3
        predictPairs.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/predict_pairsDEEPDRIM.txt"),
                                      sep='\t',header=False,index=False)

        predictPairs['rnum'] = predictPairs.index
        predictDividePos = [min(predictPairs.loc[predictPairs.Gene1==tf,'rnum']) for tf in predictEdgesFrom]
        predictDividePos.append(predictPairs.shape[0])
        np.savetxt(RunnerObj.inputDir.joinpath("DEEPDRIM/predict_pairsDEEPDRIM.txtTF_divide_pos.txt"),predictDividePos,fmt='%d')
        
def runForPredict(RunnerObj,algName,trainingPairsInputDir):
    '''
    Function to run DEEPDRIM algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
    inputDir = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),algName)
    trainingPairsInputDir = Path('/ext3/data').joinpath(trainingPairsInputDir.relative_to(Path.cwd()),algName)

    # outputs/L0/hESC/DEEPDRIM
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"
    os.makedirs(outDir, exist_ok = True)
    
    with open(RunnerObj.inputDir.joinpath(algName,'predict_pairs'+algName+'.txtTF_divide_pos.txt'), 'r') as fp:
        numBatches = len(fp.readlines()) -1

    loadFromH5 = False
    TFOrderRandom = False

    if Path(outDir).joinpath('representation_predict').exists():
        print("Representation folder for "+algName+" predict exists, removing representation predict folder...")
        shutil.rmtree(Path(outDir).joinpath('representation_predict'))
    cmdToRun = ' '.join([
        'singularity exec --no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
        'cd',str(Path('/ext3/data').joinpath(outDir)),';'
        'command time -v -o time_generate_predict.txt',
        'python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_input_realdata.py ',
        '-out_dir representation_predict',
        '-expr_file ',str(trainingPairsInputDir.joinpath('ExpressionData.csv')),
        '-pairs_for_predict_file ',str(inputDir.joinpath('predict_pairs'+algName+'.txt')),
        '-geneName_map_file ',str(trainingPairsInputDir.joinpath(algName+'_geneName_map.txt')),
        '-flag_load_from_h5 '+str(loadFromH5),
        '-flag_load_split_batch_pos True',
        '-TF_divide_pos_file '+str(inputDir.joinpath('predict_pairs'+algName+'.txtTF_divide_pos.txt')),
        '-TF_num '+str(numBatches), '-TF_order_random '+str(TFOrderRandom),
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)
    
    cv = 3
    for i in range(cv):
        model_file = trainingPairsInputDir.joinpath('training_pairs'+algName+'.txtCV_fold_divide','test_fold-'+str(i)+'_saved_models200','keras_cnn_trained_model_DeepDRIM.h5')
        cmdToRun = ' '.join([
            'singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
            '--no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
            'cd',str(Path('/ext3/data').joinpath(outDir)),';',
            'mkdir -p '+'training_pairs'+algName+'.txtCV_fold_divide/test_fold-'+str(i)+'_saved_models200_predict/; ',
            'command time -v -o time_predict.txt',
            'python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/DeepDRIM.py',
            '-to_predict 2',
            '-num_batches '+str(numBatches),
            '-data_path representation_predict/version11/',
            '-output_dir '+'training_pairs'+algName+'.txtCV_fold_divide/test_fold-'+str(i)+'_saved_models200_predict/',
            '-weight_path '+str(model_file),
        '\"'])
        print(cmdToRun)
        os.system(cmdToRun)
                                          
def run_old(RunnerObj):
    '''
    Function to run DEEPDRIM algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
    inputDir = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),'DEEPDRIM')
    fullInputDir = Path('/ext3/data').joinpath(RunnerObj.fullInputDir.relative_to(Path.cwd()),'DEEPDRIM')

    # outputs/L0/hESC/DEEPDRIM
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPDRIM/"
    os.makedirs(outDir, exist_ok = True)
    
    with open(RunnerObj.inputDir.joinpath('DEEPDRIM/predict_pairsDEEPDRIM.txtTF_divide_pos.txt'), 'r') as fp:
        numBatches = len(fp.readlines()) -1

    loadFromH5 = False
    TFOrderRandom = False

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
        '-expr_file ',str(fullInputDir.joinpath('ExpressionData.csv')),
        '-pairs_for_predict_file ',str(inputDir.joinpath('predict_pairsDEEPDRIM.txt')),
        '-geneName_map_file ',str(fullInputDir.joinpath('DEEPDRIM_geneName_map.txt')),
        '-flag_load_from_h5 '+str(loadFromH5),
        '-flag_load_split_batch_pos True',
        '-TF_divide_pos_file '+str(inputDir.joinpath('predict_pairsDEEPDRIM.txtTF_divide_pos.txt')),
        '-TF_num '+str(numBatches), '-TF_order_random '+str(TFOrderRandom),
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)
    
    cv = 3
    for i in range(cv):
        model_file = fullInputDir.joinpath('training_pairsDEEPDRIM.txtCV_fold_divide','test_fold-'+str(i)+'_saved_models200','keras_cnn_trained_model_DeepDRIM.h5')
        cmdToRun = ' '.join([
            'singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
            '--no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
            'cd',str(Path('/ext3/data').joinpath(outDir)),';',
            'mkdir -p '+'training_pairsDEEPDRIM.txtCV_fold_divide/test_fold-'+str(i)+'_saved_models200_predict/; ',
            'command time -v -o time_predict.txt',
            'python /ext3/DeepDRIM/DeepDRIM.py',
            '-to_predict True',
            '-num_batches '+str(numBatches),
            '-data_path representation_predict/version11/',
            '-output_dir '+'training_pairsDEEPDRIM.txtCV_fold_divide/test_fold-'+str(i)+'_saved_models200_predict/',
            '-weight_path '+str(model_file),
        '\"'])
        print(cmdToRun)
        os.system(cmdToRun)

def runForTrainAndPredict(RunnerObj,algName,trainingPairsInputDir):
    '''
    Function to run DEEPDRIM algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
    inputDirB = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),algName)
    trainingPairsInputDir = trainingPairsInputDir.relative_to(Path.cwd()).joinpath(algName)
    trainingPairsInputDirB = Path('/ext3/data').joinpath(trainingPairsInputDir)
    train_pairs_cv_fold_divide = 'training_pairs'+algName+'.txtCV_fold_divide'

    # outputs/L0/hESC/DEEPDRIM
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"
    outDirB = Path('/ext3/data').joinpath(outDir)
    os.makedirs(outDir, exist_ok = True)

    # Training
    with open(RunnerObj.inputDir.joinpath(algName,'training_pairs'+algName+'.txtTF_divide_pos.txt'), 'r') as fp:
        train_num_batches = len(fp.readlines()) -1
    if trainingPairsInputDir.joinpath(train_pairs_cv_fold_divide).exists():
        print("CV training folder for "+algName+" exists, removing CV training folder...")
        shutil.rmtree(trainingPairsInputDir.joinpath(train_pairs_cv_fold_divide))
        
    cmdToRun = ' '.join([
        'singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
        '--no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
        'cd',str(trainingPairsInputDirB),';'
        'command time -v -o time_train.txt',
        'python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/DeepDRIM.py ',
        '-num_batches '+str(train_num_batches),
        '-data_path representation_train/version11/',
        '-output_dir '+train_pairs_cv_fold_divide+'/',
        '-cross_validation_fold_divide_file '+train_pairs_cv_fold_divide+'.txt',
        '\"'])
    print(cmdToRun)
    os.system(cmdToRun)

    # Prediction
    loadFromH5 = False
    TFOrderRandom = False
    with open(trainingPairsInputDir.joinpath(train_pairs_cv_fold_divide+'.txt'), 'r') as fp:
        cv = len(fp.readlines()) 
    with open(RunnerObj.inputDir.joinpath(algName,'predict_pairs'+algName+'.txtTF_divide_pos.txt'), 'r') as fp:
        predict_num_batches = len(fp.readlines()) -1
    for i in range(cv):
        model_file = trainingPairsInputDirB.joinpath(train_pairs_cv_fold_divide,'test_fold-'+str(i)+'_saved_models200','keras_cnn_trained_model_DeepDRIM.h5')
        cmdToRun = ' '.join([
            'TF_FORCE_GPU_ALLOW_GROWTH=true singularity exec ',RunnerObj.singularityGPUFlag if RunnerObj.params['useGPU'] else '',
            '--no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; ',
            'cd',str(outDirB),';',
            'mkdir -p '+'training_pairs'+algName+'.txtCV_fold_divide/test_fold-'+str(i)+'_saved_models200_predict/; ',
            'command time -v -o time_predict'+str(i)+'.txt',
            'python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/DeepDRIM.py',
            '-to_predict 2',
            '-num_batches '+str(predict_num_batches),
            '-data_path '+str(trainingPairsInputDirB.joinpath('representation_predict/version11'))+'/',
            '-output_dir '+str(Path(train_pairs_cv_fold_divide).joinpath('test_fold-'+str(i)+'_saved_models200_predict'))+'/',
            '-weight_path '+str(model_file),
        '\"'])
        print(cmdToRun)
        os.system(cmdToRun)

def run(RunnerObj):
    runForTrainAndPredict(RunnerObj,'DEEPDRIM',RunnerObj.fullInputDir)

def parseOutputForPredict(RunnerObj,algName,trainingPairsInputDir):
    '''
    Function to parse outputs from DEEPDRIM. 
    For each edge to predict, select the model that was trained on the 
    training fold that contained the source node of this edge in the test set.
   
    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"
    #if not Path(outDir+'outFile.txt').exists():
    #    print(outDir+'outFile.txt'+'does not exist, skipping...')
    #    return
    
    trainingPairs = pd.read_csv(trainingPairsInputDir.joinpath(algName,'training_pairs'+algName+'.txt'),
                                sep='\t',header=None,index_col=None)
    trainingDividePos = np.genfromtxt(trainingPairsInputDir.joinpath(algName,'training_pairs'+algName+'.txtTF_divide_pos.txt'),
                                delimiter=',',dtype='int')
    trainingEdgesFrom = trainingPairs.iloc[trainingDividePos[:-1],0].reset_index(drop=True)

    predictPairs = pd.read_csv(RunnerObj.inputDir.joinpath(algName,'predict_pairs'+algName+'.txt'),
                                sep='\t',header=None,index_col=None)
    predictDividePos = np.genfromtxt(RunnerObj.inputDir.joinpath(algName,'predict_pairs'+algName+'.txtTF_divide_pos.txt'),
                                delimiter=',',dtype='int')
    predictEdgesFrom = predictPairs.iloc[predictDividePos[:-1],0].reset_index(drop=True)
    predictEdgesFromInTrain = trainingEdgesFrom.index[trainingEdgesFrom.isin(predictEdgesFrom)]

    geneNameMap = pd.read_csv(RunnerObj.inputDir.joinpath(algName,'GeneName_map.txt'),
                              sep='\t',header=None,index_col=0,names=['name_to'])
    with open(trainingPairsInputDir.joinpath(algName,'training_pairs'+algName+'.txtCV_fold_divide.txt'),'r') as ifh:
        cross_fold = []
        for line in ifh:
            line = line.strip()
            indel_list = [int(i) for i in line.split(',')]
            cross_fold.append(indel_list)

    outDF_list = []
    for indel_i,indel_list in enumerate(cross_fold):
        predict_dir = 'training_pairs'+algName+'.txtCV_fold_divide/test_fold-'+str(indel_i)+'_saved_models200_predict'
        for indel in indel_list:
            if indel in predictEdgesFromInTrain:# this TF in the predicted pair was in this CV test set in training
                #print(Path(outDir).joinpath(predict_dir))
                tf = trainingEdgesFrom.iloc[indel]
                tfInPredict = predictEdgesFrom.index[predictEdgesFrom==tf][0]
                y_predict = pd.read_csv(Path(outDir).joinpath(predict_dir,str(tfInPredict)+'end_y_predict.csv'),
                                        sep=',',header=None,index_col=None)
                z_test = pd.read_csv(Path(outDir).joinpath(predict_dir,str(tfInPredict)+'end_z_test.csv'),
                                     sep=',',header=0,names=['pair_str'],index_col=0)
                z_test[['TF','target']] = z_test['pair_str'].str.split(',', expand=True)
                
                outDF_indel=pd.DataFrame({'TF':geneNameMap.loc[z_test['TF'],'name_to'].values,
                                          'target':geneNameMap.loc[z_test['target'],'name_to'].values})
                outDF_indel['importance'] = y_predict
                outDF_indel['predict_dir'] = predict_dir
                outDF_list.append(outDF_indel)

    outDF = pd.concat(outDF_list)
    outDF = outDF.sort_values(by=['importance'],ascending=False)
    outDF.to_csv(Path(outDir).joinpath('rankedEdges.csv'),
                 sep='\t',index=False,header=['Gene1','Gene2','EdgeWeight','PredictDir'])


def parseOutputForPredictFromBestModel(RunnerObj,algName,trainingPairsInputDir,metricName='auprc'):
    '''
    Function to parse outputs from DEEPDRIM. 
    Select the model that has the best performance metric out of the CV models, 
    and use this model to calclate predicted scores for all edges.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"
    #if not Path(outDir+'outFile.txt').exists():
    #    print(outDir+'outFile.txt'+'does not exist, skipping...')
    #    return
    geneNameMap = pd.read_csv(RunnerObj.inputDir.joinpath(algName,'GeneName_map.txt'),
                              sep='\t',header=None,index_col=0,names=['name_to'])
    with open(trainingPairsInputDir.joinpath(algName,'training_pairs'+algName+'.txtCV_fold_divide.txt'),'r') as ifh:
        cross_fold = []
        for line in ifh:
            line = line.strip()
            indel_list = [int(i) for i in line.split(',')]
            cross_fold.append(indel_list)

    best_metric = -1.0
    best_i = -1
    # use performance in test set from training
    for indel_i,indel_list in enumerate(cross_fold):
        indel_i_dir = 'training_pairs'+algName+'.txtCV_fold_divide/test_fold-'+str(indel_i)+'_saved_models200'
        y_predict_train = pd.read_csv(trainingPairsInputDir.joinpath(algName,indel_i_dir,'end_y_predict.csv'),
                                      sep=',',header=None,index_col=None)
        y_test_train = pd.read_csv(trainingPairsInputDir.joinpath(algName,indel_i_dir,'end_y_test.csv'),
                                   sep=',',header=None,index_col=None)
        print('indel_list'+str(indel_i),indel_list)
        if (metricName=='auprc'):
            precision, recall, thresholds = metrics.precision_recall_curve(y_test_train,y_predict_train,pos_label=1)
            auprc = metrics.auc(recall,precision)
            print(metricName+':',auprc)
            if (auprc > best_metric):
                best_metric = auprc
                best_i = indel_i
        else:
            print('Metric',metricName,'not recognized. Exit')
            return
    if (best_i==-1):
        print('Unable to find CV model with best metric. Exit')
        return

    # Just read all the predictions from this model
    # Last item in predictDividePos is the position of last line so one more than the number of items
    predictDividePos = np.genfromtxt(RunnerObj.inputDir.joinpath(algName,'predict_pairs'+algName+'.txtTF_divide_pos.txt'),
                                     delimiter=',',dtype='int')
    predict_dir = 'training_pairs'+algName+'.txtCV_fold_divide/test_fold-'+str(best_i)+'_saved_models200_predict'
    print('Selected model:',predict_dir)
    outDF_list = []
    for gene1_i in range(len(predictDividePos)-1):
        y_predict = pd.read_csv(Path(outDir).joinpath(predict_dir,str(gene1_i)+'end_y_predict.csv'),
                                sep=',',header=None,index_col=None)
        z_test = pd.read_csv(Path(outDir).joinpath(predict_dir,str(gene1_i)+'end_z_test.csv'),
                             sep=',',header=0,names=['pair_str'],index_col=0)
        z_test[['TF','target']] = z_test['pair_str'].str.split(',', expand=True)
        
        outDF_i = pd.DataFrame({'TF':geneNameMap.loc[z_test['TF'],'name_to'].values,
                                'target':geneNameMap.loc[z_test['target'],'name_to'].values})
        outDF_i['importance'] = y_predict
        outDF_i['predict_dir'] = predict_dir
        outDF_list.append(outDF_i)

    outDF = pd.concat(outDF_list)
    outDF = outDF.sort_values(by=['importance'],ascending=False)
    outDF.to_csv(Path(outDir).joinpath('rankedEdges.csv'),
                 sep='\t',index=False,header=['Gene1','Gene2','EdgeWeight','PredictDir'])
    
def parseOutput(RunnerObj):
    parseOutputForPredict(RunnerObj,'DEEPDRIM',RunnerObj.fullInputDir)
                                          
def parseOutput_old(RunnerObj):
    '''
    Function to parse outputs from DEEPDRIM.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPDRIM/"
    #if not Path(outDir+'outFile.txt').exists():
    #    print(outDir+'outFile.txt'+'does not exist, skipping...')
    #    return
    
    trainingPairs = pd.read_csv(RunnerObj.fullInputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txt'),
                                sep='\t',header=None,index_col=None)
    trainingDividePos = np.genfromtxt(RunnerObj.fullInputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtTF_divide_pos.txt'),
                                delimiter=',',dtype='int')
    trainingEdgesFrom = trainingPairs.iloc[trainingDividePos[:-1],0].reset_index(drop=True)

    predictPairs = pd.read_csv(RunnerObj.inputDir.joinpath('DEEPDRIM/predict_pairsDEEPDRIM.txt'),
                                sep='\t',header=None,index_col=None)
    predictDividePos = np.genfromtxt(RunnerObj.inputDir.joinpath('DEEPDRIM/predict_pairsDEEPDRIM.txtTF_divide_pos.txt'),
                                delimiter=',',dtype='int')
    predictEdgesFrom = predictPairs.iloc[predictDividePos[:-1],0].reset_index(drop=True)
    predictEdgesFromInTrain = trainingEdgesFrom.index[trainingEdgesFrom.isin(predictEdgesFrom)]
    
    geneNameMap = pd.read_csv(RunnerObj.inputDir.joinpath('DEEPDRIM/GeneName_map.txt'),
                              sep='\t',header=None,index_col=0,names=['name_to'])
    with open(RunnerObj.fullInputDir.joinpath('DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide.txt'),'r') as ifh:
        cross_fold = []
        for line in ifh:
            line = line.strip()
            indel_list = [int(i) for i in line.split(',')]
            cross_fold.append(indel_list)

    outDF_list = []
    for indel_i,indel_list in enumerate(cross_fold):
        predict_dir = 'training_pairsDEEPDRIM.txtCV_fold_divide/test_fold-'+str(indel_i)+'_saved_models200_predict'
        for indel in indel_list:
            if indel in predictEdgesFromInTrain:# this TF in the predicted pair was in the test set in training
                #print(Path(outDir).joinpath(predict_dir))
                tf = trainingEdgesFrom.iloc[indel]
                tfInPredict = predictEdgesFrom.index[predictEdgesFrom==tf][0]
                y_predict = pd.read_csv(Path(outDir).joinpath(predict_dir,str(tfInPredict)+'end_y_predict.csv'),
                                        sep=',',header=None,index_col=None)
                z_test = pd.read_csv(Path(outDir).joinpath(predict_dir,str(tfInPredict)+'end_z_test.csv'),
                                     sep=',',header=0,names=['pair_str'],index_col=0)
                z_test[['TF','target']] = z_test['pair_str'].str.split(',', expand=True)
                
                outDF_indel=pd.DataFrame({'TF':geneNameMap.loc[z_test['TF'],'name_to'].values,
                                          'target':geneNameMap.loc[z_test['target'],'name_to'].values})
                outDF_indel['importance'] = y_predict
                outDF_indel['predict_dir'] = predict_dir
                outDF_list.append(outDF_indel)

    outDF = pd.concat(outDF_list)
    outDF = outDF.sort_values(by=['importance'],ascending=False)
    outDF.to_csv(Path(outDir).joinpath('rankedEdges.csv'),
                 sep='\t',index=False,header=['Gene1','Gene2','EdgeWeight','PredictDir'])

    

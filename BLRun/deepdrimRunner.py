import os
import pandas as pd
from pathlib import Path
import numpy as np

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
        
    if not RunnerObj.inputDir.joinpath("DEEPDRIM/ExpressionData.csv").exists():
        # input expression data
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        newExpressionData = ExpressionData.copy()
        
        # output expression .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)
        
        # Read file for trueEdges and create training pairs and TF divide pos file
        print('True edges: '+str(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges)))
        trueEdgesDF = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.trueEdges),
                                sep = ',', header = 0, index_col = None)
        trueEdgesDFSorted = trueEdges.sort_values(['Gene1']).reset_index(drop=True)
        # two columns: Gene1, position of Gene 1 group in trueEdgesDFSorted
        sourcePos = trueEdgesDFSorted.index.to_series().groupby(x['trueEdgesDFSorted']).first().reset_index(name='pos')

        trueEdgesDFSorted.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/TrainingPairs.txt"),
                                 sep='\t',header=False,index=False)
        sourcePos['pos'].to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/TrainingPairsDivPos.txt"),
                                sep='\t',header=False,index=False)

        # Gene name map file
        geneNameMap = pd.concat([ExpressionData.index.to_series(),ExpressionData.index.to_series()],axis=1)
        geneNameMap.to_csv(RunnerObj.inputDir.joinpath("DEEPDRIM/GeneNameMap.txt"),
                           sep='\t',header=False,index=False)

        numBatches = RunnerObj.params['num_batches']
        loadFromH5 = 'False'
        loadSplitBatchPos = 'True'
        TFOrderRandom = 'True'

        cmdToRun = ' '.join([
            'singularity exec --no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/',
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/sh -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd /ext3 ; ',
            'python python/DeepDRIM/generate_input_realdata.py ',
            '-out_dir '+RunnerObj.inputDir.joinpath("DEEPDRIM/representation"),
            '-expr_file '+RunnerObj.inputDir.joinpath("DEEPDRIM/ExpressionData.csv"),
            '-pairs_for_predict_file '+RunnerObj.inputDir.joinpath("DEEPDRIM/TrainingPairs.txt"),
            '-geneName_map_file '+RunnerObj.inputDir.joinpath("DEEPDRIM/GeneNameMap.txt"),
            '-flag_load_from_h5',loadFromH5,'-flag_load_split_batch_pos',loadSplitBatchPos,
            '-TF_divide_pos_file '+RunnerObj.inputDir.joinpath("DEEPDRIM/TrainingPairsDivPos.txt"),
            '-TF_num',numBatches, '-TF_order_random',TFOrderRandom,
            '\"'])
    print(cmdToRun)
    os.system(cmdToRun)
    
def run(RunnerObj):
    '''
    Function to run GENIE3 algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    inputDir = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/DEEPDRIM"
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPDRIM/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "data/" +  str(outDir) + 'outFile.txt'
    #cmdToRun = ' '.join(['docker run --rm -v',
    #                     str(Path.cwd())+':/data/ --expose=41269',
    #                     'arboreto:base /bin/sh -c \"time -v -o',
    #                     "data/" + str(outDir) + 'time.txt', 'python runArboreto.py --algo=GENIE3',
    #                     '--inFile='+inputPath, '--outFile='+outPath, '\"'])

    cmdToRun = ' '.join([
            'singularity exec --no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/',
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/sh -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd /ext3 ; ',
            'time -v -o',"data/" + str(outDir) + 'time.txt',
            'python python/DeepDRIM/DeepDRIM.py',
            '-num_batches',RunnerObj.params['num_batches'],
            '-data_path',inputDir.joinpath('representation/version11'),
            '-output_dir',outDir,
            'cross_validation_fold_divide_file',RunnerObj['cross_validation_fold_divide_file'],
            '\"'])

    print(cmdToRun)
    os.system(cmdToRun)

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GENIE3.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/DEEPDRIM/"

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
    

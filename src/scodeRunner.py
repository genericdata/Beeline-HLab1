import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCODE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SCODE").exists():
        print("Input folder for SCODE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCODE").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("SCODE/ExpressionData.csv").exists():
        # input file
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        # output file
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("SCODE/ExpressionData.csv"),
                             sep = '\t', header  = False, index = False)
        
    if not RunnerObj.inputDir.joinpath("SCODE/PseudoTime.csv").exists():
        # input file
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
        # output file
        PTData.to_csv(RunnerObj.inputDir.joinpath("SCODE/PseudoTime.csv"),
                             sep = '\t', header  = False)
        
    
def run(RunnerObj):
    '''
    Function to run SCODE algorithm
    '''
    inputPath = "data/"+str(RunnerObj.inputDir).split("RNMethods/")[1]+"/SCODE/"
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    nGenes = str(ExpressionData.shape[0])
    z = str(RunnerObj.params['z'])
    nCells = str(ExpressionData.shape[1])
    nIter = str(RunnerObj.params['nIter'])
    nRep = str(RunnerObj.params['nRep'])
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCODE/"
    os.makedirs(outDir, exist_ok = True)
    
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/SCODE/data/  scode:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 'ruby run_R.rb',
                    inputPath +'ExpressionData.csv', inputPath + 'PseudoTime.csv', 
                         ' data/'+outDir,
                         nGenes, z, nCells, nIter, nRep, '\"'])
    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCODE.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCODE/"
    if not Path(outDir+'meanA.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'meanA.txt', sep = '\t', header = None)
    
    # Sort values in a matrix using code from:
    # https://stackoverflow.com/questions/21922806/sort-values-of-matrix-in-python
    OutMatrix = np.abs(OutDF.values)
    idx = np.argsort(OutMatrix, axis = None)[::-1]
    rows, cols = np.unravel_index(idx, OutDF.shape)    
    DFSorted = OutMatrix[rows, cols]
    
    # read input file for list of gene names
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    GeneList = list(ExpressionData.index)
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for row, col, val in zip(rows, cols, DFSorted):
        outFile.write('\t'.join([GeneList[row],GeneList[col],str(val)])+'\n')
    outFile.close()
    
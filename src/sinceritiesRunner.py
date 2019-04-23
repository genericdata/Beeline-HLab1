import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SINCERITIES.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SINCERITIES").exists():
        print("Input folder for SINCERITIES does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SINCERITIES").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("SINCERITIES/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        newExpressionData = ExpressionData.T.copy()
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
        # make sure the indices are strings for both dataframes
        newExpressionData.index = newExpressionData.index.map(str) 
        PTData.index = PTData.index.map(str) 
        newExpressionData['Time'] = PTData['Time']
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("SINCERITIES/ExpressionData.csv"),
                             sep = ',', header  = True, index = False)
    
def run(RunnerObj):
    '''
    Function to run SINCERITIES algorithm
    '''
    inputPath = "data/" + str(RunnerObj.inputDir).split("RNMethods/")[1] + \
                    "/SINCERITIES/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SINCERITIES/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "data/" + str(outDir) + 'outFile.txt'
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/SINCERITIES/data/ sincerities:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 'Rscript MAIN.R',
                         inputPath, outPath, '\"'])
    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SINCERITIES.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SINCERITIES/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = ',', header = 0)
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.iterrows():
        outFile.write('\t'.join([row['SourceGENES'],row['TargetGENES'],str(row['Interaction'])])+'\n')
    outFile.close()
    

import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for CICT.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("CICT").exists():
        print("Input folder for CICT does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("CICT").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("CICT/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        
        newExpressionData = ExpressionData.copy()
        
        # Write .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("CICT/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)
    
def run(RunnerObj):
    '''
    Function to run CICT algorithm
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/CICT/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/CICT/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "data/" +  str(outDir) + 'outFile.txt'
    #cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data/ CICT:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', 'Rscript runCICT.R',
    #                     inputPath, outPath, '\"'])

    cmdToRun = ' '.join([
        'singularity exec --no-home',
        '-B ' + str(Path.cwd())+':/ext3/data/',
        '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
        str(RunnerObj.singularityImage),
        '/bin/sh -c \"source /ext3/env.sh; conda activate CICT; cd /ext3 ; time -v -o', "data/" + str(outDir) + 'time.txt', 'Rscript runCICT.R',
        inputPath, outPath, '\"'])


    print(cmdToRun)
    os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from CICT.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/CICT/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = 0)
    # edges with significant p-value
    part1 = OutDF.loc[OutDF['pValue'] <= float(RunnerObj.params['pVal'])]
    part1 = part1.assign(absCorVal = part1['corVal'].abs())
    # edges without significant p-value
    part2 = OutDF.loc[OutDF['pValue'] > float(RunnerObj.params['pVal'])]
    
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in part1.sort_values('absCorVal', ascending = False).iterrows():
        outFile.write('\t'.join([row['Gene1'],row['Gene2'],str(row['corVal'])])+'\n')
    
    for idx, row in part2.iterrows():
        outFile.write('\t'.join([row['Gene1'],row['Gene2'],str(0)])+'\n')
    outFile.close()
    

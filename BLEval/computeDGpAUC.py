import os
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc, average_precision_score
from itertools import product, permutations, combinations, combinations_with_replacement
from operator import itemgetter
from tqdm import tqdm
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector,pandas2ri

def computeScoresByAdjMat(uniqNodes,trueEdgesDF,predEdgeDF,directed,selfEdges):
                                 
    if directed:        
        # initialize adjacency matrix and dictionary of node name to index
        nodeIndexDict = {val:idx for idx, val in enumerate(uniqNodes)}# gene name->index
        numNodes = len(uniqNodes)
        TrueEdgeAdj = np.zeros([numNodes,numNodes],dtype=int)
        PredEdgeAdj = np.zeros([numNodes,numNodes],dtype=float)

        # Compute TrueEdgeAdj adjacency matrix
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        trueEdgesSelect = (trueEdgesDF['Gene1'].isin(nodeIndexDict)) & (trueEdgesDF['Gene2'].isin(nodeIndexDict))
        trueEdgesFrom = trueEdgesDF.loc[trueEdgesSelect,'Gene1']
        trueEdgesFromIdx = itemgetter(*trueEdgesFrom)(nodeIndexDict)
        trueEdgesTo = trueEdgesDF.loc[trueEdgesSelect,'Gene2']
        trueEdgesToIdx = itemgetter(*trueEdgesTo)(nodeIndexDict)
        TrueEdgeAdj[(trueEdgesFromIdx,trueEdgesToIdx)] = 1

        # Get the predicted edge scores
        predEdgeDF['EdgeWeightAbs'] = predEdgeDF['EdgeWeight'].abs()
        predEdgeMax = predEdgeDF.loc[(predEdgeDF['Gene1'].isin(nodeIndexDict)) &
                                     (predEdgeDF['Gene2'].isin(nodeIndexDict))]\
                                     .groupby(['Gene1','Gene2'])['EdgeWeightAbs'].first().reset_index()
        #predEdgeMax = predEdgeDF.groupby(['Gene1','Gene2'])['EdgeWeightAbs'].first().reset_index() 
        predEdgeFrom = predEdgeMax.loc[:,'Gene1']
        predEdgeFromIdx = itemgetter(*predEdgeFrom)(nodeIndexDict)
        predEdgeTo = predEdgeMax.loc[:,'Gene2']
        predEdgeToIdx = itemgetter(*predEdgeTo)(nodeIndexDict)
        PredEdgeAdj[(predEdgeFromIdx,predEdgeToIdx)] = predEdgeMax.loc[:,'EdgeWeightAbs']

    # if undirected
    else:
        # initialize adjacency matrix and dictionary of node name to index
        nodeIndexDict = {val:idx for idx, val in enumerate(uniqNodes)}# gene name->index
        numNodes = len(uniqNodes)
        TrueEdgeAdj = np.zeros([numNodes,numNodes],dtype=int)
        PredEdgeAdjF = np.zeros([numNodes,numNodes],dtype=float)
        PredEdgeAdjR = np.zeros([numNodes,numNodes],dtype=float)

        # Compute TrueEdgeAdj adjacency matrix
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        trueEdgesSelect = (trueEdgesDF['Gene1'].isin(nodeIndexDict)) & (trueEdgesDF['Gene2'].isin(nodeIndexDict))
        trueEdgesFrom = trueEdgesDF.loc[trueEdgesSelect,'Gene1']
        trueEdgesFromIdx = itemgetter(*trueEdgesFrom)(nodeIndexDict)
        trueEdgesTo = trueEdgesDF.loc[trueEdgesSelect,'Gene2']

        # Get the predicted edge scores
        predEdgeDF['EdgeWeightAbs'] = predEdgeDF['EdgeWeight'].abs()
        predEdgeMax = predEdgeDF.loc[(predEdgeDF['Gene1'].isin(nodeIndexDict)) &
                                     (predEdgeDF['Gene2'].isin(nodeIndexDict))]\
                                     .groupby(['Gene1','Gene2'])['EdgeWeightAbs'].max().reset_index()
        #predEdgeMax = predEdgeDF.groupby(['Gene1','Gene2'])['EdgeWeightAbs'].max().reset_index()
        predEdgeFromF = predEdgeMax.loc[:,'Gene1']
        predEdgeFromIdxF = itemgetter(*predEdgeFromF)(nodeIndexDict)
        predEdgeToF = predEdgeMax.loc[:,'Gene2']
        predEdgeToIdxF = itemgetter(*predEdgeToF)(nodeIndexDict)
        PredEdgeAdjF[(predEdgeFromIdxF,predEdgeToIdxF)] = predEdgeMax.loc[:,'EdgeWeightAbs']
        PredEdgeAdjF[(predEdgeToIdxF,predEdgeFromIdxF)] = predEdgeMax.loc[:,'EdgeWeightAbs']
        predEdgeFromR = predEdgeMax.loc[:,'Gene2']
        predEdgeFromIdxR = itemgetter(*predEdgeFromR)(nodeIndexDict)
        predEdgeToR = predEdgeMax.loc[:,'Gene1']
        predEdgeToIdxR = itemgetter(*predEdgeToR)(nodeIndexDict)
        PredEdgeAdjR[(predEdgeFromIdxR,predEdgeToIdxR)] = predEdgeMax.loc[:,'EdgeWeightAbs']
        PredEdgeAdjR[(predEdgeToIdxR,predEdgeFromIdxR)] = predEdgeMax.loc[:,'EdgeWeightAbs']
        PredEdgeAdj = np.maximum(PredEdgeAdjF,PredEdgeAdjR)

    if selfEdges:
        TrueEdgeArr = TrueEdgeAdj.reshape(-1)
        PredEdgeArr = PredEdgeAdj.reshape(-1)
    else:
        TrueEdgeArr = TrueEdgeAdj[~np.eye(TrueEdgeAdj.shape[0],dtype=bool)]
        PredEdgeArr = PredEdgeAdj[~np.eye(PredEdgeAdj.shape[0],dtype=bool)]

    print('Edge label and score matrices constructed')
    # Combine into one dataframe
    # to pass it to sklearn
    print('Computing PR curve and partial area in R')
    precrec = importr('precrec')
    PredEdgeFV = FloatVector(PredEdgeArr)
    TrueEdgeFV = FloatVector(TrueEdgeArr)
    print('   Convert to R obj done; start calculation')

    ssCurve = precrec.evalmod(scores=PredEdgeFV,labels=TrueEdgeFV)
    ssCurvePart = precrec.part(ssCurve, xlim = FloatVector([0.0, 0.2]))
    ssCurvePAUC = pandas2ri.ri2py(precrec.pauc(ssCurvePart))
    ssCurvePAUC['metric'] = 'AU' + ssCurvePAUC['curvetypes'] # AUPRC or AUROC
    ssCurvePartX = precrec.part(ssCurve, ylim = FloatVector([0.8, 1.0]))
    ssCurvePAUCX = pandas2ri.ri2py(precrec.pauc(ssCurvePartX))
    ssCurvePAUCX['metric'] = 'AU' + ssCurvePAUC['curvetypes'] + 'X' # AUPRCX or AUROCX

    pAUCDF = pd.concat([pd.melt(ssCurvePAUC,id_vars=['modnames','dsids','curvetypes','metric']),
                        pd.melt(ssCurvePAUCX,id_vars=['modnames','dsids','curvetypes','metric'])])
    metric_prefix = pd.Series({'paucs':'p','spaucs':'sp'})
    
    pAUCDF['metric'] = metric_prefix[pAUCDF['variable']].values + pAUCDF['metric'].values
    pAUCDict = dict(zip(pAUCDF['metric'],pAUCDF['value']))
    return pAUCDict
    
    
def pPRROC3(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
    '''
    Computes areas under the precision-recall and ROC curves
    for a given dataset for each algorithm.
    

    :param directed:   A flag to indicate whether to treat predictions as directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
    :param selfEdges:   A flag to indicate whether to includeself-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
    :param plotFlag:   A flag to indicate whether or not to save PR and ROC plots.
    :type plotFlag: bool
        
    :returns:
            - pAUC: A dictionary containing for each algorithm, a dictionary of various partial AUC values
    '''
    
    # Read file for trueEdges
    print('True edges: '+str(inputSettings.datadir)+'/'+ dataDict['name'] +'/' +dataDict['trueEdges'])
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    # Initialize data dictionaries
    pAUC = {}
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    print('Outdir: '+outDir)
    
    if directed:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):
            
            print(algo)
            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion
                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)
                print(trueEdgesDF.head())
                print(predDF.head())
                if predDF.shape[0] > 0:    
                    pAUC[algo[0]] = computeScores3(trueEdgesDF, predDF, directed = True, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')                    
            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/PRplot'
            ROCName = '/ROCplot'
    else:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            print(algo)
            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion
                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)
                print(trueEdgesDF.head())
                print(predDF.head())
                if predDF.shape[0] > 0:
                    pAUC[algo[0]] = computeScores3(trueEdgesDF, predDF, directed = False, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')
            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                  ' does not exist. Skipping...')

    return pAUC


def computeScores3(trueEdgesDF, predEdgeDF, 
                  directed = True, selfEdges = True):
    '''        
    Computes precision-recall and ROC curves
    using scikit-learn for a given set of predictions in the 
    form of a DataFrame.
    
    :param trueEdgesDF:   A pandas dataframe containing the true classes.The indices of this dataframe are all possible edgesin a graph formed using the genes in the given dataset. This dataframe only has one column to indicate the classlabel of an edge. If an edge is present in the reference network, it gets a class label of 1, else 0.
    :type trueEdgesDF: DataFrame
        
    :param predEdgeDF:   A pandas dataframe containing the edge ranks from the prediced network. The indices of this dataframe are all possible edges.This dataframe only has one column to indicate the edge weightsin the predicted network. Higher the weight, higher the edge confidence.
    :type predEdgeDF: DataFrame
    
    :param directed:   A flag to indicate whether to treat predictionsas directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
        
    :param selfEdges:   A flag to indicate whether to includeself-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
        
    :returns:
            - pAUCDict: A dictionary containing for each algorithm, a dictionary of various partial AUC values
    '''
    uniqNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
    return computeScoresByAdjMat(uniqNodes,trueEdgesDF, predEdgeDF,
                                 directed, selfEdges)



def pPRROC4(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
    '''
    Computes areas under the precision-recall and ROC curves
    for a given dataset for each algorithm.
    

    :param directed:   A flag to indicate whether to treat predictions as directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
    :param selfEdges:   A flag to indicate whether to includeself-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
    :param plotFlag:   A flag to indicate whether or not to save PR and ROC plots.
    :type plotFlag: bool
        
    :returns:
            - pAUCDict: A dictionary containing for each algorithm, a DataFrame of various partial AUC values
    '''
    # Read file for expression data
    print('Expression data: '+ str(inputSettings.datadir)+'/'+ dataDict['name'] +'/' +dataDict['exprData'])
    expressionData = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name']
                                 +'/' +dataDict['exprData'],
                                 header = 0, index_col = 0)
    
    # Read file for trueEdges
    print('True edges: '+str(inputSettings.datadir)+'/'+ dataDict['name'] +'/' +dataDict['trueEdges'])
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    # Initialize data dictionaries
    pAUC = {}
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    print('Outdir: '+outDir)
    
    if directed:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):
            
            print(algo)
            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion
                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)
                print(trueEdgesDF.head())
                print(predDF.head())
                if predDF.shape[0] > 0:    
                    pAUC[algo[0]] = computeScores4(expressionData,trueEdgesDF, predDF, directed = True, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')                    
            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/PRplot'
            ROCName = '/ROCplot'
    else:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            print(algo)
            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion
                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)
                print(trueEdgesDF.head())
                print(predDF.head())
                if predDF.shape[0] > 0:
                    pAUC[algo[0]] = computeScores4(expressionData, trueEdgesDF, predDF, directed = False, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')
            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                  ' does not exist. Skipping...')
                
    return pAUC

def computeScores4(expressionData, trueEdgesDF, predEdgeDF, 
                  directed = True, selfEdges = True):
    '''        
    Computes precision-recall and ROC curves
    using scikit-learn for a given set of predictions in the 
    form of a DataFrame.

    :param expressionData:   A pandas dataframe containing the expression data used in the inference.
    :type expressionData: DataFrame

    :param trueEdgesDF:   A pandas dataframe containing the true classes.The indices of this dataframe are all possible edgesin a graph formed using the genes in the given dataset. This dataframe only has one column to indicate the classlabel of an edge. If an edge is present in the reference network, it gets a class label of 1, else 0.
    :type trueEdgesDF: DataFrame
        
    :param predEdgeDF:   A pandas dataframe containing the edge ranks from the prediced network. The indices of this dataframe are all possible edges.This dataframe only has one column to indicate the edge weightsin the predicted network. Higher the weight, higher the edge confidence.
    :type predEdgeDF: DataFrame
    
    :param directed:   A flag to indicate whether to treat predictionsas directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
        
    :param selfEdges:   A flag to indicate whether to includeself-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
        
    :returns:
            - pAUCDict: A dictionary containing for each algorithm, a dictionary of various partial AUC values
    '''
    uniqNodes = np.unique(expressionData.index)
    return computeScoresByAdjMat(uniqNodes,trueEdgesDF, predEdgeDF,
                                 directed, selfEdges)

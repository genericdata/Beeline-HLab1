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
from rpy2.robjects import FloatVector


def PRROC(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
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
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''
    
    # Read file for trueEdges
    print('True edges: '+str(inputSettings.datadir)+'/'+ dataDict['name'] +'/' +dataDict['trueEdges'])
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    print('Outdir: '+outDir)
    
    if directed:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion
                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)

                print(algo)
                print(trueEdgesDF.head())
                print(predDF.head())
                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores(trueEdgesDF, predDF, directed = True, selfEdges = selfEdges)

            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/' + runDir + 'PRplot'
            ROCName = '/' + runDir + 'ROCplot'
    else:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):
            
            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion

                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)

                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores(trueEdgesDF, predDF, directed = False, selfEdges = selfEdges)

            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
            PRName = '/' + runDir + 'uPRplot'
            ROCName = '/' + runDir + 'uROCplot'
    if (plotFlag):
         ## Make PR curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(recallDict[key],precisionDict[key], ci=None)
            legendList.append(key + ' (AUPRC = ' + str("%.2f" % (AUPRC[key]))+')')
        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(legendList) 
        plt.savefig(outDir+PRName+'.pdf')
        plt.savefig(outDir+PRName+'.png')
        plt.clf()

        ## Make ROC curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(FPRDict[key],TPRDict[key], ci=None)
            legendList.append(key + ' (AUROC = ' + str("%.2f" % (AUROC[key]))+')')

        plt.plot([0, 1], [0, 1], linewidth = 1.5, color = 'k', linestyle = '--')

        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.legend(legendList) 
        plt.savefig(outDir+ROCName+'.pdf')
        plt.savefig(outDir+ROCName+'.png')
        plt.clf()
    return AUPRC, AUROC


def computeScores(trueEdgesDF, predEdgeDF, 
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
            - prec: A list of precision values (for PR plot)
            - recall: A list of precision values (for PR plot)
            - fpr: A list of false positive rates (for ROC plot)
            - tpr: A list of true positive rates (for ROC plot)
            - AUPRC: Area under the precision-recall curve
            - AUROC: Area under the ROC curve
    '''

    if directed:        
        # Initialize dictionaries with all 
        # possible edges
        if selfEdges:
            possibleEdges = list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         repeat = 2))
        else:
            possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         r = 2))
        
        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        
        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                   (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
                    TrueEdgeDict[key] = 1
                
        for key in PredEdgeDict.keys():
            subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])]
            if len(subDF)>0:
                PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])

    # if undirected
    else:
        
        # Initialize dictionaries with all 
        # possible edges
        if selfEdges:
            possibleEdges = list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2))
        else:
            possibleEdges = list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2))
        TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth

        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])) |
                              ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene1'] == key.split('|')[1]))]) > 0:
                TrueEdgeDict[key] = 1  

        # Compute PredEdgeDict Dictionary
        # from predEdgeDF

        for key in PredEdgeDict.keys():
            subDF = predEdgeDF.loc[((predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])) |
                              ((predEdgeDF['Gene2'] == key.split('|')[0]) &
                               (predEdgeDF['Gene1'] == key.split('|')[1]))]
            if len(subDF)>0:
                PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))

                
                
    # Combine into one dataframe
    # to pass it to sklearn
    outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    outDF.columns = ['TrueEdges','PredEdges']
    prroc = importr('PRROC')
    prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(outDF['PredEdges'].values)), 
              weights_class0 = FloatVector(list(outDF['TrueEdges'].values)))

    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)
    
    return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr)

def PRROC2(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
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
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''
    
    # Read file for trueEdges
    print('True edges: '+os.path.abspath(str(inputSettings.datadir)+'/'+ dataDict['name'] +'/' +dataDict['trueEdges']))
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' +dataDict['name']
    print('Outdir: '+os.path.abspath(outDir))
    
    if directed:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion
                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)

                print(algo)
                print(trueEdgesDF.head())
                print(predDF.head())
                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores2(trueEdgesDF, predDF, directed = True, selfEdges = selfEdges)

            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/' + runDir + 'PRplot'
            ROCName = '/' + runDir + 'ROCplot'
    else:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            runDir = algo[1]['run_dir']+'/' if 'run_dir' in algo[1] else ''
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion

                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)

                precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]] = computeScores2(trueEdgesDF, predDF, directed = False, selfEdges = selfEdges)

            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
            PRName = '/' + runDir + 'uPRplot'
            ROCName = '/' + runDir + 'uROCplot'
    if (plotFlag):
         ## Make PR curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(recallDict[key],precisionDict[key], ci=None)
            legendList.append(key + ' (AUPRC = ' + str("%.2f" % (AUPRC[key]))+')')
        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(legendList) 
        plt.savefig(outDir+PRName+'.pdf')
        plt.savefig(outDir+PRName+'.png')
        plt.clf()

        ## Make ROC curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(FPRDict[key],TPRDict[key], ci=None)
            legendList.append(key + ' (AUROC = ' + str("%.2f" % (AUROC[key]))+')')

        plt.plot([0, 1], [0, 1], linewidth = 1.5, color = 'k', linestyle = '--')

        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.legend(legendList) 
        plt.savefig(outDir+ROCName+'.pdf')
        plt.savefig(outDir+ROCName+'.png')
        plt.clf()
    return AUPRC, AUROC

def computeScores2(trueEdgesDF, predEdgeDF, 
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
            - prec: A list of precision values (for PR plot)
            - recall: A list of precision values (for PR plot)
            - fpr: A list of false positive rates (for ROC plot)
            - tpr: A list of true positive rates (for ROC plot)
            - AUPRC: Area under the precision-recall curve
            - AUROC: Area under the ROC curve
    '''
    if directed:        
        # Initialize dictionaries with all 
        # possible edges
        if selfEdges:
            #print('self started')
            #possibleEdges = list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                             repeat = 2))
            possibleEdges = pd.DataFrame(list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         repeat = 2)))
            #print('self end')
        else:
            #print('not self started')
            #possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                             r = 2))
            possibleEdges = pd.DataFrame(list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                              r = 2)))
            #print('not self end')

        #print(len(possibleEdges))
        #TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        #PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        edgeDict = dict((k, [0, 0]) for k in possibleEdges[0]+'|'+possibleEdges[1])# first element is true edge, second element is predicted edge
        #print(len(TrueEdgeDict))
        
        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        #for key in TrueEdgeDict.keys():
        #    if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
        #           (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
        #            TrueEdgeDict[key] = 1
        trueEdgeKeys = trueEdgesDF['Gene1'].astype(str) + '|' + trueEdgesDF['Gene2']
        for key in trueEdgeKeys:
            if key in edgeDict:
                edgeDict[key][0] = 1
        #print('done TrueEdgeDict')
        
        #print(predEdgeDF.EdgeWeight.values[0])
        #for key in PredEdgeDict.keys():
        #    subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
        #                       (predEdgeDF['Gene2'] == key.split('|')[1])]
        #    if len(subDF)>0:
        #        PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
        #        print(subDF.EdgeWeight.values)
        predEdgeKeys = predEdgeDF['Gene1'].astype(str) + '|' + predEdgeDF['Gene2']
        predEdgeValues =  np.abs(predEdgeDF.EdgeWeight)
        for k,v in zip(predEdgeKeys,predEdgeValues):
            if k in edgeDict:
                edgeDict[k][1] = v
        #print('done PredEdgeDict')
    # if undirected
    else:
        
        # Initialize dictionaries with all 
        # possible edges
        if selfEdges:
            #possibleEdges = list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                                                   r = 2))
            possibleEdges = pd.DataFrame(list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2)))
        else:
            #possibleEdges = list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                                                   r = 2))
            possibleEdges = pd.DataFrame(list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2)))
        #TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        #PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        edgeDict = dict((k, [0, 0]) for k in possibleEdges[0]+'|'+possibleEdges[1])# first element is true edge, second element is predicted edge

        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth

        #for key in TrueEdgeDict.keys():
        #    if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
        #                   (trueEdgesDF['Gene2'] == key.split('|')[1])) |
        #                      ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
        #                   (trueEdgesDF['Gene1'] == key.split('|')[1]))]) > 0:
        #        TrueEdgeDict[key] = 1  
        trueEdgeKeys = pd.concat([trueEdgesDF['Gene1'].astype(str) + '|' + trueEdgesDF['Gene2'],
                                  trueEdgesDF['Gene2'].astype(str) + '|' + trueEdgesDF['Gene1']],
                                 ignore_index=True)
        for key in trueEdgeKeys:
            if key in edgeDict:
                edgeDict[key][0] = 1
        #print('done TrueEdgeDict')

        # Compute PredEdgeDict Dictionary
        # from predEdgeDF

        #for key in PredEdgeDict.keys():
        #    subDF = predEdgeDF.loc[((predEdgeDF['Gene1'] == key.split('|')[0]) &
        #                       (predEdgeDF['Gene2'] == key.split('|')[1])) |
        #                      ((predEdgeDF['Gene2'] == key.split('|')[0]) &
        #                       (predEdgeDF['Gene1'] == key.split('|')[1]))]
        #    if len(subDF)>0:
        #        PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))
        
        # take max egde weights of two directions
        predEdgeDF_BiDir = pd.concat(predEdgeDF,
                                     predEdgeDF.rename(columns={"Gene1": "Gene2", "Gene2": "Gene1"},
                                                       copy=TRUE),
                                     ignore_index=True)
        predEdgeDF_BiDir['EdgeWeight'] = np.abs(predEdgeDF_BiDir['EdgeWeight'])
        predEdgeDF_max = predEdgeDF_BiDir.groupby(['Gene1','Gene2'])['EdgeWeight'].max().reset_index()
        predEdgeKeys = predEdgeDF_max['Gene1'].astype(str) + '|' + predEdgeDF_max['Gene2']
        predEdgeValues =  predEdgesDF_max.EdgeWeight.values
        for k,v in zip(predEdgeKeys,predEdgeValues):
            if k in edgeDict:
                edgeDict[k][1] = v
        #print('done PredEdgeDict')


    print('Edge label and score matrices constructed')
    # Combine into one dataframe
    # to pass it to sklearn
    #outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    #outDF.columns = ['TrueEdges','PredEdges']
    outDF = pd.DataFrame.from_dict(edgeDict,orient='index',columns=['TrueEdges','PredEdges'])
    print('Computing PR curve and area')
    prroc = importr('PRROC')
    prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(outDF['PredEdges'].values)), 
              weights_class0 = FloatVector(list(outDF['TrueEdges'].values)))

    print('Computing ROC curve') 
    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    print('Computing PR curve') 
    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)
    
    return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr)


def computeScores_sshv1(trueEdgesDF, predEdgeDF, 
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
            - prec: A list of precision values (for PR plot)
            - recall: A list of precision values (for PR plot)
            - fpr: A list of false positive rates (for ROC plot)
            - tpr: A list of true positive rates (for ROC plot)
            - AUPRC: Area under the precision-recall curve
            - AUROC: Area under the ROC curve
    '''

    if directed:        
        # Initialize dictionaries with all 
        # possible edges
        if selfEdges:
            #print('self started')
            #possibleEdges = list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                             repeat = 2))
            possibleEdges = pd.DataFrame(list(product(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                         repeat = 2)))
            #print('self end')
        else:
            #print('not self started')
            #possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                             r = 2))
            possibleEdges = pd.DataFrame(list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                              r = 2)))
            #print('not self end')

        #print(len(possibleEdges))
        #TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        #PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        TrueEdgeDict = dict.fromkeys(possibleEdges[0]+'|'+possibleEdges[1],0)
        PredEdgeDict = TrueEdgeDict.copy()
        #print(len(TrueEdgeDict))
        
        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth
        #for key in TrueEdgeDict.keys():
        #    if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
        #           (trueEdgesDF['Gene2'] == key.split('|')[1])])>0:
        #            TrueEdgeDict[key] = 1
        trueEdgeKeys = trueEdgesDF['Gene1'].astype(str) + '|' + trueEdgesDF['Gene2']
        for key in trueEdgeKeys:
            if key in TrueEdgeDict:
                TrueEdgeDict[key] = 1
        #print('done TrueEdgeDict')
        
        #print(predEdgeDF.EdgeWeight.values[0])
        #for key in PredEdgeDict.keys():
        #    subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
        #                       (predEdgeDF['Gene2'] == key.split('|')[1])]
        #    if len(subDF)>0:
        #        PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])
        #        print(subDF.EdgeWeight.values)
        predEdgeKeys = predEdgeDF['Gene1'].astype(str) + '|' + predEdgeDF['Gene2']
        predEdgeValues =  np.abs(predEdgeDF.EdgeWeight)
        for k,v in zip(predEdgeKeys,predEdgeValues):
            if k in PredEdgeDict:
                PredEdgeDict[k] = v
        #print('done PredEdgeDict')
    # if undirected
    else:
        
        # Initialize dictionaries with all 
        # possible edges
        if selfEdges:
            #possibleEdges = list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                                                   r = 2))
            possibleEdges = pd.DataFrame(list(combinations_with_replacement(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2)))
        else:
            #possibleEdges = list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
            #                                                   r = 2))
            possibleEdges = pd.DataFrame(list(combinations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
                                                               r = 2)))
        #TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        #PredEdgeDict = {'|'.join(p):0 for p in possibleEdges}
        TrueEdgeDict = dict.fromkeys(possibleEdges[0]+'|'+possibleEdges[1],0)
        PredEdgeDict = TrueEdgeDict.copy()


        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth
        # 0 if edge is not present in the ground-truth

        #for key in TrueEdgeDict.keys():
        #    if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
        #                   (trueEdgesDF['Gene2'] == key.split('|')[1])) |
        #                      ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
        #                   (trueEdgesDF['Gene1'] == key.split('|')[1]))]) > 0:
        #        TrueEdgeDict[key] = 1  
        trueEdgeKeys = pd.concat([trueEdgesDF['Gene1'].astype(str) + '|' + trueEdgesDF['Gene2'],
                                  trueEdgesDF['Gene2'].astype(str) + '|' + trueEdgesDF['Gene1']],
                                 ignore_index=True)
        for key in trueEdgeKeys:
            if key in TrueEdgeDict:
                TrueEdgeDict[key] = 1
        #print('done TrueEdgeDict')

        # Compute PredEdgeDict Dictionary
        # from predEdgeDF

        #for key in PredEdgeDict.keys():
        #    subDF = predEdgeDF.loc[((predEdgeDF['Gene1'] == key.split('|')[0]) &
        #                       (predEdgeDF['Gene2'] == key.split('|')[1])) |
        #                      ((predEdgeDF['Gene2'] == key.split('|')[0]) &
        #                       (predEdgeDF['Gene1'] == key.split('|')[1]))]
        #    if len(subDF)>0:
        #        PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))
        # take max egde weights of two directions
        predEdgeDF_BiDir = pd.concat(predEdgeDF,
                                     predEdgeDF.rename(columns={"Gene1": "Gene2", "Gene2": "Gene1"},
                                                       copy=TRUE),
                                     ignore_index=True)
        predEdgeDF_BiDir['EdgeWeight'] = np.abs(predEdgeDF_BiDir['EdgeWeight'])
        predEdgeDF_max = predEdgeDF_BiDir.groupby(['Gene1','Gene2'])['EdgeWeight'].max().reset_index()
        predEdgeKeys = predEdgeDF_max['Gene1'].astype(str) + '|' + predEdgeDF_max['Gene2']
        predEdgeValues =  predEdgesDF_max.EdgeWeight.values
        for k,v in zip(predEdgeKeys,predEdgeValues):
            if k in PredEdgeDict:
                PredEdgeDict[k] = v
        #print('done PredEdgeDict')

                
    # Combine into one dataframe
    # to pass it to sklearn
    #outDF = pd.DataFrame([TrueEdgeDict,PredEdgeDict]).T
    #outDF.columns = ['TrueEdges','PredEdges']
    outDF = pd.concat([pd.DataFrame.from_dict(TrueEdgeDict,orient='index',columns=['TrueEdges']),
                       pd.DataFrame.from_dict(PredEdgeDict,orient='index',columns=['PredEdges'])],
                      axis=1, join="outer")
    print(outDF.iloc[0:10])
    prroc = importr('PRROC')
    prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(outDF['PredEdges'].values)), 
              weights_class0 = FloatVector(list(outDF['TrueEdges'].values)))

    fpr, tpr, thresholds = roc_curve(y_true=outDF['TrueEdges'],
                                     y_score=outDF['PredEdges'], pos_label=1)

    prec, recall, thresholds = precision_recall_curve(y_true=outDF['TrueEdges'],
                                                      probas_pred=outDF['PredEdges'], pos_label=1)
    
    return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr)


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
    print('Computing PR curve and area in R')
    prroc = importr('PRROC')
    PredEdgeFV = FloatVector(PredEdgeArr)
    TrueEdgeFV = FloatVector(TrueEdgeArr)
    print('   Convert to R obj done; start calculation')
    prCurve = prroc.pr_curve(scores_class0 = PredEdgeFV, 
                             weights_class0 = TrueEdgeFV)

    print('Computing ROC curve') 
    fpr, tpr, thresholds = roc_curve(y_true=TrueEdgeArr,
                                     y_score=PredEdgeArr, pos_label=1)

    print('Computing PR curve')
    prec, recall, thresholds = precision_recall_curve(y_true=TrueEdgeArr,
                                                      probas_pred=PredEdgeArr, pos_label=1)

    #print('Computing average precision score')
    #average_precision = average_precision_score(y_true=TrueEdgeArr,
    #                                            y_score=PredEdgeArr, pos_label=1)

    
    return prec, recall, fpr, tpr, prCurve[2][0], auc(fpr, tpr), auc(recall, prec)


def PRROC3(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
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
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''
    
    # Read file for trueEdges
    print('True edges: '+str(inputSettings.datadir)+'/'+ dataDict['name'] +'/' +dataDict['trueEdges'])
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ dataDict['name'] +
                                '/' +dataDict['trueEdges'],
                                sep = ',', 
                                header = 0, index_col = None)
            
    # Initialize data dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    AveP = {}
    
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
                    precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]], AveP[algo[0]] = computeScores3(trueEdgesDF, predDF, directed = True, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')                    
            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/' + runDir + 'PRplot'
            ROCName = '/' + runDir + 'ROCplot'
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
                    precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]], AveP[algo[0]] = computeScores3(trueEdgesDF, predDF, directed = False, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')
            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
            PRName = '/' + runDir + 'uPRplot'
            ROCName = '/' + runDir + 'uROCplot'
    if (plotFlag):
         ## Make PR curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(recallDict[key],precisionDict[key], ci=None)
            legendList.append(key + ' (AUPRC = ' + str("%.2f" % (AUPRC[key]))+')')
        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(legendList) 
        plt.savefig(outDir+PRName+'.pdf')
        plt.savefig(outDir+PRName+'.png')
        plt.clf()

        ## Make ROC curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(FPRDict[key],TPRDict[key], ci=None)
            legendList.append(key + ' (AUROC = ' + str("%.2f" % (AUROC[key]))+')')

        plt.plot([0, 1], [0, 1], linewidth = 1.5, color = 'k', linestyle = '--')

        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.legend(legendList) 
        plt.savefig(outDir+ROCName+'.pdf')
        plt.savefig(outDir+ROCName+'.png')
        plt.clf()
    return AUPRC, AUROC, AveP


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
            - prec: A list of precision values (for PR plot)
            - recall: A list of precision values (for PR plot)
            - fpr: A list of false positive rates (for ROC plot)
            - tpr: A list of true positive rates (for ROC plot)
            - AUPRC: Area under the precision-recall curve
            - AUROC: Area under the ROC curve
    '''
    uniqNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
    return computeScoresByAdjMat(uniqNodes,trueEdgesDF, predEdgeDF,
                                 directed, selfEdges)



def PRROC4(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False):
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
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
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
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    AveP = {}
    
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
                    precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]], AveP[algo[0]] = computeScores4(expressionData,trueEdgesDF, predDF, directed = True, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')                    
            else:
                print(outDir + '/' +algo[0] + '/' + runDir +'rankedEdges.csv', \
                      ' does not exist. Skipping...')
            PRName = '/' + runDir + 'PRplot'
            ROCName = '/' + runDir + 'ROCplot'
    else:
        for algo in tqdm(inputSettings.algorithms, 
                         total = len(inputSettings.algorithms), unit = " Algorithms"):

            print(algo)
            # check if the output rankedEdges file exists
            if Path(outDir + '/' +algo[0]+ '/' + runDir + 'rankedEdges.csv').exists():
                 # Initialize Precsion
                predDF = pd.read_csv(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                                            sep = '\t', header =  0, index_col = None)
                print(trueEdgesDF.head())
                print(predDF.head())
                if predDF.shape[0] > 0:
                    precisionDict[algo[0]], recallDict[algo[0]], FPRDict[algo[0]], TPRDict[algo[0]], AUPRC[algo[0]], AUROC[algo[0]], AveP[algo[0]] = computeScores4(expressionData, trueEdgesDF, predDF, directed = False, selfEdges = selfEdges)
                else:
                    print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                          ' is empty. Skipping...')
            else:
                print(outDir + '/' +algo[0] + '/' + runDir + 'rankedEdges.csv', \
                  ' does not exist. Skipping...')
            
            PRName = '/' + runDir + 'uPRplot'
            ROCName = '/' + runDir + 'uROCplot'
    if (plotFlag):
         ## Make PR curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(recallDict[key],precisionDict[key], ci=None)
            legendList.append(key + ' (AUPRC = ' + str("%.2f" % (AUPRC[key]))+')')
        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(legendList) 
        plt.savefig(outDir+PRName+'.pdf')
        plt.savefig(outDir+PRName+'.png')
        plt.clf()

        ## Make ROC curves
        legendList = []
        for key in recallDict.keys():
            sns.lineplot(FPRDict[key],TPRDict[key], ci=None)
            legendList.append(key + ' (AUROC = ' + str("%.2f" % (AUROC[key]))+')')

        plt.plot([0, 1], [0, 1], linewidth = 1.5, color = 'k', linestyle = '--')

        plt.xlim(0,1)    
        plt.ylim(0,1)
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.legend(legendList) 
        plt.savefig(outDir+ROCName+'.pdf')
        plt.savefig(outDir+ROCName+'.png')
        plt.clf()
    return AUPRC, AUROC, AveP


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
            - prec: A list of precision values (for PR plot)
            - recall: A list of precision values (for PR plot)
            - fpr: A list of false positive rates (for ROC plot)
            - tpr: A list of true positive rates (for ROC plot)
            - AUPRC: Area under the precision-recall curve
            - AUROC: Area under the ROC curve
    '''
    uniqNodes = np.unique(expressionData.index)
    return computeScoresByAdjMat(uniqNodes,trueEdgesDF, predEdgeDF,
                                 directed, selfEdges)

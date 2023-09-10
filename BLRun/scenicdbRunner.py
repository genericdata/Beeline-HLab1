import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil

add_bind = ''
RANKING_DBS_FNAMES = {'feather_v1':
                      {'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr':'resources/cistarget/databases/feather_v1/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather',
                       'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr':'resources/cistarget/databases/feather_v1/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather',
                      'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr':'resources/cistarget/databases/feather_v1/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather',
                       'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr':'resources/cistarget/databases/feather_v1/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'
                      }
                      }
MOTIF_ANNOTATIONS_FNAMES = {'motifs-v9-nr.hgnc-m0.001-o0.0.tbl':'resources/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
                            'motifs-v9-nr.mgi-m0.001-o0.0.tbl':'resources/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl'}
TFS_FNAMES = {'allTFs_hg38':'resources/cistarget/tf_lists/allTFs_hg38.txt',
             'allTFs_mm':'resources/cistarget/tf_lists/allTFs_mm.txt'}

MATCH_PRIORITIES_EXCLUDE = 99 # exclude the match types
HGNC_MATCH_TYPE_PRIORITIES = {'Approved symbol': 1, 'Previous symbol': 2, 'Alias symbol': 3,
                              'Entry withdrawn': 4, 'Unmatched':MATCH_PRIORITIES_EXCLUDE }
MGI_INPUT_TYPE_PRIORITIES = {'current symbol':1, 'old symbol': 2, 'MGI': 3, 'synonym': 4, 'related synonym': 5,
                             'Genbank':6, 'human symbol':MATCH_PRIORITIES_EXCLUDE,'human synonym':MATCH_PRIORITIES_EXCLUDE,
                             'rat symbol':MATCH_PRIORITIES_EXCLUDE,'rat synonym':MATCH_PRIORITIES_EXCLUDE,
                             'zebrafish symbol':MATCH_PRIORITIES_EXCLUDE}

def parseIdMapFile(idMapFile,idMapType):
    if idMapType=='HGNCSymbolCheck230830.csv':
        idMapDf = pd.read_csv(idMapFile,sep=',',skiprows=1)
        idMapDf['InputID'] = idMapDf['Input']
        idMapDf['UnifiedID'] = idMapDf['Approved symbol']
        idMapDf['InputIDType'] = idMapDf['Match type']
        idMapDf = idMapDf[idMapDf['InputIDType'].notnull()]
        idMapDf['MatchPriority'] = idMapDf['InputIDType'].map(HGNC_MATCH_TYPE_PRIORITIES)
    elif idMapType=='MGIBatchReport230830.tsv':
        idMapDf = pd.read_csv(idMapFile,sep='\t')
        idMapDf['InputID'] = idMapDf['Input']
        idMapDf['UnifiedID'] = idMapDf['Symbol']
        idMapDf['InputIDType'] = idMapDf['Input Type']
        idMapDf = idMapDf[idMapDf['InputIDType'].notnull()]
        idMapDf['MatchPriority'] = idMapDf['InputIDType'].map(MGI_INPUT_TYPE_PRIORITIES)
    else:
        print('Do not know how to parse ID map type for expression data:',RunnerObj.idMapAllExprDataType,'. Exiting')
        sys.exit(1)

    idMapDfOut = idMapDf[['InputID','InputIDType','UnifiedID','MatchPriority']].drop_duplicates()
    idMapDfOut = idMapDfOut.loc[idMapDfOut.groupby(['InputID'])['MatchPriority'].idxmin()]# take the lowest priority matches for the same InputID
    idMapDfOut = idMapDfOut.loc[idMapDfOut.groupby(['UnifiedID'])['MatchPriority'].idxmin()]# take the lowest priority matches for each UnifiedID
    return idMapDfOut

def processExpressionData(exprDataFile,exprDataIdMapFile,exprDataIdMapType,
                          dbIdMapFiles,dbIdMapType,tfFile,mapTFIds):
                          
    print('Expression data ID map: ',str(exprDataIdMapFile))
    print('Expression data ID map type: ',exprDataIdMapType)
    idMapExprData = parseIdMapFile(exprDataIdMapFile,exprDataIdMapType)

    print('DB ID map type: ',dbIdMapType)
    idMapDbList = []
    for dbFile in dbIdMapFiles:
        idMapDb0 = parseIdMapFile(dbFile,dbIdMapType)
        idMapDbList.append(idMapDb0)
    idMapDb = pd.concat(idMapDbList).drop_duplicates()

    print('Expression data: ',str(exprDataFile))
    ExpressionData = pd.read_csv(exprDataFile,header = 0, index_col = 0)
    ExpressionData['ExprDataID'] = ExpressionData.index

    # convert to Unified ID then to the ID in the DB
    ExpressionData_Map1 = pd.merge(ExpressionData,idMapExprData[idMapExprData['UnifiedID'].notnull()],
                                   left_on='ExprDataID',right_on='InputID')
    ExpressionData_Map2 = pd.merge(ExpressionData_Map1,idMapDb[idMapDb['UnifiedID'].notnull()],
                                   left_on='UnifiedID',right_on='UnifiedID',
                                   suffixes=['_ExprData','_Db'])
    ExpresssionData_Map_Out = ExpressionData_Map2[['UnifiedID','ExprDataID','InputID_ExprData','InputIDType_ExprData','MatchPriority_ExprData',
                                                   'InputID_Db','InputIDType_Db','MatchPriority_Db']]
    print('Exclude these IDs from ExpressionData that have no acceptable match types. See ID map file for details.')
    print(ExpressionData_Map2[(ExpressionData_Map2.MatchPriority_ExprData==MATCH_PRIORITIES_EXCLUDE) |
                              (ExpressionData_Map2.MatchPriority_Db==MATCH_PRIORITIES_EXCLUDE)]['ExprDataID'])
    ExpressionData_Map3 = ExpressionData_Map2[(ExpressionData_Map2.MatchPriority_ExprData!=MATCH_PRIORITIES_EXCLUDE) &
                                              (ExpressionData_Map2.MatchPriority_Db!=MATCH_PRIORITIES_EXCLUDE)]
    ExpressionData_Out = ExpressionData_Map3.set_index('InputID_Db').drop(
        columns=['UnifiedID','ExprDataID','InputID_ExprData',
                 'InputIDType_ExprData','InputIDType_Db',
                 'MatchPriority_ExprData','MatchPriority_Db'])

    tfNames = pd.read_csv(tfFile,names=['TFID'])
    if (mapTFIds):
        tfNames_Map1 = pd.merge(tfNames,idMapExprData[idMapExprData['UnifiedID'].notnull()],
                                left_on='TFID',right_on='InputID')
        tfNames_Map2 = pd.merge(tfNames_Map1,idMapDb[idMapDb['UnifiedID'].notnull()],
                                left_on='UnifiedID',right_on='UnifiedID',
                                suffixes=['_TFName','_Db'])
        tfNames_Map_Out = tfNames_Map2[['TFID','InputID_TFName','InputID_Db','UnifiedID']]
        tfNames_Map2 = tfNames_Map2.set_index('InputID_Db')
        tfNames_Out = tfNames_Map2.index.to_series()
    else:
        tfNames_Map_Out = None
        tfNames_Out = None
    
    return ExpresssionData_Map_Out,ExpressionData_Out,tfNames_Map_Out,tfNames_Out


def generateInputs(RunnerObj):

    algName = RunnerObj.name
    datasetName = RunnerObj.datasetName
    runDir = RunnerObj.params.get('run_dir','')
    dbVersion = RunnerObj.dataset_params['rankingDbVersion']
    dbList = [RANKING_DBS_FNAMES[dbVersion][dbName] for dbName in RunnerObj.dataset_params['rankingDbNames'][datasetName].split(',')]

    if not RunnerObj.inputDir.joinpath(algName,runDir).exists():
        print("Input folder for "+algName+" "+runDir+" does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath(algName,runDir).mkdir(exist_ok = False)

    print('TF names data: ',TFS_FNAMES[RunnerObj.dataset_params['tfFName'][datasetName]])
    print('mapTFIds:',RunnerObj.dataset_params['mapTFIds'][datasetName])
    exprDataIdMapOut,exprDataOut,tfNamesMapOut,tfNamesOut =\
        processExpressionData(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                              RunnerObj.inputDir.joinpath(RunnerObj.exprDataIdMap),
                              RunnerObj.exprDataIdMapType,
                              [dbFile+'.genes.'+RunnerObj.dataset_params['rankingDbIdMapType'][datasetName] for dbFile in dbList],
                              RunnerObj.dataset_params['rankingDbIdMapType'][datasetName],
                              TFS_FNAMES[RunnerObj.dataset_params['tfFName'][datasetName]],      
                              RunnerObj.dataset_params['mapTFIds'][datasetName])
    exprDataIdMapOut.to_csv(RunnerObj.inputDir.joinpath(algName,"ExpressionData-IDMap.csv"),
                            sep=',',header=True,index=False)
    exprDataOut.to_csv(RunnerObj.inputDir.joinpath(algName,"ExpressionData.csv"),
                       sep = ',', header  = True, index = True,index_label=False)
    if (RunnerObj.dataset_params['mapTFIds'][datasetName]):
        tfNamesMapOut.to_csv(RunnerObj.inputDir.joinpath(algName,"TFs-IDMap.csv"),
                             sep=',',header=True,index=False)
        tfNamesOut.to_csv(RunnerObj.inputDir.joinpath(algName,"TFs.txt"),
                   sep = ',', header=False, index =False)
    else:
        shutil.copyfile(TFS_FNAMES[RunnerObj.dataset_params['tfFName'][datasetName]],
                        RunnerObj.inputDir.joinpath(algName,'TFs.txt'))

def run(RunnerObj):
    '''
    Function to run SCENIC algorithm using the cisTarget databases

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # current path will be mounted as /ext3/data
    # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
    algName = RunnerObj.name
    datasetName = RunnerObj.datasetName
    dbVersion = RunnerObj.dataset_params['rankingDbVersion']
    dbList = [RANKING_DBS_FNAMES[dbVersion][dbName] for dbName in RunnerObj.dataset_params['rankingDbNames'][datasetName].split(',')]
    dbs_param = ' '.join(dbList)
    motif_fname = MOTIF_ANNOTATIONS_FNAMES[RunnerObj.dataset_params['motifFName'][datasetName]]

    random_seed_list = [int(seed) for seed in RunnerObj.params['randomSeedList'].split(',')]
    for seed in random_seed_list:
            
        runDir = 'random_seed_'+str(seed)
        algName = RunnerObj.name
        inputDirB = Path('/ext3/data').joinpath(RunnerObj.inputDir.relative_to(Path.cwd()),algName)

        # outputs/L0/hESC/DEEPDRIM
        # make output dirs if they do not exist:
        outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"+runDir+"/"
        outDirB = Path('/ext3/data').joinpath(outDir)
        os.makedirs(outDir, exist_ok = True)

        cmdToRun = ' '.join([
            'singularity exec --no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate SCENIC; ',
            'cd',str(Path('/ext3/data')),';'
            'command time -v -o '+str(outDirB.joinpath('time_grn.txt')),
            'pyscenic grn',
            str(inputDirB.joinpath('ExpressionData.csv')),
            str(inputDirB.joinpath('TFs.txt')),
            '-m grnboost2',
            '--seed '+str(seed),
            '-t',
            '-o ',str(outDirB.joinpath('adjacencies.tsv')),
            '--num_workers '+RunnerObj.params['numWorkers'],';',
            'command time -v -o '+str(outDirB.joinpath('time_ctx.txt')),
            'pyscenic ctx',
            str(outDirB.joinpath('adjacencies.tsv')),
            dbs_param,
            '--annotations_fname '+motif_fname,
            '--expression_mtx_fname '+str(inputDirB.joinpath('ExpressionData.csv'))+' -t',
            '--rank_threshold '+RunnerObj.params['rankThreshold'] if RunnerObj.params['rankThreshold']!='None' else '',
            '--nes_threshold '+RunnerObj.params['NESThreshold'] if RunnerObj.params['NESThreshold']!='None' else '',
            '--output '+str(outDirB.joinpath('motifs.csv')),
            '--num_workers '+RunnerObj.params['numWorkers'],
            '\"'])
        print(cmdToRun)
        os.system(cmdToRun)

def processRankedEdgesDb(rankedEdgesDbFile,idMapExprData,rankedEdgesFile):

    rankedEdgesDb = pd.read_csv(rankedEdgesDbFile,sep='\t',header=0)
    rankedEdgesOut1 = pd.merge(rankedEdgesDb,idMapExprData[['InputID_Db','ExprDataID']],
                                  left_on='Gene1',right_on='InputID_Db')
    rankedEdgesOut2 = pd.merge(rankedEdgesOut1,idMapExprData[['InputID_Db','ExprDataID']],
                               left_on='Gene2',right_on='InputID_Db',
                               suffixes=['1','2'])
    rankedEdgesOut = rankedEdgesOut2[['ExprDataID1','ExprDataID2','EdgeWeight']].copy()
    rankedEdgesOut = rankedEdgesOut.sort_values(by=['EdgeWeight','ExprDataID2','ExprDataID1'], ascending=[False,True,True])

    rankedEdgesOut.to_csv(rankedEdgesFile,
                          sep='\t',index=False,header=['Gene1','Gene2','EdgeWeight'])

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCENIC. 

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    random_seed_list = [int(seed) for seed in RunnerObj.params['randomSeedList'].split(',')]

    for seed in random_seed_list:

        runDir = 'random_seed_'+str(seed)
        # current path will be mounted as /ext3/data
        # '/ext3/data/inputs/L0/hESC/DEEPDRIM'
        algName = RunnerObj.name

        # outputs/L0/hESC/DEEPDRIM
        # make output dirs if they do not exist:
        outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/"+algName+"/"+runDir+"/"
        outDirB = Path('/ext3/data').joinpath(outDir)

        cmdToRun = ' '.join([
            'singularity exec --no-home',
            '-B ' + str(Path.cwd())+':/ext3/data/'+add_bind,
            '--overlay ' + str(RunnerObj.singularityOverlay) + ':ro',
            str(RunnerObj.singularityImage),
            '/bin/bash -c \"source /ext3/env.sh; conda activate SCENIC; ',
            'cd',str(Path('/ext3/data')),';'
            'command time -v -o '+str(outDirB.joinpath('time_rankededges.txt')),
            'python Algorithms/SCENICDB/scenicdf_to_rankededges.py',
            '--scenic_motifs_fname='+str(outDirB.joinpath('motifs.csv')),
            '--ranked_edges_fname='+str(outDirB.joinpath('rankedEdges_Db.csv')),
            '\"'])
        print(cmdToRun)
        os.system(cmdToRun)

        processRankedEdgesDb(Path(outDir).joinpath('rankedEdges_Db.csv'),
                             pd.read_csv(RunnerObj.inputDir.joinpath(algName,'ExpressionData-IDMap.csv'),sep=',',header=0),
                             Path(outDir).joinpath('rankedEdges.csv'))

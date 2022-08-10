import glob
import os
import sys
import datetime
import pandas as pd
import numpy as np
import shutil
import yaml

BEELINE_NETWORKS_DIR='/scratch/ch153/packages/BEELINE/BEELINE-Networks'
BEELINE_DATA_DIR='/scratch/ch153/packages/BEELINE/BEELINE-data'
SERGIO_DATASETS_DIR='/scratch/ch153/packages/SERGIO/PayamDiba/SERGIO/data_sets'

##########################
### Datasets
##########################
### BEELINE real datasets
#DATASETS = ['hESC','hHEP','mDC','mESC','mHSC-E']
#DATASETS_NS = ['hESC','hHEP','mDC','mESC','mHSC-E']
#DATASET_DIRS = ['L0','L1','L2']
#DATASET_DIRS_NS = ['L0_ns','L1_ns','L2_ns']
#DATASET_DIRS_LOFGOF = ['L0_lofgof','L1_lofgof','L2_lofgof']

ALGORITHMS = ['CICT','DEEPDRIM','GENIE3','GRISLI','GRNBOOST2','GRNVBEM',
              'LEAP','PIDC','PPCOR','RANDOM','SCNS','SCODE','SCRIBE',
              'SINCERITIES','SINGE']

# for each ground-truth network
DATASET_PARAMS = {'cs':{'dataset':['hESC','hHep','mDC','mESC','mHSC-E','mHSC-GM','mHSC-L'],
                       'netfile':['human/hESC-ChIP-seq-network.csv',
                                  'human/HepG2-ChIP-seq-network.csv',
                                  'mouse/mDC-ChIP-seq-network.csv',
                                  'mouse/mESC-ChIP-seq-network-fixed.csv',
                                  'mouse/mHSC-ChIP-seq-network.csv',
                                  'mouse/mHSC-ChIP-seq-network.csv',
                                  'mouse/mHSC-ChIP-seq-network.csv'],
                       'tffile':['human-tfs-upper.csv','human-tfs-upper.csv',
                                 'mouse-tfs.csv','mouse-tfs.csv','mouse-tfs.csv',
                                 'mouse-tfs.csv','mouse-tfs.csv']},
                  'ns': {'dataset':['hESC','hHep','mDC','mESC','mHSC-E','mHSC-GM','mHSC-L'],
                         'netfile':['human/Non-specific-ChIP-seq-network.csv',
                                    'human/Non-specific-ChIP-seq-network.csv',
                                    'mouse/Non-Specific-ChIP-seq-network.csv',
                                    'mouse/Non-Specific-ChIP-seq-network.csv',
                                    'mouse/Non-Specific-ChIP-seq-network.csv',
                                    'mouse/Non-Specific-ChIP-seq-network.csv',
                                    'mouse/Non-Specific-ChIP-seq-network.csv'],
                         'tffile':['human-tfs-upper.csv','human-tfs-upper.csv',
                                   'mouse-tfs.csv','mouse-tfs.csv','mouse-tfs.csv',
                                   'mouse-tfs.csv','mouse-tfs.csv']},
                  'lofgof':{'dataset':['mESC'],
                            'netfile':['mouse/mESC-lofgof-network.csv'],
                            'tffile':['mouse-tfs.csv']}
                  }
                  
exp_param_list = []
EXP_LEVELS = {'L0':{'num_genes':500,'num_rand_genes':0},
              'L1':{'num_genes':1000,'num_rand_genes':0},
              'L2':{'num_genes':500,'num_rand_genes':500}
              }
for level_name,level_settings in EXP_LEVELS.items():
    for gt_name,gt_settings in DATASET_PARAMS.items():
        dataset_df = pd.DataFrame(data=gt_settings)
        exp_dir = level_name if gt_name=='cs' else level_name+'_'+gt_name
        dataset_df['exp_dir'] = exp_dir
        dataset_df['num_genes'] = level_settings['num_genes']
        dataset_df['num_rand_genes'] = level_settings['num_rand_genes']
        exp_param_list.append(dataset_df)
EXP_PARAM_DF = pd.concat(exp_param_list)
#print(EXP_PARAM_DF)

def get_exp_file(exp_param_df,exp_dir,dataset,upath):
    ds_df = exp_param_df[(exp_param_df.exp_dir==exp_dir) & (exp_param_df.dataset==dataset)].iloc[[0]]
    print(ds_df)
    filepath = expand(upath,u=ds_df.itertuples())
    return filepath

rule beeline_exp_inputs:
    input: expfile=os.path.join(BEELINE_DATA_DIR,'inputs/scRNA-Seq/{dataset}/ExpressionData-upper.csv'),\
           geneorderingfile=os.path.join(BEELINE_DATA_DIR,'inputs/scRNA-Seq/{dataset}/GeneOrdering-upper.csv'),\
           netfile=lambda wc: get_exp_file(EXP_PARAM_DF,wc.exp_dir,wc.dataset,\
                                           os.path.join(BEELINE_NETWORKS_DIR,'Networks/{u.netfile}')),\
           tffile=lambda wc: get_exp_file(EXP_PARAM_DF,wc.exp_dir,wc.dataset,\
                                          os.path.join(BEELINE_NETWORKS_DIR,'Networks/{u.tffile}')),\
           ptfile=os.path.join(BEELINE_DATA_DIR,'inputs/scRNA-Seq/{dataset}/PseudoTime.csv')
    output: expout='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
            netout='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-network.csv',\
            ptout='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-PseudoTime.csv'
    log: '{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}.log'
    params: overlay="conda_greene/overlay-5GB-200K-beeline20211104.ext3",\
            sif="conda_greene/centos-8.2.2004.sif",\
            D="{exp_input_dir}/{exp_dir}/{dataset}",\
            p='-p 0.01',c='-c',t='-t',\
            n=lambda wc: get_exp_file(EXP_PARAM_DF,wc.exp_dir,wc.dataset,'{u.num_genes}'),\
            r=lambda wc: get_exp_file(EXP_PARAM_DF,wc.exp_dir,wc.dataset,'{u.num_rand_genes}'),\
            outprefix="{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}",\
            jobname="blg_{exp_dir}-{dataset}",\
            clog_prefix="{exp_dir}"
    threads: 1
    shell: """
    mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
    singularity exec \
    --overlay {params.overlay}:ro \
    {params.sif} \
    /bin/bash -c \"source /ext3/env.sh; conda activate BEELINE; \
    python generateExpInputs.py -e {input.expfile} -g {input.geneorderingfile} -f {input.netfile} -i {input.tffile} {params.p} {params.c} {params.t} -n {params.n} -r {params.r} -o {params.outprefix} > {log} 2>&1 \"
    cp {input.ptfile} {output.ptout}
    rm -f {params.D}/ExpressionData.csv; ln -rs {output.expout} {params.D}/ExpressionData.csv
    rm -f {params.D}/refNetwork.csv; ln -rs {output.netout} {params.D}/refNetwork.csv
    rm -f {params.D}/PseudoTime.csv; ln -rs {output.ptout} {params.D}/PseudoTime.csv
    """
rule beeline_exp_inputs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/{ds.exp_dir}-ExpressionData.csv',\
                 ds=EXP_PARAM_DF.itertuples())

### SERGIO data from github
### DS1, DS2, DS3 are steady state
### DS4, DS5, DS6, DS7 are differentiation datasets
###      formatted by notebooks/ipynb/format_sergio_2022-07-05.ipynb
SERGIO_DATASET_PARAMS = {
    'exp_name':['DS1','DS2','DS3'],
    'exp_dir':['SERGIO_DS1','SERGIO_DS2','SERGIO_DS3'],
    'source_dir':['De-noised_100G_9T_300cPerT_4_DS1',
                  'De-noised_400G_9T_300cPerT_5_DS2',
                  'De-noised_1200G_9T_300cPerT_6_DS3'],
    'expfile':['simulated_noNoise_{network_i}.csv',
               'simulated_noNoise_{network_i}.csv',
               'simulated_noNoise_{network_i}.csv'],
    'ncells_per_cl':[300,300,300],
    'nclusters':[9,9,9]}

SERGIO_PARAMS_DF = pd.DataFrame(data=SERGIO_DATASET_PARAMS)
def get_sergio_file(sergio_param_df,exp_name,upath):
    ds_df = sergio_param_df[(sergio_param_df.exp_name==exp_name)].iloc[[0]]
    print(ds_df)
    filepath = expand(upath,u=ds_df.itertuples())
    return filepath

rule beeline_sergio_dsdata:# data common for all simulated files in each dataset
    input: netfile=lambda wc: get_sergio_file(SERGIO_PARAMS_DF,wc.exp_name,\
                                           os.path.join(SERGIO_DATASETS_DIR,'{u.source_dir}/gt_GRN.csv'))[0]
    params: ncells_per_cl=lambda wc: int(get_sergio_file(SERGIO_PARAMS_DF,wc.exp_name,'{u.ncells_per_cl}')[0]),\
            nclusters=lambda wc: int(get_sergio_file(SERGIO_PARAMS_DF,wc.exp_name,'{u.nclusters}')[0])
    output: netout='{exp_input_dir}/SERGIO_{exp_name}/refNetwork.csv',\
            clout='{exp_input_dir}/SERGIO_{exp_name}/ClusterIds.csv',\
            ptout='{exp_input_dir}/SERGIO_{exp_name}/PseudoTime.csv'
    threads: 1
    run:
        grn = pd.read_csv(input.netfile, names=['Gene1','Gene2'],index_col=None)
        grn['Gene1'] = 'G' + grn['Gene1'].astype(str)
        grn['Gene2'] = 'G' + grn['Gene2'].astype(str)
        grn['Type'] = 1
        grn.to_csv(output.netout,sep=',',index=None)
        
        cell_type = np.repeat(np.array([i for i in range(params.nclusters)]),params.ncells_per_cl)
        cell_ids = ['E'+str(i) for i in range(params.nclusters*params.ncells_per_cl)]
        cluster_ids = pd.DataFrame(data={'cl':cell_type},index=cell_ids)
        cluster_ids.to_csv(output.clout,sep=',')

        cell_id_per_type = np.tile(np.array([float(i) for i in range(params.ncells_per_cl)]),params.nclusters)
        pt = pd.DataFrame(data={'Time':cell_id_per_type,'cl':cell_type},index=cell_ids)
        pt_grp = pt.groupby('cl').Time
        pt['PseudoTime'] = (pt.Time - pt_grp.transform('min'))/pt_grp.transform(np.ptp)
        pt.to_csv(output.ptout,sep=',',columns=['PseudoTime','Time'],header=True,index_label='Cell ID')
        
rule beeline_sergio_dsdata_out:
   input: expand('inputs_beeline2/SERGIO_{ds.exp_name}/refNetwork.csv',\
                 ds=SERGIO_PARAMS_DF.itertuples())

rule beeline_sergio_bulk_dsdata:# data common for alal simulated files in each dataset
    input: netfile='{exp_input_dir}/SERGIO_{exp_name}/refNetwork.csv'
    output: netout='{exp_input_dir}/SERGIO_{exp_name}_bulk/refNetwork.csv'
    threads: 1
    shell: """
        cp {input.netfile} {output.netout}
    """
rule beeline_sergio_bulk_dsdata_out:
   input: expand('inputs_beeline2/SERGIO_{ds.exp_name}_bulk/refNetwork.csv',\
                 ds=SERGIO_PARAMS_DF.itertuples())

rule beeline_sergio_inputs:
    input: expfile=lambda wc: get_sergio_file(SERGIO_PARAMS_DF,wc.exp_name,\
                                              os.path.join(SERGIO_DATASETS_DIR,'{u.source_dir}',\
                                                           SERGIO_PARAMS_DF[SERGIO_PARAMS_DF.exp_name==wc.exp_name].expfile.iloc[0].format(network_i=wc.network_i)))[0],\
           clfile='{exp_input_dir}/SERGIO_{exp_name}/ClusterIds.csv',\
           netfile='{exp_input_dir}/SERGIO_{exp_name}/refNetwork.csv',\
           ptfile='{exp_input_dir}/SERGIO_{exp_name}/PseudoTime.csv'
    params: ncells_per_cl=lambda wc: int(get_sergio_file(SERGIO_PARAMS_DF,wc.exp_name,'{u.ncells_per_cl}')[0]),\
            nclusters=lambda wc: int(get_sergio_file(SERGIO_PARAMS_DF,wc.exp_name,'{u.nclusters}')[0])
    output: expout_csv='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/ExpressionData.csv',\
            expoutT_csv='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/ExpressionDataT.csv',\
            netout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/refNetwork.csv',\
            ptout='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/PseudoTime.csv'
    threads: 1
    run:
        exp_from = pd.read_csv(input.expfile, index_col=0)
        colnamemap = {colname:'E'+str(colind) for colind,colname in enumerate(exp_from.columns)}
        exp = exp_from.rename(index=lambda s: "G"+str(s),columns=colnamemap)
        exp.to_csv(output.expout_csv,sep=',')
        expT = exp.T
        expT.to_csv(output.expoutT_csv,index=None)

        shutil.copy(input.netfile,output.netout)

        clusterdf = pd.read_csv(input.clfile,index_col=0)
        numclusters = clusterdf.cl.unique()
        pt_in = pd.read_csv(input.ptfile,index_col=0)
        pt = pt_in.join(clusterdf)
        pt_grp = pt.groupby('cl').PseudoTime
        pt_norm = (pt.PseudoTime - pt_grp.transform('min'))/pt_grp.transform(np.ptp)

        pt_out = pd.DataFrame(index=exp.columns,
                              columns = ['PseudoTime' + str(1+i) for i in range(len(numclusters))])
        for (cid,row), pt in zip(pt_out.iterrows(), pt_norm):
            pt_out.loc[cid]['PseudoTime'+str(clusterdf.loc[cid]['cl']+1)] = pt
        pt_out.to_csv(output.ptout,na_rep='NA')
            
rule beeline_sergio_inputs_out:
   input: expand('inputs_beeline2/SERGIO_{ds.exp_name}/net{network_i}/ExpressionData.csv',\
                 ds=SERGIO_PARAMS_DF.itertuples(),network_i=range(1))

rule beeline_sergio_bulk_inputs:
    input: expfileT='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/ExpressionDataT.csv'
    output: expout_tsv='{exp_input_dir}/SERGIO_{dataset_name,[^_]+}_bulk/net{network_i,[0-9]+}/ExpressionData.csv'
    threads: 1
    run:
        grn = pd.read_csv(input.expfile,sep=',')
        cell_type = np.repeat(np.array([i for i in range(9)]),300)
        grn['cell_type'] = cell_type
        print(grn.head())
        grn = grn.groupby('cell_type').mean().reset_index().drop(columns=['cell_type'])
        grn.to_csv(output.expout_tsv,sep='\t', index = None)
        
rule beeline_sergio_bulk_inputs_out:
   input: expand('inputs_beeline2/SERGIO_{ds.exp_name}_bulk/net{network_i}/simulated_noNoise_bulk_{network_i}.tsv',\
                 ds=SERGIO_PARAMS_DF.itertuples(),network_i=range(15))

rule beeline_sergio_netdata_inputs:
    input: netfile='{exp_input_dir}/SERGIO_{dataset_name}/refNetwork.csv'
    output: netout='{exp_input_dir}/SERGIO_{dataset_name}/net{network_i,[0-9]+}/refNetwork.csv'
    threads: 1
    shell: """
        cp {input.netfile} {output.netout}
    """
rule beeline_sergio_netdata_inputs_out:
   input: expand('inputs_beeline2/SERGIO_{ds.exp_name}/net{network_i}/refNetwork.tsv',\
                 ds=SERGIO_PARAMS_DF.itertuples(),network_i=range(15)) +\
          expand('inputs_beeline2/SERGIO_{ds.exp_name}_bulk/net{network_i}/refNetwork.tsv',\
                 ds=SERGIO_PARAMS_DF.itertuples(),network_i=range(15))
   
##########################
### BEELINE run 
##########################
# GRISLI 240000
def get_run_mem_mb(wildcards):
    switcher = {
        "L0": 48000,
        "L1": 48000,
        "L2": 48000,
        "L0_ns": 48000,
        "L1_ns": 48000,
        "L2_ns": 48000,
        "L0_lofgof": 32000,
        "L1_lofgof": 32000,
        "L2_lofgof": 32000,
        "SERGIO_DS4":32000,        
        "SERGIO_DS5":32000,
        "SERGIO_DS6":32000,
        "SERGIO_DS7":32000
        }
    return switcher.get(wildcards.exp_dir, 8000)

def get_run_time(wildcards):
    switcher = {
        "L0": "3:00:00",
        "L1": "3:00:00",
        "L2": "3:00:00",
        "L0_ns": "3:00:00",
        "L1_ns": "3:00:00",
        "L2_ns": "3:00:00",
        "L0_lofgof": "3:00:00",
        "L1_lofgof": "3:00:00",
        "L2_lofgof": "3:00:00",
        "SERGIO_DS4": "3:00:00",
        "SERGIO_DS5": "3:00:00",
        "SERGIO_DS6": "3:00:00",
        "SERGIO_DS7": "3:00:00"
        }
    return switcher.get(wildcards.exp_dir, "03:00:00")

def get_run_gpu(wildcards):
    switcher = {
        'DEEPDRIM': '--gres=gpu:1'
        }
    return switcher.get(wildcards.algorithm, '')

# split config
#parallel --dry-run "python generateSplitConfigs.py --config config-files/config_{}.yaml --config_split_dir config-files-split/config_{}_split --replace_existing" ::: L0 L0_ns L0_lofgof L1 L1_ns L1_lofgof L2 L2_ns L2_lofgof SERGIO_DS4 SERGIO_DS5 SERGIO_DS6 SERGIO_DS7
rule beeline_run:
    input: 'config-files-split/config_{exp_dir}_split/{dataset}/{algorithm}/config.yaml'
    output: 'outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.log'
    params: D="outputs/{exp_dir}/{dataset}/{algorithm}",\
            jobname="blr_{exp_dir}-{dataset}-{algorithm}",\
            clog_prefix="rankedEdges"
    threads: 1
    resources: mem_mb=get_run_mem_mb, time=get_run_time, gpu=get_run_gpu
    envmodules: "anaconda3/2020.07"
    shell: """
        mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
        python BLRunner.py --config {input} > {output} 2>&1 || true
    """
rule beeline_run_out:
   input: expand('outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.log',\
                 exp_dir=['L0','L1','L2'],\
                 dataset=DATASET_PARAMS['cs']['dataset'],\
                 algorithm=['DEEPDRIM'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_out
rule beeline_run_ns_out:
   input: expand('outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.log',\
                 exp_dir=['L0_ns','L1_ns','L2_ns'],\
                 dataset=DATASET_PARAMS['ns']['dataset'],\
                 algorithm=['DEEPDRIM'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_ns_out
rule beeline_run_lofgof_out:
   input: expand('outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.log',\
                 exp_dir=['L0_lofgof','L1_lofgof','L2_lofgof'],\
                 dataset=DATASET_PARAMS['lofgof']['dataset'],\
                 algorithm=['DEEPDRIM'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_lofgof_out
rule beeline_run_sergio_out:
   input: expand('outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.log',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 dataset=['net0'],\
                 algorithm=['DEEPDRIM'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_sergio_out

def get_eval_mem_mb(wildcards):
    switcher = {
        "epr": 16000,
        "auc": 16000,
        "auc3": 16000,
        "pauc4": 16000,
        "auc4": 16000
        }
    return switcher.get(wildcards.metric, 8000)

def get_eval_time(wildcards):
    switcher = {
        "L0": "3:00:00",
        "L1": "3:00:00",
        "L2": "3:00:00",
        "L0_ns": "10:00:00",
        "L1_ns": "10:00:00",
        "L2_ns": "10:00:00",
        "L0_lofgof": "5:00:00",
        "L1_lofgof": "5:00:00",
        "L2_lofgof": "5:00:00",
        "SERGIO_DS4": "1:00:00",
        "SERGIO_DS5": "1:00:00",
        "SERGIO_DS6": "1:00:00",
        "SERGIO_DS7": "1:00:00"
        }
    return switcher.get(wildcards.exp_dir, "03:00:00")

rule beeline_eval:
    input: 'config-files-split/config_{exp_dir}_split/{dataset}/{algorithm}/config.yaml'
    output: 'outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out'
    log: 'outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.log'
    params: overlay="conda_greene/overlay-5GB-200K-beeline20211104.ext3",\
            sif="conda_greene/centos-8.2.2004.sif",\
            D="outputs_eval/{exp_dir}/{dataset}/{algorithm}",\
            output_dir="outputs_eval",\
            output_prefix="{dataset}/{algorithm}/{exp_dir}",\
            jobname="ble_{exp_dir}-{dataset}-{algorithm}",\
            clog_prefix="{exp_dir}-{metric}"
    threads: 1
    resources: mem_mb=get_eval_mem_mb, time=get_eval_time,gpu=""
    shell: """
    mkdir -p {params.D}
    singularity exec \
    --overlay {params.overlay}:ro \
    {params.sif} \
    /bin/bash -c \"source /ext3/env.sh; conda activate BEELINE; \
    python BLEvaluator.py --config {input} --{wildcards.metric} --output_dir {params.output_dir} --output_prefix {params.output_prefix} > {log} 2>&1 \" 
    touch {output}
    """
#; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete    
rule beeline_eval_out:
   input: expand('outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['L0','L1','L2'],\
                 dataset=DATASET_PARAMS['cs']['dataset'],\
                 algorithm=ALGORITHMS,\
                 metric=['auc'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_out
rule beeline_eval_ns_out:
   input: expand('outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['L0_ns','L1_ns','L2_ns'],\
                 dataset=DATASET_PARAMS['ns']['dataset'],\
                 algorithm=ALGORITHMS,\
                 metric=['auc'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_ns_out
rule beeline_eval_ns_out2:
   shell: "echo {rules.beeline_eval_ns_out.input}"
rule beeline_eval_lofgof_out:
   input: expand('outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['L0_lofgof','L1_lofgof','L2_lofgof'],\
                 dataset=DATASET_PARAMS['lofgof']['dataset'],\
                 algorithm=ALGORITHMS,\
                 metric=['auc','auc3'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_lofgof_out
rule beeline_eval_sergio_out:
   input: expand('outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 dataset=['net0'],\
                 algorithm=ALGORITHMS,\
                 metric=['auc','auc3']) 
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_sergio_out
rule beeline_eval_sergio_out2:
   shell: "echo {rules.beeline_eval_sergio_out.input}"
   

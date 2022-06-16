import glob
import os
import sys
import datetime
import pandas as pd
import yaml

BEELINE_NETWORKS_DIR='/scratch/ch153/packages/BEELINE/BEELINE-Networks'
BEELINE_DATA_DIR='/scratch/ch153/packages/BEELINE/BEELINE-data'

DATASETS = ['dream5_1','dream5_3','dream5_4','hESC','hHEP','mDC','mESC','mHSC-E']
DATASETS_NS = ['hESC','hHEP','mDC','mESC','mHSC-E']
ALGORITHMS = ['GENIE3','GRISLI','GRNBOOST2','GRNVBEM',
               'LEAP','PIDC','PPCOR','SCNS','SCODE','SCRIBE',
               'SINCERITIES','SINGE']
DATASET_DIRS = ['L0','L1','L2']
DATASET_DIRS_NS = ['L0_ns','L1_ns','L2_ns']

DATASET_PARAMS = {'dataset':['hESC','hHep','mDC','mESC','mHSC-E'],
                  'netfile':['human/hESC-ChIP-seq-network.csv',
                             'human/HepG2-ChIP-seq-network.csv',
                             'mouse/mDC-ChIP-seq-network.csv',
                             'mouse/mESC-ChIP-seq-network-fixed.csv',
                             'mouse/mHSC-ChIP-seq-network.csv'],
                  'tffile':['human-tfs-upper.csv','human-tfs-upper.csv',
                            'mouse-tfs.csv','mouse-tfs.csv','mouse-tfs.csv']}
DATASET_PARAMS_NS = {'dataset':['hESC','hHep','mDC','mESC','mHSC-E'],
                     'netfile':['human/Non-specific-ChIP-seq-network.csv',
                                'human/Non-specific-ChIP-seq-network.csv',
                                'mouse/Non-Specific-ChIP-seq-network.csv',
                                'mouse/Non-Specific-ChIP-seq-network.csv',
                                'mouse/Non-Specific-ChIP-seq-network.csv'],
                     'tffile':['human-tfs-upper.csv','human-tfs-upper.csv',
                               'mouse-tfs.csv','mouse-tfs.csv','mouse-tfs.csv']}
exp_param_list = []
for exp in ['L0','L0_ns']:
    dataset_df = pd.DataFrame(data=DATASET_PARAMS)
    dataset_df['exp_dir'] = exp
    dataset_df['numgenes'] = 500
    exp_param_list.append(dataset_df)
for exp in ['L1','L1_ns']:
    dataset_df = pd.DataFrame(data=DATASET_PARAMS_NS)
    dataset_df['exp_dir'] = exp
    dataset_df['numgenes'] = 1000
    exp_param_list.append(dataset_df)
EXP_PARAM_DF = pd.concat(exp_param_list)

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
                                          os.path.join(BEELINE_NETWORKS_DIR,'Networks/{u.tffile}'))
    output: expout='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
            netout='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-network.csv'
    log: '{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}.log'
    params: overlay="conda_greene/overlay-5GB-200K-beeline20211104.ext3",\
            sif="conda_greene/centos-8.2.2004.sif",\
            D="{exp_input_dir}/{exp_dir}/{dataset}",\
            p='-p 0.01',c='-c',t='-t',\
            n=lambda wc: get_exp_file(EXP_PARAM_DF,wc.exp_dir,wc.dataset,'{u.numgenes}'),\
            outprefix="{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}",\
            jobname="blg_{exp_dir}-{dataset}",\
            clog_prefix="{exp_dir}"
    threads: 1
#    resources: mem_mb=get_eval_mem_mb, time="00:10:00",
    shell: """
    mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
    singularity exec \
    --overlay {params.overlay}:ro \
    {params.sif} \
    /bin/bash -c \"source /ext3/env.sh; conda activate BEELINE; \
    python generateExpInputs.py -e {input.expfile} -g {input.geneorderingfile} -f {input.netfile} -i {input.tffile} {params.p} {params.c} {params.t} -n={params.n} -o={params.outprefix} > {log} 2>&1 \"
    """
rule beeline_exp_inputs_out:
   input: expand('inputs_beeline/{ds.exp_dir}/{ds.dataset}/{ds.exp_dir}-ExpressionData.csv',\
                 ds=EXP_PARAM_DF.itertuples())

def get_run_mem_mb(wildcards):
    switcher = {
        "L0": 64000,
        "L1": 32000,
        "L2": 64000,
        "L0_ns": 64000,
        "L1_ns": 64000,
        "L2_ns": 64000        
        }
    return switcher.get(wildcards.dataset_dir, 8000)

def get_run_time(wildcards):
    switcher = {
        "L0": "08:00:00",
        "L1": "08:00:00",
        "L2": "08:00:00",
        "L0_ns": "08:00:00",
        "L1_ns": "08:00:00",
        "L2_ns": "08:00:00"
        }
    return switcher.get(wildcards.dataset_dir, "03:00:00")

rule beeline_run:
    input: 'config-files-split/config_{dataset_dir}_split/{dataset}/{algorithm}/config.yaml'
    output: 'outputs/{dataset_dir}/{dataset}/{algorithm}/rankedEdges.log'
#    log: 'outputs/{dataset_dir}/{dataset}/{algorithm}/rankedEdges.log'
    params: D="outputs/{dataset_dir}/{dataset}/{algorithm}",\
            jobname="blr_{dataset_dir}-{dataset}-{algorithm}",\
            clog_prefix="rankedEdges"
    threads: 1
    resources: mem_mb=get_run_mem_mb, time=get_run_time,
    envmodules: "anaconda3/2020.07"
    shell: """
        mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
        python BLRunner.py --config {input} > {output} 2>&1 || true
    """
rule beeline_run_out:
   input: expand('outputs/{dataset_dir}/{dataset}/{algorithm}/rankedEdges.log',\
                 dataset_dir=['L2'],\
                 dataset=DATASETS,\
                 algorithm=ALGORITHMS)
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_out
rule beeline_run_ns_out:
   input: expand('outputs/{dataset_dir}/{dataset}/{algorithm}/rankedEdges.log',\
                 dataset_dir=['L2_ns'],\
                 dataset=DATASETS_NS,\
                 algorithm=ALGORITHMS)
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_ns_out

def get_eval_mem_mb(wildcards):
    switcher = {
        "epr": 16000,
        "auc3": 200000,
        "auc4": 16000
        }
    return switcher.get(wildcards.metric, 8000)

rule beeline_eval:
    input: 'config-files-split/config_{dataset_dir}_split/{dataset}/{algorithm}/config.yaml'
    output: 'outputs_eval/{dataset_dir}/{dataset}/{algorithm}/{dataset_dir}-{metric}.out'
    log: 'outputs_eval/{dataset_dir}/{dataset}/{algorithm}/{dataset_dir}-{metric}.log'
    params: overlay="conda_greene/overlay-5GB-200K-beeline20211104.ext3",\
            sif="conda_greene/centos-8.2.2004.sif",\
            D="outputs_eval/{dataset_dir}/{dataset}/{algorithm}",\
            output_dir="outputs_eval",\
            output_prefix="{dataset}/{algorithm}/{dataset_dir}",\
            jobname="ble_{dataset_dir}-{dataset}-{algorithm}",\
            clog_prefix="{dataset_dir}-{metric}"
    threads: 1
    resources: mem_mb=get_eval_mem_mb, time="00:10:00",
    shell: """
    mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
    singularity exec \
    --overlay {params.overlay}:ro \
    {params.sif} \
    /bin/bash -c \"source /ext3/env.sh; conda activate BEELINE; \
    python BLEvaluator.py --config {input} --{wildcards.metric} --output_dir {params.output_dir} --output_prefix {params.output_prefix} > {log} 2>&1 \" 
    touch {output}
    """
rule beeline_eval_l0_out:
   input: expand('outputs_eval/{dataset_dir}/{dataset}/{algorithm}/{dataset_dir}-{metric}.out',\
                 dataset_dir=['L0'],\
                 dataset=DATASETS,\
                 algorithm=ALGORITHMS,\
                 metric=['epr','auc4'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_l0_out
rule beeline_eval_ns_out:
   input: expand('outputs/{dataset_dir}/{dataset}/{algorithm}/{dataset_dir}-{metric}.log',\
                 dataset_dir=DATASET_DIRS_NS,\
                 dataset=DATASETS_NS,\
                 algorithm=ALGORITHMS,\
                 metric=['epr','auc4']) 
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_ns_out


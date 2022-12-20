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

ALGORITHMS = ['CICT','DEEPDRIM','DEEPDRIM5','DEEPDRIM6','DEEPDRIM7','GENIE3','GRISLI',
              'GRNBOOST2','GRNVBEM',
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
EXP_LEVELS = {'L0':{'num_genes':500,'num_rand_genes':0,'num_exclude_genes':500},
              'L1':{'num_genes':1000,'num_rand_genes':0,'num_exclude_genes':1000},
              'L2':{'num_genes':500,'num_rand_genes':500,'num_exclude_genes':1000}
              }
for level_name,level_settings in EXP_LEVELS.items():
    for gt_name,gt_settings in DATASET_PARAMS.items():
        dataset_df = pd.DataFrame(data=gt_settings)
        exp_dir = level_name if gt_name=='cs' else level_name+'_'+gt_name
        dataset_df['gt_name'] = gt_name
        dataset_df['exp_dir'] = exp_dir
        dataset_df['num_genes'] = level_settings['num_genes']
        dataset_df['num_rand_genes'] = level_settings['num_rand_genes']
        dataset_df['num_exclude_genes'] = level_settings['num_exclude_genes']
        exp_param_list.append(dataset_df)
EXP_PARAM_DF = pd.concat(exp_param_list)
print(EXP_PARAM_DF)

def get_exp_file(exp_param_df,exp_dir,dataset,upath):
    ds_df = exp_param_df[(exp_param_df.exp_dir==exp_dir) & (exp_param_df.dataset==dataset)].iloc[[0]]
    #print(ds_df)
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
            x=lambda wc: get_exp_file(EXP_PARAM_DF,wc.exp_dir,wc.dataset,'{u.num_exclude_genes}'),\
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
                       python generateExpInputs.py -e {input.expfile} -g {input.geneorderingfile} -f {input.netfile} -i {input.tffile} {params.p} {params.c} {params.t} -n {params.n} -r {params.r} -x {params.x} -o {params.outprefix} > {log} 2>&1 \"
    cp {input.ptfile} {output.ptout}
    rm -f {params.D}/ExpressionData.csv; ln -rs {output.expout} {params.D}/ExpressionData.csv
    rm -f {params.D}/refNetwork.csv; ln -rs {output.netout} {params.D}/refNetwork.csv
    rm -f {params.D}/PseudoTime.csv; ln -rs {output.ptout} {params.D}/PseudoTime.csv
    """
rule beeline_exp_inputs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/{ds.exp_dir}-ExpressionData.csv',\
                 ds=EXP_PARAM_DF.itertuples())

rule beeline_exp_tfs:
    input: expfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
           netfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-network.csv',\
           tffile=lambda wc: get_exp_file(EXP_PARAM_DF,wc.exp_dir,wc.dataset,\
                                          os.path.join(BEELINE_NETWORKS_DIR,'Networks/{u.tffile}'))
    output: exptf1='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData-TFs.csv',\
            exptf2='{exp_input_dir}/{exp_dir}/{dataset}/ExpressionData-TFs.csv',\
            nettf1='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-network-TFs.csv',\
            nettf2='{exp_input_dir}/{exp_dir}/{dataset}/refNetwork-TFs.csv'
    params: D="{exp_input_dir}/{exp_dir}/{dataset}"
    threads: 1
    shell: """
        awk -v FS=',' -v OFS=',' '(FNR>1) && (FNR==NR) {{ a[$1]; next }} $1 in a {{ print }}' {input.expfile} {input.tffile} > {output.exptf1}
        awk -v FS=',' -v OFS=',' '(FNR>1) && (FNR==NR) {{ a[$1]; next }} $1 in a {{ print }}' {input.netfile} {input.tffile} > {output.nettf1}
        ln -rs {output.exptf1} {output.exptf2}
        ln -rs {output.nettf1} {output.nettf2}
    """
rule beeline_exp_tfs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/{ds.exp_dir}-ExpressionData-TFs.csv',\
                 ds=EXP_PARAM_DF.itertuples())

# Generate pairs data for DEEPDRIM using L0/L1/L2 CICT training datasets
rule beeline_exp_cictpairs:
   input: netfile_train='outputs/{exp_dir}/{dataset}/CICT/train.csv',\
          netfile_test='outputs/{exp_dir}/{dataset}/CICT/test.csv'
   output: netout_train='{exp_input_dir}/{exp_dir}/{dataset}/cictTrain.csv',\
           netout_test='{exp_input_dir}/{exp_dir}/{dataset}/cictTest.csv'
   shell: """
        cp {input.netfile_train} {output.netout_train}
        cp {input.netfile_test} {output.netout_test}
   """
rule beeline_exp_cictpairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/cict{train_test}.csv',\
                 train_test=['Train','Test'],\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())

# rsync -Pavz /scratch/as15096/eric/outputs/L2_lofgof/mESC/mESC_lofgof_training_sets/ /scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs_cict_learn/L2_lofgof/mESC
rule beeline_exp_cictpairs_multi:
   input: netfile_train='outputs_cict_learn/{exp_dir}/{dataset}/{train_i}/training.csv',\
          netfile_test='outputs_cict_learn/{exp_dir}/{dataset}/{train_i}/test.csv'
   output: netout_train='{exp_input_dir}/{exp_dir}/{dataset}/CICT/run_{train_i}/train.csv',\
           netout_test='{exp_input_dir}/{exp_dir}/{dataset}/CICT/run_{train_i}/test.csv'
   shell: """
        awk -v FS=',' -v OFS='\\t' '{{ if (NR>1) {{ $1=toupper($1); $2=toupper($2) }} print $1,$2,$3 }}' {input.netfile_train} > {output.netout_train}
        awk -v FS=',' -v OFS='\\t' '{{ if (NR>1) {{ $1=toupper($1); $2=toupper($2) }} print $1,$2,$3 }}' {input.netfile_test} > {output.netout_test}
   """
rule beeline_exp_cictpairs_multi_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/CICT/run_{train_i}/{train_test}.csv',\
                 train_test=['train','test'],\
                 train_i=range(1,11),\
                 ds=EXP_PARAM_DF[(EXP_PARAM_DF.exp_dir=='L2_lofgof') & (EXP_PARAM_DF.dataset=='mESC')].groupby(['exp_dir','dataset']).count().reset_index().itertuples())


# Generate training pairs data for DEEPDRIM using full BEELINE datasets
rule beeline_exp_deepdrim_pairs:
   input: expfile=os.path.join(BEELINE_DATA_DIR,'inputs/scRNA-Seq/{dataset}/ExpressionData-upper.csv'),\
          netfile=lambda wc: get_exp_file(EXP_PARAM_DF,\
                                          'L0' if wc.gt_name=='cs' else 'L0'+'_'+wc.gt_name,\
                                          wc.dataset,\
                                          os.path.join(BEELINE_NETWORKS_DIR,'Networks/{u.netfile}'))
   output: expout="{exp_input_dir}/ALL_{gt_name}/{dataset}/DEEPDRIM/ExpressionData.csv",\
           netout="{exp_input_dir}/ALL_{gt_name}/{dataset}/DEEPDRIM/PositivePairs.txt",\
           pairout="{exp_input_dir}/ALL_{gt_name}/{dataset}/DEEPDRIM/training_pairsDEEPDRIM.txt",\
           cvout="{exp_input_dir}/ALL_{gt_name}/{dataset}/DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/ALL_{gt_name}/{dataset}/DEEPDRIM/training_pairsDEEPDRIM.log"
   params: overlay="images_singularity_overlay/DEEPDRIM.ext3",\
           sif="images_singularity/DEEPDRIM.sif",\
           D="{exp_input_dir}/ALL_{gt_name}/{dataset}/DEEPDRIM",\
           label='DEEPDRIM',\
           jobname='blg_{dataset}_deepdrim',\
           clog_prefix='training_pairsDEEPDRIM'
   threads: 1
   resources: mem_mb='16000', time='08:00:00', gpu=''
   shell: """
     mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
     cp {input.expfile} {output.expout}
     tail -n +2 {input.netfile} > {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /ext3/DeepDRIM/generate_pairs_realdata.py -expr_file ExpressionData.csv -pos_pair_file PositivePairs.txt -label {params.label} > training_pairsDEEPDRIM.log 2>&1; \
                       python /ext3/DeepDRIM/generate_cvfold_realdata.py -TF_divide_pos_file training_pairsDEEPDRIM.txtTF_divide_pos.txt -cross_validation_fold_divide_file training_pairsDEEPDRIM.txtCV_fold_divide.txt >> training_pairsDEEPDRIM.log 2>&1 \"
   """
rule beeline_exp_deepdrim_pairs_out:
   input: expand('inputs_beeline2/ALL_{ds.gt_name}/{ds.dataset}/DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide.txt',\
                 ds=EXP_PARAM_DF.groupby(['gt_name','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim_pairs_out

# Generate training pairs data for DEEPDRIM using L0/L1/L2 datasets
rule beeline_exp_deepdrim5_pairs:
   input: expfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
          netfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-network.csv'
   output: expout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM5/ExpressionData.csv",\
           netout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM5/PositivePairs.txt",\
           pairout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM5/training_pairsDEEPDRIM5.txt",\
           cvout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM5/training_pairsDEEPDRIM5.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM5/training_pairsDEEPDRIM5.log"
   params: overlay="images_singularity_overlay/DEEPDRIM5.ext3",\
           sif="images_singularity/DEEPDRIM5.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM5",\
           label='DEEPDRIM5',\
           jobname='blg_{dataset}_deepdrim5',\
           clog_prefix='training_pairsDEEPDRIM5'
   threads: 1
   resources: mem_mb='16000', time='02:00:00', gpu=''
   shell: """
     mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
     cp {input.expfile} {output.expout}
     tail -n +2 {input.netfile} > {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /ext3/DeepDRIM/generate_pairs_realdata.py -expr_file ExpressionData.csv -pos_pair_file PositivePairs.txt -label {params.label} > training_pairsDEEPDRIM5.log 2>&1; \
                       python /ext3/DeepDRIM/generate_cvfold_realdata.py -TF_divide_pos_file training_pairsDEEPDRIM5.txtTF_divide_pos.txt -cross_validation_fold_divide_file training_pairsDEEPDRIM5.txtCV_fold_divide.txt >> training_pairsDEEPDRIM5.log 2>&1 \"
   """
rule beeline_exp_deepdrim5_pairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM5/training_pairsDEEPDRIM5.txtCV_fold_divide.txt',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim5_pairs_out   

# Generate training pairs data for DEEPDRIM using L0/L1/L2 CICT training datasets
rule beeline_exp_deepdrim6_pairs:
   input: expfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
          netfile="{exp_input_dir}/{exp_dir}/{dataset}/cictTrain.csv"
   output: expout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM6/ExpressionData.csv",\
           netout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM6/cictTrain.csv",\
           pairout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM6/training_pairsDEEPDRIM6.txt",\
           cvout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM6/training_pairsDEEPDRIM6.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM6/training_pairsDEEPDRIM6.log"
   params: overlay="images_singularity_overlay/DEEPDRIM6.ext3",\
           sif="images_singularity/DEEPDRIM6.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM6",\
           label='DEEPDRIM6',\
           jobname='blg_{dataset}_deepdrim6',\
           clog_prefix='training_pairsDEEPDRIM6'
   threads: 1
   resources: mem_mb='16000', time='02:00:00', gpu=''
   shell: """
     mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
     cp {input.expfile} {output.expout}
     cp {input.netfile} {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pairs_cictpairs.py -expr_file ExpressionData.csv -cict_pair_file cictTrain.csv -label {params.label} > training_pairsDEEPDRIM6.log 2>&1 
                       python /ext3/DeepDRIM/generate_cvfold_realdata.py -TF_divide_pos_file training_pairsDEEPDRIM6.txtTF_divide_pos.txt -cross_validation_fold_divide_file training_pairsDEEPDRIM6.txtCV_fold_divide.txt >> training_pairsDEEPDRIM6.log 2>&1 \"
   """
rule beeline_exp_deepdrim6_pairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM6/training_pairsDEEPDRIM6.txt',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim6_pairs_out

# Generate training pairs data for DEEPDRIM7 using L0/L1/L2 CICT training+testing datasets
rule beeline_exp_deepdrim7_pairs:
   input: expfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
          netfile_cicttrain="{exp_input_dir}/{exp_dir}/{dataset}/cictTrain.csv",\
          netfile_cicttest="{exp_input_dir}/{exp_dir}/{dataset}/cictTest.csv"
   output: expout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM7/ExpressionData.csv",\
           netout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM7/cictLearn.csv",\
           pairout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM7/training_pairsDEEPDRIM7.txt",\
           cvout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM7/training_pairsDEEPDRIM7.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM7/training_pairsDEEPDRIM7.log"
   params: overlay="images_singularity_overlay/DEEPDRIM7.ext3",\
           sif="images_singularity/DEEPDRIM7.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM7",\
           label='DEEPDRIM7',random_state='12',\
           jobname='blg_{dataset}_deepdrim7',\
           clog_prefix='training_pairsDEEPDRIM7'
   threads: 1
   resources: mem_mb='16000', time='02:00:00', gpu=''
   shell: """
     mkdir -p {params.D}
     cp {input.expfile} {output.expout}
     cp {input.netfile_cicttrain} {output.netout}
     tail -n +2 {input.netfile_cicttest} >> {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pairs_cictpairs2.py -expr_file ExpressionData.csv -cict_pair_file cictLearn.csv -label {params.label} -random_state {params.random_state} > training_pairsDEEPDRIM7.log 2>&1 \"
   """
rule beeline_exp_deepdrim7_pairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM7/training_pairsDEEPDRIM7.txt',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim7_pairs_out


# Generate training pairs data for DEEPDRIM8 using multiple sets of L2_lofgof CICT training+testing sets
rule beeline_exp_deepdrim8_pairs:
   input: expfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
          netfile_cicttrain="{exp_input_dir}/{exp_dir}/{dataset}/CICT/run_{train_i}/train.csv",\
          netfile_cicttest="{exp_input_dir}/{exp_dir}/{dataset}/CICT/run_{train_i}/test.csv"
   output: expout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM8/run_{train_i}/ExpressionData.csv",\
           netout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM8/run_{train_i}/cictLearn.csv",\
           pairout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM8/run_{train_i}/training_pairsDEEPDRIM8.txt",\
           cvout="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM8/run_{train_i}/training_pairsDEEPDRIM8.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM8/run_{train_i}/training_pairsDEEPDRIM8.log"
   params: overlay="images_singularity_overlay/DEEPDRIM8.ext3",\
           sif="images_singularity/DEEPDRIM8.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/DEEPDRIM8/run_{train_i}",\
           label='DEEPDRIM8',random_state='12',\
           jobname='blg_{dataset}_deepdrim8',\
           clog_prefix='training_pairsDEEPDRIM8'
   threads: 1
   resources: mem_mb='16000', time='02:00:00', gpu=''
   shell: """
     mkdir -p {params.D}
     cp {input.expfile} {output.expout}
     cp {input.netfile_cicttrain} {output.netout}
     tail -n +2 {input.netfile_cicttest} >> {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pairs_cictpairs2.py -expr_file ExpressionData.csv -cict_pair_file cictLearn.csv -label {params.label} -random_state {params.random_state} > training_pairsDEEPDRIM8.log 2>&1 \"
   """
rule beeline_exp_deepdrim8_pairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM8/run_{train_i}/training_pairsDEEPDRIM8.txt',\
                 ds=EXP_PARAM_DF[(EXP_PARAM_DF.exp_dir=='L2_lofgof') & (EXP_PARAM_DF.dataset=='mESC')].groupby(['exp_dir','dataset']).count().reset_index().itertuples(),\
                 train_i=range(1,11))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim8_pairs_out

# Generate prediction pairs data for DEEPDRIM,DEEPDRIM5,DEEPDRIM6 using outgoing edges for TFs in the expression data
rule beeline_exp_deepdrimXpred_pairs:
   input: expfile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData.csv',\
          tffile='{exp_input_dir}/{exp_dir}/{dataset}/{exp_dir}-ExpressionData-TFs.csv',\
          config='config-files-split/config_{exp_dir}_split/{dataset}/{dd_version}/config.yaml'
   output: pairout="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEPDRIM(5|6)?}/predict_pairs{dd_version}.txt",\
           tfout="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEPDRIM(5|6)?}/predict_pairs{dd_version}.txtTF_divide_pos.txt",\
           nameout="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEPDRIM(5|6)?}/predict_pairs{dd_version}.txtGeneName_map.txt"
   log: "{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/predict_pairs{dd_version}.log"
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}",\
           label='{dd_version}',\
           tf_name_col=0,
           jobname='blg_{dataset}_{dd_version}',\
           clog_prefix='predict_pairs{dd_version}'
   threads: 1
   resources: mem_mb='8000', time='00:10:00', gpu=''
   shell: """
     mkdir -p {params.D}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pred_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pred_pairs_tf.py -expr_file ExpressionData.csv -tf_file {input.tffile} -tf_name_col {params.tf_name_col} -label {params.label} > predict_pairs{wildcards.dd_version}.log 2>&1 \"
   """
rule beeline_exp_deepdrimXpred_pairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/{dd_version}/predict_pairs{dd_version}.txt',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples(),\
                 dd_version=['DEEPDRIM','DEEPDRIM5','DEEPDRIM6'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim7pred_pairs_out

# Generate prediction pairs data for DEEPDRIM7,DEEPDRIM8 using all pairwise interaction bwtween genes in the expression data
rule beeline_exp_deepdrimYpred_pairs:
   input: expfile='{exp_input_dir}/{exp_dir}/{dataset}/ExpressionData.csv',\
          config='config-files-split/config_{exp_dir}_split/{dataset}/{dd_version}/config.yaml'
   output: expout='{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEDRPIM(7|8)}/ExpressionData.csv',\
           pairout="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEPDRIM(7|8)}/predict_pairs{dd_version}.txt",\
           tfout="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEPDRIM(7|8)}/predict_pairs{dd_version}.txtTF_divide_pos.txt",\
           nameout="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEPDRIM(7|8)}/predict_pairs{dd_version}.txtGeneName_map.txt"
   log: "{exp_input_dir}/{exp_dir}/{dataset}/{dd_version,DEEPDRIM(7|8)}/predict_pairs{dd_version}.log"
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}",\
           label='{dd_version}',\
           num_batches=lambda wildcards, input: yaml.safe_load(open(input.config,'r'))['input_settings']['algorithms'][0]['params']['numBatches'][0],\
           jobname='blg_{dataset}_{dd_version}',\
           clog_prefix='predict_pairs{dd_version}'
   threads: 1
   resources: mem_mb='8000', time='00:10:00', gpu=''
   shell: """
     mkdir -p {params.D}
     cp {input.expfile} {output.expout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pred_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pred_pairs_all.py -expr_file ExpressionData.csv -num_batches {params.num_batches} -label {params.label} > {params.clog_prefix}.log 2>&1 \"
   """
rule beeline_exp_deepdrim7pred_pairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM7/predict_pairsDEEPDRIM7.txt',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim7pred_pairs_out
rule beeline_exp_deepdrim8pred_pairs_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM8/predict_pairsDEEPDRIM8.txt',\
                 ds=EXP_PARAM_DF[(EXP_PARAM_DF.exp_dir=='L2_lofgof') & (EXP_PARAM_DF.dataset=='mESC')].groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim8pred_pairs_out

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
    output: netout='{exp_input_dir}/SERGIO_{exp_name,DS[1-3]}/refNetwork.csv',\
            clout='{exp_input_dir}/SERGIO_{exp_name,DS[1-3]}/ClusterIds.csv',\
            ptout='{exp_input_dir}/SERGIO_{exp_name,DS[1-3]}/PseudoTime.csv'
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

rule beeline_sergio_bulk_dsdata:# data common for all simulated files in each dataset
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
    output: expout_csv='{exp_input_dir}/SERGIO_{exp_name,DS[1-3]}/net{network_i}/ExpressionData.csv',\
            expoutT_csv='{exp_input_dir}/SERGIO_{exp_name,DS[1-3]}/net{network_i}/ExpressionDataT.csv',\
            netout='{exp_input_dir}/SERGIO_{exp_name,DS[1-3]}/net{network_i}/refNetwork.csv',\
            ptout='{exp_input_dir}/SERGIO_{exp_name,DS[1-3]}/net{network_i}/PseudoTime.csv'
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
                 ds=SERGIO_PARAMS_DF.itertuples(),network_i=range(15))

rule beeline_sergio_tfs:
    input: expfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/ExpressionData.csv',\
           netfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/refNetwork.csv'
    output: exptf1='{exp_input_dir}/SERGIO_{exp_name,DS[4-7]}/net{network_i}/ExpressionData-TFs.csv',\
            nettf1='{exp_input_dir}/SERGIO_{exp_name,DS[4-7]}/net{network_i}/refNetwork-TFs.csv'
    params: D="{exp_input_dir}/SERGIO_{exp_name,DS[4-7]}/net{network_i}"
    threads: 1
    shell: """
        awk -v FS=',' -v OFS=',' '{{ if (NR>1) {{ print $1 }} }}' {input.expfile} | sort -u > {output.exptf1}
        awk -v FS=',' -v OFS=',' '{{ if (NR>1) {{ print $1 }} }}' {input.netfile} | sort -u > {output.nettf1}
    """
rule beeline_sergio_tfs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/ExpressionData-TFs.csv',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
   
rule beeline_sergio_bulk_inputs:
    input: expfileT='{exp_input_dir}/SERGIO_{dataset_name}/net{network_i}/ExpressionDataT.csv'
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

   
# Generate pairs data for DEEPDRIM using SERGIO CICT training datasets
# should be using beeline_exp_cictpairs
rule beeline_sergio_cictpairs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/cict{train_test}.csv',\
                 train_test=['Train','Test'],\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
   
# generate pairs data for DEEPDRIM using SERGIO full datasets 
rule beeline_sergio_deepdrim_pairs:
   input: expfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/ExpressionData.csv',\
          netfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/refNetwork.csv'
   output: expout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM[5]?}/ExpressionData.csv',\
           netout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,,DEEPDRIM[5]?}/PositivePairs.txt',\
           pairout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM[5]?}/training_pairs{dd_version}.txt',\
           cvout="{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM[5]?}/training_pairs{dd_version}.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version}/training_pairs{dd_version}.log"
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version}",\
           label='{dd_version}',\
           jobname='blg_{exp_name}_net{network_i}_deepdrim',\
           clog_prefix='training_pairs{dd_version}'
   threads: 1
   resources: mem_mb='16000', time='00:15:00', gpu=''
   shell: """
     mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
     cp {input.expfile} {output.expout}
     tail -n +2 {input.netfile} | cut -d, -f 1-2 > {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /ext3/DeepDRIM/generate_pairs_realdata.py -expr_file ExpressionData.csv -pos_pair_file PositivePairs.txt -label {params.label} > training_pairs{wildcards.dd_version}.log 2>&1; \
                       python /ext3/DeepDRIM/generate_cvfold_realdata.py -TF_divide_pos_file training_pairs{wildcards.dd_version}.txtTF_divide_pos.txt -cross_validation_fold_divide_file training_pairs{wildcards.dd_version}.txtCV_fold_divide.txt >> training_pairs{wildcards.dd_version}.log 2>&1 \"
   """
rule beeline_sergio_deepdrim_pairs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM/training_pairsDEEPDRIM.txt',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim_pairs_out
rule beeline_sergio_deepdrim5_pairs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM5/training_pairsDEEPDRIM5.txt',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim5_pairs_out

# Generate pairs data for DEEPDRIM using CICT training datasets
rule beeline_sergio_deepdrim6_pairs:
   input: expfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/ExpressionData.csv',\
          netfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/cictTrain.csv'
   output: expout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM6}/ExpressionData.csv',\
           netout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,,DEEPDRIM6}/cictTrain.csv',\
           pairout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM6}/training_pairs{dd_version}.txt',\
           cvout="{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM6}/training_pairs{dd_version}.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version,DEEPDRIM6}/training_pairs{dd_version}.log"
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version}",\
           label='{dd_version}',\
           jobname='blg_{exp_name}_net{network_i}_{dd_version}',\
           clog_prefix='training_pairs{dd_version}'
   threads: 1
   resources: mem_mb='16000', time='00:15:00', gpu=''
   shell: """
     mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
     cp {input.expfile} {output.expout}
     cp {input.netfile} {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pairs_cictpairs.py -expr_file ExpressionData.csv -cict_pair_file cictTrain.csv -label {params.label} > training_pairs{wildcards.dd_version}.log 2>&1
                       python /ext3/DeepDRIM/generate_cvfold_realdata.py -TF_divide_pos_file training_pairs{wildcards.dd_version}.txtTF_divide_pos.txt -cross_validation_fold_divide_file training_pairs{wildcards.dd_version}.txtCV_fold_divide.txt >> training_pairs{wildcards.dd_version}.log 2>&1 \"
   """
rule beeline_sergio_deepdrim6_pairs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM6/training_pairsDEEPDRIM6.txt',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim6_pairs_out

# Generate pairs data for DEEPDRIM7 using SERGIO CICT training+testing datasets
rule beeline_sergio_deepdrim7_pairs:
   input: expfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/ExpressionData.csv',\
          netfile_cicttrain='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/cictTrain.csv',\
          netfile_cicttest='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/cictTest.csv'
   output: expout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM7}/ExpressionData.csv',\
           netout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,,DEEPDRIM7}/cictLearn.csv',\
           pairout='{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM7}/training_pairs{dd_version}.txt',\
           cvout="{exp_input_dir}/SERGIO_{exp_name,[^_]+}/net{network_i}/{dd_version,DEEPDRIM7}/training_pairs{dd_version}.txtCV_fold_divide.txt"
   log: "{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version,DEEPDRIM7}/training_pairs{dd_version}.log"
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version}",\
           label='{dd_version}',\
           jobname='blg_{exp_name}_net{network_i}_{dd_version}',\
           clog_prefix='training_pairs{dd_version}'
   threads: 1
   resources: mem_mb='16000', time='00:15:00', gpu=''
   shell: """
     mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
     cp {input.expfile} {output.expout}
     cp {input.netfile_cicttrain} {output.netout}
     tail -n +2 {input.netfile_cicttest} >> {output.netout}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pairs_cictpairs2.py -expr_file ExpressionData.csv -cict_pair_file cictLearn.csv -label {params.label} > training_pairs{wildcards.dd_version}.log 2>&1 \"
   """
rule beeline_sergio_deepdrim7_pairs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM7/training_pairsDEEPDRIM7.txt',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim7_pairs_out


# Generate prediction pairs data for DEEPDRIM,DEEPDRIM5,DEEPDRIM6 using outgoing edges for TFs in the expression data
rule beeline_sergio_deepdrimXpred_pairs:
   input: expfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/ExpressionData.csv',\
          tffile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/ExpressionData-TFs.csv',\
          config='config-files-split/config_SERGIO_{exp_name}_split/net{network_i}/{dd_version}/config.yaml'
   output: pairout="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version,DEEPDRIM(5|6)?}/predict_pairs{dd_version}.txt",\
           tfout="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version,DEEPDRIM(5|6)?}/predict_pairs{dd_version}.txtTF_divide_pos.txt",\
           nameout="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version,DEEPDRIM(5|6)?}/predict_pairs{dd_version}.txtGeneName_map.txt"
   log: "{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version}/predict_pairs{dd_version}.log"
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/{dd_version}",\
           label='{dd_version}',\
           tf_name_col=0,
           jobname='blg_net{network_i}_{dd_version}',\
           clog_prefix='predict_pairs{dd_version}'
   threads: 1
   resources: mem_mb='8000', time='00:10:00', gpu=''
   shell: """
     mkdir -p {params.D}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pred_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pred_pairs_tf.py -expr_file ExpressionData.csv -tf_file {input.tffile} -tf_name_col {params.tf_name_col} -label {params.label} > predict_pairs{wildcards.dd_version}.log 2>&1 \"
   """
rule beeline_sergio_deepdrimXpred_pairs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/{dd_version}/predict_pairs{dd_version}.txt',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)],\
                 dd_version=['DEEPDRIM5'])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrimXpred_pairs_out

rule beeline_sergio_deepdrim7pred_pairs:
   input: expfile='{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/ExpressionData.csv',\
          config='config-files-split/config_SERGIO_{exp_name}_split/net{network_i}/DEEPDRIM7/config.yaml'
   output: pairout="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/DEEPDRIM7/predict_pairsDEEPDRIM7.txt",\
           tfout="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/DEEPDRIM7/predict_pairsDEEPDRIM7.txtTF_divide_pos.txt",\
           nameout="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/DEEPDRIM7/predict_pairsDEEPDRIM7.txtGeneName_map.txt"
   log: "{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/DEEPDRIM7/predict_pairsDEEPDRIM7.log"
   params: overlay="images_singularity_overlay/DEEPDRIM7.ext3",\
           sif="images_singularity/DEEPDRIM7.sif",\
           D="{exp_input_dir}/SERGIO_{exp_name}/net{network_i}/DEEPDRIM7",\
           num_batches=lambda wildcards, input: yaml.safe_load(open(input.config,'r'))['input_settings']['algorithms'][0]['params']['numBatches'][0],\
           label='DEEPDRIM7',\
           jobname='blg_{exp_name}_net{network_i}_DEEPDRIM7',\
           clog_prefix='predict_pairsDEEPDRIM7'
   threads: 1
   resources: mem_mb='16000', time='00:15:00', gpu=''
   shell: """
     mkdir -p {params.D}
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_pred_pairs.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_pred_pairs_all.py -expr_file ExpressionData.csv -num_batches {params.num_batches} -label {params.label} > predict_pairsDEEPDRIM7.log 2>&1 \"
   """
rule beeline_sergio_deepdrim7pred_pairs_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM7/predict_pairsDEEPDRIM7.txt',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim7_pairs_out

ruleorder: beeline_sergio_deepdrim_pairs > beeline_exp_deepdrim5_pairs
ruleorder: beeline_sergio_deepdrim6_pairs > beeline_exp_deepdrim6_pairs
ruleorder: beeline_sergio_deepdrim7_pairs > beeline_exp_deepdrim7_pairs
ruleorder: beeline_sergio_deepdrim7pred_pairs > beeline_exp_deepdrimXpred_pairs 

def get_tf_num(wildcards,tfdiv_file):
    with open(tfdiv_file.format(**wildcards),'r') as fp:
        tf_num = len(fp.readlines()) - 1
    return tf_num

# generate representation of training datasets for DEEPDRIM, DEEPDRIM5, DEEPDRIM6, DEEPDRIM7, DEEPDRIM8
rule beeline_deepdrim_repr:
   input: expr_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i}/ExpressionData.csv",\
          pairs_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i}/training_pairs{dd_version}.txt",\
          tfdiv_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i}/training_pairs{dd_version}.txtTF_divide_pos.txt",\
          genename_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i}/{dd_version}_geneName_map.txt"
   output: '{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i,(run_[0-9]+)?}/representation_train.out'
   log: '{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i,(run_[0-9]+)?}/representation_train.log'
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i}",\
           load_from_h5=False,load_split_batch_pos=True,tf_order_random=False,\
           tf_num=lambda wc: get_tf_num(wc,"{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/{run_i}/training_pairs{dd_version}.txtTF_divide_pos.txt"),\
           jobname='blr_{dataset}_{dd_version}',\
           clog_prefix='representation_train'
   threads: 1
   resources: mem_mb='16000', time='2:00:00', gpu=''
   shell: """
     mkdir -p {params.D}; rm -rf {params.D}/representation_train
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_repr.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_input_realdata.py -out_dir representation_train -expr_file ExpressionData.csv -pairs_for_predict_file training_pairs{wildcards.dd_version}.txt -geneName_map_file {wildcards.dd_version}_geneName_map.txt -flag_load_from_h5 {params.load_from_h5} -flag_load_split_batch_pos {params.load_split_batch_pos} -TF_divide_pos_file training_pairs{wildcards.dd_version}.txtTF_divide_pos.txt -TF_num {params.tf_num} -TF_order_random {params.tf_order_random} > {params.clog_prefix}.log 2>&1 \"
     echo "Complete" > {output}
   """

rule beeline_deepdrim_repr_out:
   input: expand('inputs_beeline2/ALL_{ds.gt_name}/{ds.dataset}/DEEPDRIM/representation_train/version11/0_xdata.npy',\
                 ds=EXP_PARAM_DF.groupby(['gt_name','dataset']).count().reset_index().itertuples()) +\
        expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM/representation_train/version11/0_xdata.npy',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_deepdrim_repr_out
rule beeline_deepdrim5_repr_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM5/representation_train/version11/0_xdata.npy',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples()) +\
          expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM5/representation_train/version11/0_xdata.npy',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_deepdrim5_repr_out
rule beeline_deepdrim6_repr_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM6/representation_train/version11/0_xdata.npy',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples()) +\
          expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM6/representation_train/version11/0_xdata.npy',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_deepdrim6_repr_out
rule beeline_deepdrim7_repr_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM7/representation_train.out',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples()) +\
          expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM7/representation_train.out',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_deepdrim7_repr_out
rule beeline_deepdrim8_repr_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM8/run_{train_i}/representation_train.out',\
                 ds=EXP_PARAM_DF[(EXP_PARAM_DF.exp_dir=='L2_lofgof') & (EXP_PARAM_DF.dataset=='mESC')].groupby(['exp_dir','dataset']).count().reset_index().itertuples(),\
                 train_i=range(1,11))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_deepdrim8_repr_out

# generate representation of prediction datasets for DEEPDRIM7/8
rule beeline_deepdrim_pred_repr:
   input: expr_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/ExpressionData.csv",\
          pairs_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/predict_pairs{dd_version}.txt",\
          tfdiv_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/predict_pairs{dd_version}.txtTF_divide_pos.txt",\
          genename_file="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/predict_pairs{dd_version}.txtGeneName_map.txt"
   output: '{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/representation_predict.out'
   log: '{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/representation_predict.log'
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}",\
           tf_num=lambda wc: get_tf_num(wc,"{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/predict_pairs{dd_version}.txtTF_divide_pos.txt"),\
           load_from_h5=False,load_split_batch_pos=True,tf_order_random=False,\
           jobname='blg_{dataset}_{dd_version}',\
           clog_prefix='representation_predict'
   threads: 1
   resources: mem_mb='64000', time='24:00:00', gpu=''
   shell: """
     mkdir -p {params.D}; rm -rf {params.D}/representation_predict
     singularity exec \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                       command time -v -o time_generate_predict.txt python /scratch/ch153/packages/DeepDRIM/hlab1/DeepDRIM/generate_input_realdata.py -out_dir representation_predict -expr_file ExpressionData.csv -pairs_for_predict_file predict_pairs{wildcards.dd_version}.txt -geneName_map_file predict_pairs{wildcards.dd_version}.txtGeneName_map.txt -flag_load_from_h5 {params.load_from_h5} -flag_load_split_batch_pos {params.load_split_batch_pos} -TF_divide_pos_file predict_pairs{wildcards.dd_version}.txtTF_divide_pos.txt -TF_num {params.tf_num} -TF_order_random {params.tf_order_random} > {params.clog_prefix}.log 2>&1 \"
     echo "Complete" > {output}
   """
rule beeline_exp_deepdrim7pred_repr_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM7/representation_predict.out',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim7pred_repr_out

rule beeline_exp_deepdrim8pred_repr_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM8/representation_predict.out',\
                 ds=EXP_PARAM_DF[(EXP_PARAM_DF.exp_dir=='L2_lofgof') & (EXP_PARAM_DF.dataset=='mESC')].groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim8pred_repr_out

rule beeline_sergio_deepdrim7pred_repr_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM7/representation_predict.out',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim7pred_repr_out


### NEED TO REMOVE START: DEEPDRIM TRAINING is now run as part of BLRun####
def get_deepdrim_train_mem_mb(wildcards):
    switcher = {
        "ALL_cs": 280000,
        "ALL_ns": 350000,
        "ALL_lofgof": 48000,
        "L0": 64000,
        "L1": 64000,
        "L2": 64000,
        "L0_ns": 64000,
        "L1_ns": 64000,
        "L2_ns": 64000,
        "L0_lofgof": 128000,
        "L1_lofgof": 128000,
        "L2_lofgof": 128000,        
        "SERGIO_DS4":16000,
        "SERGIO_DS5":16000,
        "SERGIO_DS6":16000,
        "SERGIO_DS7":16000
        }
    return switcher.get(wildcards.exp_dir, 8000)

def get_deepdrim_train_time(wildcards):
    switcher = {
        "ALL_cs": "12:00:00",
        "ALL_ns": "8:00:00",
        "ALL_lofgof": "4:00:00",
        "L0": "2:00:00",
        "L1": "2:00:00",
        "L2": "2:00:00",
        "L0_ns": "2:00:00",
        "L1_ns": "2:00:00",
        "L2_ns": "2:00:00",
        "L0_lofgof": "12:00:00",
        "L1_lofgof": "12:00:00",
        "L2_lofgof": "12:00:00",
        "SERGIO_DS4": "2:00:00",
        "SERGIO_DS5": "2:00:00",
        "SERGIO_DS6": "2:00:00",
        "SERGIO_DS7": "2:00:00"
        }
    return switcher.get(wildcards.exp_dir, "03:00:00")

# train DEEPDRIM (using full datasets), DEEPDRIM5 (L0/L1/L2_* full datasets), DEEPDRIM6 (L0/L1/L2_* CICT sample datesets)
rule beeline_deepdrim_train:
   input: repr_file='{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/representation_train/version11/0_xdata.npy',\
          cv_file='{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/training_pairs{dd_version}.txtCV_fold_divide.txt'
   output: '{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/training_pairs{dd_version}.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5'
   log: '{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/training_pairs{dd_version}.txtCV_fold_divide.log'
   params: overlay="images_singularity_overlay/{dd_version}.ext3",\
           sif="images_singularity/{dd_version}.sif",\
           D="{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}",\
           load_from_h5=False,load_split_batch_pos=True,tf_order_random=False,\
           tf_num=lambda wc: get_tf_num(wc,"{exp_input_dir}/{exp_dir}/{dataset}/{dd_version}/training_pairs{dd_version}.txtTF_divide_pos.txt"),\
           jobname='deepdrimtr_{dataset}',\
           clog_prefix='training_pairs{dd_version}.txtCV_fold_divide'
   threads: 3
   resources: mem_mb=get_deepdrim_train_mem_mb, time=get_deepdrim_train_time, gpu='--gres=gpu:1'
   shell: """
     mkdir -p {params.D}; rm -rf {params.D}/training_pairs{wildcards.dd_version}.txtCV_fold_divide
     singularity exec \
        --nv \
        --no-home \
        --overlay {params.overlay}:ro \
        {params.sif} \
        /bin/bash -c \"source /ext3/env.sh; conda activate DEEPDRIM; cd {params.D}; \
                      command time -v -o time_train.txt python /ext3/DeepDRIM/DeepDRIM.py -num_batches {params.tf_num} -data_path representation_train/version11/ -output_dir training_pairs{wildcards.dd_version}.txtCV_fold_divide -cross_validation_fold_divide_file training_pairs{wildcards.dd_version}.txtCV_fold_divide.txt > training_pairs{wildcards.dd_version}.txtCV_fold_divide.log 2>&1 \"
     """

rule beeline_exp_deepdrim_train_out:# non cell-specific
   input: expand('inputs_beeline2/ALL_{ds.gt_name}/{ds.dataset}/DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 ds=EXP_PARAM_DF[EXP_PARAM_DF.gt_name=='ns'].groupby(['gt_name','dataset']).count().reset_index().itertuples())

rule beeline_exp_deepdrim_train2_out:# cell-specific
   input: expand('inputs_beeline2/ALL_{ds.gt_name}/{ds.dataset}/DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 ds=EXP_PARAM_DF[(EXP_PARAM_DF.gt_name=='cs')].groupby(['gt_name','dataset']).count().reset_index().itertuples())

rule beeline_sergio_deepdrim_train_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM/training_pairsDEEPDRIM.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim_train_out

rule beeline_exp_deepdrim5_train_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM5/training_pairsDEEPDRIM5.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim5_train_out

rule beeline_sergio_deepdrim5_train_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM5/training_pairsDEEPDRIM5.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15))
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim5_train_out

rule beeline_exp_deepdrim6_train_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM6/training_pairsDEEPDRIM6.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim6_train_out

rule beeline_sergio_deepdrim6_train_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM6/training_pairsDEEPDRIM6.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim6_train_out

rule beeline_exp_deepdrim7_train_out:
   input: expand('inputs_beeline2/{ds.exp_dir}/{ds.dataset}/DEEPDRIM7/training_pairsDEEPDRIM7.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 ds=EXP_PARAM_DF.groupby(['exp_dir','dataset']).count().reset_index().itertuples())
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_exp_deepdrim7_train_out

rule beeline_sergio_deepdrim7_train_out:
   input: expand('inputs_beeline2/{exp_dir}/net{network_i}/DEEPDRIM7/training_pairsDEEPDRIM7.txtCV_fold_divide/test_fold-0_saved_models200/keras_cnn_trained_model_DeepDRIM.h5',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i) for i in range(15)]+[str(i)+'_sh6.5_perc80' for i in range(15)])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_sergio_deepdrim7_train_out
### NEED TO REMOVE END: DEEPDRIM TRAINING is now run as part of BLRun####

##############################################
### organize CICT SERGIO results for BEELINE
##############################################
# rsync -Pavz /scratch/as15096/eric/outputs/cict_par/sens_sparsity_ranked_edges /scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs_cict_par/
rule beeline_sergio_cict_out:
    input: sergio_cict_csv='outputs_cict_par/sens_sparsity_ranked_edges/sens_sparsity.csv'
    params: sergio_cict_dir='outputs_cict_par/sens_sparsity_ranked_edges',\
            D='outputs'
    output: sergio_cict_out='outputs_cict_par/sens_sparsity_ranked_edges_format.out'
    log: 'outputs_cict_par/sens_sparsity_ranked_edges_format.log'
#    envmodules: "anaconda3/2020.07"
    run:
        sergio_cict_df = pd.read_csv(input.sergio_cict_csv)
        with open(output.sergio_cict_out,'w') as ofh:
            for index, row in sergio_cict_df.iterrows():
                exp_dir = row['dataset.dir']
                dataset = row['dataset']
                ranked_edges_file = '%d_rankedEdges.csv'%(row['job_idx'])
                ranked_edges_from = os.path.join(params.sergio_cict_dir,ranked_edges_file)
                out_dir = os.path.join(params.D,exp_dir,dataset,'CICT')
                ranked_edges_to = os.path.join(out_dir,'rankedEdges.csv')
                if (os.path.exists(ranked_edges_from)):
                    os.makedirs(out_dir,exist_ok=True)
                    shutil.copy(ranked_edges_from,ranked_edges_to)
                    ofh.write('From %s to %s\n'%(ranked_edges_from,ranked_edges_to))
                else:
                    print('%s does not exist.'%(ranked_edges_from))
                          

##########################
### BEELINE run 
##########################
# GRISLI 240000
def get_run_mem_mb(wildcards):
    switcher = {
        "L0": 200000,
        "L1": 200000,
        "L2": 200000,
        "L0_ns": 200000,
        "L1_ns": 200000,
        "L2_ns": 200000,
        "L0_lofgof": 200000,
        "L1_lofgof": 200000,
        "L2_lofgof": 200000,
        "SERGIO_DS4":16000,
        "SERGIO_DS5":16000,
        "SERGIO_DS6":16000,
        "SERGIO_DS7":16000
        }
    return switcher.get(wildcards.exp_dir, 8000)

def get_run_time(wildcards):
    switcher = {
        "L0": "6:00:00",
        "L1": "24:00:00",
        "L2": "6:00:00",
        "L0_ns": "6:00:00",
        "L1_ns": "6:00:00",
        "L2_ns": "6:00:00",
        "L0_lofgof": "6:00:00",
        "L1_lofgof": "6:00:00",
        "L2_lofgof": "6:00:00",
        "SERGIO_DS4": "0:30:00",
        "SERGIO_DS5": "0:30:00",
        "SERGIO_DS6": "0:30:00",
        "SERGIO_DS7": "0:30:00"
        }
    return switcher.get(wildcards.exp_dir, "03:00:00")

def get_run_gpu(wildcards):
    switcher = {
        'DEEPDRIM': '--gres=gpu:1',
        'DEEPDRIM4': '--gres=gpu:1',
        'DEEPDRIM7': '--gres=gpu:rtx8000:1'
        }
    return switcher.get(wildcards.algorithm, '')

# split config
#parallel --dry-run "python generateSplitConfigs.py --config config-files/config_{}.yaml --config_split_dir config-files-split/config_{}_split --replace_existing" ::: L0 L0_ns L0_lofgof L1 L1_ns L1_lofgof L2 L2_ns L2_lofgof SERGIO_DS4 SERGIO_DS5 SERGIO_DS6 SERGIO_DS7
rule beeline_run:#\\b(?!CICT).*\\b matches anything that does not start with CICT since CICT results were obtained from Abbas and thus no need to run CICT here
    input: config='config-files-split/config_{exp_dir}_split/{dataset}/{algorithm}/{run_i}/config.yaml',\
           train_repr='inputs/{exp_dir}/{dataset}/{algorithm}/{run_i}/representation_train.out',\
           predict_repr='inputs/{exp_dir}/{dataset}/{algorithm}/representation_predict.out'
    output: '{outputs_dir}/{exp_dir}/{dataset}/{algorithm,\\b(?!CICT).*\\b}/{run_i,(run_[0-9]+)?}/rankedEdges.csv'
    log: '{outputs_dir}/{exp_dir}/{dataset}/{algorithm}/{run_i}/rankedEdges.log'
    params: D="{outputs_dir}/{exp_dir}/{dataset}/{algorithm}/{run_i}",\
            jobname="blr_{exp_dir}-{dataset}-{algorithm}",\
            clog_prefix="rankedEdges"
    threads: 1
    resources: mem_mb=get_run_mem_mb, time=get_run_time, gpu=get_run_gpu
    envmodules: "anaconda3/2020.07"
    shell: """
        mkdir -p {params.D}; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
        python BLRunner.py --config {input.config} > {log} 2>&1
    """
#            python BLRunner.py --config {input} > {log} 2>&1
# find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete
#ALGORITHMS = ['CICT','DEEPDRIM','DEEPDRIM5','GENIE3','GRISLI','GRNBOOST2','GRNVBEM',
#              'LEAP','PIDC','PPCOR','RANDOM','SCNS','SCODE','SCRIBE',
#              'SINCERITIES','SINGE']

rule beeline_run_out:
   input: expand('outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.csv',\
                 exp_dir=['L0','L1','L2'],\
                 dataset=DATASET_PARAMS['ns']['dataset'],\
                 algorithm=['DEEPDRIM7'])
                 #dataset=DATASET_PARAMS['cs']['dataset'],\
                 #algorithm=[alg for alg in ALGORITHMS if alg not in ['CICT','DEEPDRIM','DEEPDRIM5','SCNS','GRNVBEM']])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_out
rule beeline_run_ns_out:
   input: expand('outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.csv',\
                 exp_dir=['L0_ns','L1_ns','L2_ns'],\
                 dataset=DATASET_PARAMS['ns']['dataset'],\
                 algorithm=['DEEPDRIM7'])
                 #algorithm=[alg for alg in ALGORITHMS if alg not in ['CICT','DEEPDRIM','DEEPDRIM5','DEEPDRIM6','SCNS','GRNVBEM']])
                 #                 dataset=['mESC'],\

#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_ns_out
rule beeline_run_lofgof_out:
   input: expand('outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.csv',\
                 exp_dir=['L0_lofgof','L1_lofgof','L2_lofgof'],\
                 dataset=DATASET_PARAMS['lofgof']['dataset'],\
                 algorithm=['DEEPDRIM7']) + \
          expand('outputs/{exp_dir}/{dataset}/{algorithm}/run_{train_i}/rankedEdges.csv',\
                 exp_dir=['L2_lofgof'],\
                 dataset=DATASET_PARAMS['lofgof']['dataset'],\
                 algorithm=['DEEPDRIM8'],\
                 train_i=range(1,11)) 
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_lofgof_out
rule beeline_run_sergio_out:
   input: expand('outputs/{exp_dir}/net{network_i}/{algorithm}/rankedEdges.csv',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15),\
                 algorithm=['DEEPDRIM7']) +\
          expand('outputs/{exp_dir}/net{network_i}/{algorithm}/rankedEdges.csv',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i)+'_sh6.5_perc80' for i in range(15)],\
                 algorithm=['DEEPDRIM7'])
#                 algorithm=[alg for alg in ALGORITHMS if alg not in ['CICT','DEEPDRIM','DEEPDRIM5']])
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_run_sergio_out

def get_eval_mem_mb(wildcards):
    switcher = {
        "epr": 16000,
        "auc": 16000,
        "auc2": 16000,
        "auc3": 16000,
        "auc4": 16000,
        "pauc4": 16000
        }
    return switcher.get(wildcards.metric, 8000)

def get_eval_time(wildcards):
    switcher = {
        "L0": "4:00:00",
        "L1": "4:00:00",
        "L2": "4:00:00",
        "L0_ns": "4:00:00",
        "L1_ns": "4:00:00",
        "L2_ns": "4:00:00",
        "L0_lofgof": "12:00:00",
        "L1_lofgof": "12:00:00",
        "L2_lofgof": "12:00:00",
        "SERGIO_DS4": "0:10:00",
        "SERGIO_DS5": "0:10:00",
        "SERGIO_DS6": "0:10:00",
        "SERGIO_DS7": "0:10:00"
        }
    return switcher.get(wildcards.exp_dir, "03:00:00")

def beeline_eval_input(wildcards):
    netout = 'outputs/{exp_dir}/{dataset}/{algorithm}/rankedEdges.csv'.format(**wildcards)
    if (os.path.exists(netout)):
        return netout
    else:
        return []
    
rule beeline_eval:
    input: netout=beeline_eval_input
    output: 'outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out'
    log: 'outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.log'
    params: config='config-files-split/config_{exp_dir}_split/{dataset}/{algorithm}/config.yaml',\
            overlay="conda_greene/overlay-5GB-200K-beeline20211104.ext3",\
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
    python BLEvaluator.py --config {params.config} --{wildcards.metric} --output_dir {params.output_dir} --output_prefix {params.output_prefix} > {log} 2>&1 \" 
    touch {output}
    """
#; find {params.D} -type f ! -name {params.clog_prefix}.slurm-out -delete

rule beeline_eval_out:
   input: expand('outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['L0','L1','L2'],\
                 dataset=DATASET_PARAMS['ns']['dataset'],\
                 algorithm=['DEEPDRIM7'],\
                 metric=['epr','auc','auc3','auc4','pauc4'])
#                 algorithm=[alg for alg in ALGORITHMS if alg not in ['DEEPDRIM','DEEPDRIM5','DEEPDRIM6','GRNVBEM']],\
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_out
rule beeline_eval_ns_out:
   input: expand('outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['L0_ns','L1_ns','L2_ns'],\
                 dataset=DATASET_PARAMS['ns']['dataset'],\
                 algorithm=['DEEPDRIM7'],\
                 metric=['epr','auc','auc3','auc4','pauc4'])
#                 algorithm=[alg for alg in ALGORITHMS if alg not in ['DEEPDRIM','DEEPDRIM5','GRNVBEM']],\
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_ns_out
rule beeline_eval_ns_out2:
   shell: "echo {rules.beeline_eval_ns_out.input}"
rule beeline_eval_lofgof_out:
   input: expand('outputs_eval/{exp_dir}/{dataset}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['L0_lofgof','L1_lofgof','L2_lofgof'],\
                 dataset=DATASET_PARAMS['lofgof']['dataset'],\
                 algorithm=['DEEPDRIM7'],\
                 metric=['epr','auc','auc3','auc4','pauc4'])
#                 algorithm=[alg for alg in ALGORITHMS if alg not in ['DEEPDRIM','DEEPDRIM5','DEEPDRIM6','SCNS','GRNVBEM']],\   
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_lofgof_out
rule beeline_eval_sergio_out:
   input: expand('outputs_eval/{exp_dir}/net{network_i}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=range(15),\
                 algorithm=['DEEPDRIM7'],\
                 metric=['epr','auc','auc3','auc4','pauc4'])+\
          expand('outputs_eval/{exp_dir}/net{network_i}/{algorithm}/{exp_dir}-{metric}.out',\
                 exp_dir=['SERGIO_DS4','SERGIO_DS5','SERGIO_DS6','SERGIO_DS7'],\
                 network_i=[str(i)+'_sh6.5_perc80' for i in range(15)],\
                 algorithm=['DEEPDRIM7'],\
                 metric=['epr','auc','auc3','auc4','pauc4'])
                 #algorithm=[alg for alg in ALGORITHMS if alg not in ['DEEPDRIM','DEEPDRIM5']],\
#snakemake -s Snakefile --profile snakefiles/profiles/slurm beeline_eval_sergio_out
rule beeline_eval_sergio_out2:
   shell: "echo {rules.beeline_eval_sergio_out.input}"
   

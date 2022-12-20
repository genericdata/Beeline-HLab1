import os
import argparse
from pathlib import Path
import yaml

import BLRun as br
yaml.warnings({'YAMLLoadWarning': False})

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Split one config file into a set of config files, each of which contains one dataset and one algorithm.')

    parser.add_argument('--config', default='config.yaml',
        help='Path to input config file to be split')

    parser.add_argument('--config_split_dir', default='config_split',
        help='Output directory of split configs')

    parser.add_argument('--replace_existing', action='store_true',
        help='Whether to replace existing config files')

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main():
    opts = parse_arguments()
    config_file = opts.config

    with open(config_file, 'r') as conf:
        #evaluation = br.ConfigParser.parse(conf)
        config_map = yaml.load(conf)
        input_settings_map = config_map['input_settings']
        output_settings_map = config_map['output_settings']
        algorithms_list = input_settings_map['algorithms']
        datasets = input_settings_map['datasets']
        
    config_map_out = config_map.copy()
    for dataset in datasets:
        for idx in range(len(algorithms_list)):
            #print(dataset)
            #print(idx)
            config_map_out['input_settings']['datasets'] = [dataset]
            config_map_out['input_settings']['algorithms'] = [algorithms_list[idx]]
            config_map_out['input_settings']['algorithms'][0]['params']['should_run'] = [True]

            # create one config for this algorithm on this dataset
            config_file_dir = os.path.join(opts.config_split_dir,dataset['name'],algorithms_list[idx]['name'])
            Path(config_file_dir).mkdir(parents=True,exist_ok=True)
            config_file_out = os.path.join(config_file_dir,'config.yaml')
            print(config_file_out)
            if os.path.exists(config_file_out):
                if opts.replace_existing:
                    with open(config_file_out,'w') as ofh:
                        yaml.dump(config_map_out,ofh,default_flow_style=False)
            else:
                with open(config_file_out,'w') as ofh:
                    yaml.dump(config_map_out,ofh,default_flow_style=False)

            # if more than 1 run is specified, create one config for each run of this algorithm on this dataset
            num_runs = int(config_map_out['input_settings']['algorithms'][0]['params'].get('numRuns',[0])[0])
            if (num_runs>0):
                for run_i in range(1,num_runs+1):
                    config_map_out['input_settings']['algorithms'][0]['params']['run_dir'] = ['run_'+str(run_i)]
                    config_file_dir = os.path.join(opts.config_split_dir,dataset['name'],algorithms_list[idx]['name'],'run_'+str(run_i))
                    Path(config_file_dir).mkdir(parents=True,exist_ok=True)
                    config_file_out = os.path.join(config_file_dir,'config.yaml')
                    print(config_file_out)
                    if os.path.exists(config_file_out):
                        if opts.replace_existing:
                            with open(config_file_out,'w') as ofh:
                                yaml.dump(config_map_out,ofh,default_flow_style=False)
                    else:
                        with open(config_file_out,'w') as ofh:
                            yaml.dump(config_map_out,ofh,default_flow_style=False)
                    

if __name__ == '__main__':
  main()

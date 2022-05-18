#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import yaml
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
import multiprocessing
from multiprocessing import Pool, cpu_count
import concurrent.futures
import yaml
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
import multiprocessing
from multiprocessing import Pool, cpu_count
import concurrent.futures
from typing import Dict, List
from BLRun.runner import Runner
import os
import pandas as pd

import BLRun as br
yaml.warnings({'YAMLLoadWarning': False})


def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml',
        help='Path to config file')

    parser.add_argument('--dataset_names', default='dream5_1',
        help='Comma delimited list of datasets from the config file to run')

    parser.add_argument('--algorithm_names', default='GENIE3',
        help='Comma delimited list of algorithms from the config file to run; \
              overrides the should_run parameter in the config file for the algorithms')
    
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
        evaluation = br.ConfigParser.parse(conf,opts)
    print('Evaluation started')

    for idx in range(len(evaluation.runners)):
        evaluation.runners[idx].generateInputs()

    for idx in range(len(evaluation.runners)):
        evaluation.runners[idx].run()

    for idx in range(len(evaluation.runners)):
        evaluation.runners[idx].parseOutput()

    print('Evaluation complete')


if __name__ == '__main__':
  main()

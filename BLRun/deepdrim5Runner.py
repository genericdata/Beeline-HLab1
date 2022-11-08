import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from sklearn.model_selection import GroupKFold
from BLRun import deepdrimRunner

def generateInputs(RunnerObj):
    deepdrimRunner.generateInputsForPredict(RunnerObj,'DEEPDRIM5',RunnerObj.inputDir)

def run(RunnerObj):
    deepdrimRunner.runForPredict(RunnerObj,'DEEPDRIM5',RunnerObj.inputDir)
                                          
def parseOutput(RunnerObj):
    deepdrimRunner.parseOutputForPredict(RunnerObj,'DEEPDRIM5',RunnerObj.inputDir)
                                          

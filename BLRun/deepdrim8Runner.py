import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from sklearn.model_selection import GroupKFold
from BLRun import deepdrimRunner

def generateInputs(RunnerObj):
    deepdrimRunner.generateInputsForPredict(RunnerObj,'DEEPDRIM8',RunnerObj.inputDir)

def run(RunnerObj):
    deepdrimRunner.runForTrainAndPredict(RunnerObj,'DEEPDRIM8',RunnerObj.inputDir)
                                          
def parseOutput(RunnerObj):
    deepdrimRunner.parseOutputForPredictFromBestModel(RunnerObj,'DEEPDRIM8',RunnerObj.inputDir,'auprc')
                                          

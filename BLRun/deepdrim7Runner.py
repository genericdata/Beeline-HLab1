import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from sklearn.model_selection import GroupKFold
from BLRun import deepdrimRunner

def generateInputs(RunnerObj):
    deepdrimRunner.generateInputsForPredict(RunnerObj,'DEEPDRIM7',RunnerObj.inputDir)

def run(RunnerObj):
    deepdrimRunner.runForTrainAndPredict(RunnerObj,'DEEPDRIM7',RunnerObj.inputDir)
                                          
def parseOutput(RunnerObj):
    deepdrimRunner.parseOutputForPredictFromBestModel(RunnerObj,'DEEPDRIM7',RunnerObj.inputDir,'auprc')
                                          

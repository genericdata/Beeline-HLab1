import os
import sys
import pandas as pd
from pathlib import Path
import numpy as np
import shutil
from sklearn.model_selection import GroupKFold
from BLRun import deepdrimRunner

def generateInputs(RunnerObj):
    algName = RunnerObj.name
    deepdrimRunner.generateInputsForPredict(RunnerObj,algName,RunnerObj.inputDir)

def run(RunnerObj):
    algName = RunnerObj.name
    deepdrimRunner.runForTrainAndPredict(RunnerObj,algName,RunnerObj.inputDir)
                                          
def parseOutput(RunnerObj):
    algName = RunnerObj.name
    deepdrimRunner.parseOutputForPredictFromBestModel(RunnerObj,algName,RunnerObj.inputDir,'auprc')
                                          

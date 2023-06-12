import BLRun.scodeRunner as SCODE
import BLRun.scnsRunner as SCNS
import BLRun.sinceritiesRunner as SINCERITIES
import BLRun.pidcRunner as PIDC
import BLRun.grnvbemRunner as GRNVBEM
import BLRun.genie3Runner as GENIE3
import BLRun.grnboost2Runner as GRNBOOST2
import BLRun.leapRunner as LEAP
import BLRun.jump3Runner as JUMP3
import BLRun.ppcorRunner as PPCOR
import BLRun.grisliRunner as GRISLI
import BLRun.singeRunner as SINGE
import BLRun.scribeRunner as SCRIBE
import BLRun.cictRunner as CICT
import BLRun.randomRunner as RANDOM
import BLRun.deepdrimRunner as DEEPDRIM
import BLRun.deepdrim4Runner as DEEPDRIM4
import BLRun.deepdrim5Runner as DEEPDRIM5
import BLRun.deepdrim6Runner as DEEPDRIM6
import BLRun.deepdrim7Runner as DEEPDRIM7
import BLRun.deepdrim8Runner as DEEPDRIM8
import BLRun.cnncRunner as CNNC
import BLRun.inferelator3Runner as INFERELATOR3

from pathlib import Path


InputMapper = {'SCODE':SCODE.generateInputs,
               'SINCERITIES':SINCERITIES.generateInputs,
               'SCNS':SCNS.generateInputs,
               'PIDC':PIDC.generateInputs,
               'GRNVBEM':GRNVBEM.generateInputs,
               'GENIE3':GENIE3.generateInputs,
               'GRNBOOST2':GRNBOOST2.generateInputs,
               'LEAP':LEAP.generateInputs,
               'JUMP3':JUMP3.generateInputs,
               'PPCOR':PPCOR.generateInputs,
               'GRISLI':GRISLI.generateInputs,
               'SINGE':SINGE.generateInputs,
               'SCRIBE':SCRIBE.generateInputs,
               'CICT':CICT.generateInputs,
               'RANDOM':RANDOM.generateInputs,
               'DEEPDRIM':DEEPDRIM.generateInputs,
               'DEEPDRIM4':DEEPDRIM4.generateInputs,
               'DEEPDRIM5':DEEPDRIM5.generateInputs,
               'DEEPDRIM6':DEEPDRIM6.generateInputs,
               'DEEPDRIM7':DEEPDRIM7.generateInputs,
               'DEEPDRIM8':DEEPDRIM8.generateInputs,
               'DEEPDRIM72':DEEPDRIM7.generateInputs,
               'DEEPDRIM7_v2':DEEPDRIM7.generateInputs,
               'DEEPDRIM72_v2':DEEPDRIM7.generateInputs,
               'CNNC':CNNC.generateInputs,
               'INFERELATOR31':INFERELATOR3.generateInputs,
               'INFERELATOR32':INFERELATOR3.generateInputs,
               'INFERELATOR33':INFERELATOR3.generateInputs,
               'INFERELATOR34':INFERELATOR3.generateInputs,
               'INFERELATOR35':INFERELATOR3.generateInputs,
               'INFERELATOR36':INFERELATOR3.generateInputs,
               'INFERELATOR37':INFERELATOR3.generateInputs,
               'INFERELATOR38':INFERELATOR3.generateInputs,
               'INFERELATOR31_v2':INFERELATOR3.generateInputs,
               'INFERELATOR32_v2':INFERELATOR3.generateInputs,
               'INFERELATOR33_v2':INFERELATOR3.generateInputs,
               'INFERELATOR34_v2':INFERELATOR3.generateInputs,
               'INFERELATOR35_v2':INFERELATOR3.generateInputs,
               'INFERELATOR36_v2':INFERELATOR3.generateInputs,
               'INFERELATOR37_v2':INFERELATOR3.generateInputs,
               'INFERELATOR38_v2':INFERELATOR3.generateInputs,
               'INFERELATOR39_v2':INFERELATOR3.generateInputs,
               'INFERELATOR310_v2':INFERELATOR3.generateInputs}






AlgorithmMapper = {'SCODE':SCODE.run,
            'SINCERITIES':SINCERITIES.run,
            'SCNS':SCNS.run,
            'PIDC':PIDC.run,
            'GRNVBEM':GRNVBEM.run,
            'GENIE3':GENIE3.run,
            'GRNBOOST2':GRNBOOST2.run,
            'LEAP':LEAP.run,
            'JUMP3':JUMP3.run,
            'PPCOR':PPCOR.run,
            'GRISLI':GRISLI.run,
            'SINGE':SINGE.run,
            'SCRIBE':SCRIBE.run,
            'CICT':CICT.run,
            'RANDOM':RANDOM.run,
            'DEEPDRIM':DEEPDRIM.run,
            'DEEPDRIM4':DEEPDRIM4.run,
            'DEEPDRIM5':DEEPDRIM5.run,
            'DEEPDRIM6':DEEPDRIM6.run,
            'DEEPDRIM7':DEEPDRIM7.run,
            'DEEPDRIM72':DEEPDRIM7.run,
            'DEEPDRIM8':DEEPDRIM8.run,
            'DEEPDRIM7_v2':DEEPDRIM7.run,
            'DEEPDRIM72_v2':DEEPDRIM7.run,
            'CNNC':CNNC.run,
            'INFERELATOR31':INFERELATOR3.run,
            'INFERELATOR32':INFERELATOR3.run,
            'INFERELATOR33':INFERELATOR3.run,
            'INFERELATOR34':INFERELATOR3.run,
            'INFERELATOR35':INFERELATOR3.run,
            'INFERELATOR36':INFERELATOR3.run,
            'INFERELATOR37':INFERELATOR3.run,
            'INFERELATOR38':INFERELATOR3.run,                   
            'INFERELATOR31_v2':INFERELATOR3.run,
            'INFERELATOR32_v2':INFERELATOR3.run,
            'INFERELATOR33_v2':INFERELATOR3.run,
            'INFERELATOR34_v2':INFERELATOR3.run,
            'INFERELATOR35_v2':INFERELATOR3.run,
            'INFERELATOR36_v2':INFERELATOR3.run,
            'INFERELATOR37_v2':INFERELATOR3.run,
            'INFERELATOR38_v2':INFERELATOR3.run,
            'INFERELATOR39_v2':INFERELATOR3.run,
            'INFERELATOR310_v2':INFERELATOR3.run}





OutputParser = {'SCODE':SCODE.parseOutput, 
            'SINCERITIES':SINCERITIES.parseOutput,
            'SCNS':SCNS.parseOutput,
            'PIDC':PIDC.parseOutput,
            'GRNVBEM':GRNVBEM.parseOutput,
            'GENIE3':GENIE3.parseOutput,
            'GRNBOOST2':GRNBOOST2.parseOutput,
            'LEAP': LEAP.parseOutput,
            'JUMP3': JUMP3.parseOutput,
            'PPCOR':PPCOR.parseOutput,
            'GRISLI':GRISLI.parseOutput,
            'SINGE':SINGE.parseOutput,
            'SCRIBE':SCRIBE.parseOutput,
            'CICT':CICT.parseOutput,
            'RANDOM':RANDOM.parseOutput,
            'DEEPDRIM':DEEPDRIM.parseOutput,
            'DEEPDRIM4':DEEPDRIM4.parseOutput,
            'DEEPDRIM5':DEEPDRIM5.parseOutput,
            'DEEPDRIM6':DEEPDRIM6.parseOutput,
            'DEEPDRIM7':DEEPDRIM7.parseOutput,
            'DEEPDRIM72':DEEPDRIM7.parseOutput,
            'DEEPDRIM8':DEEPDRIM8.parseOutput,
            'DEEPDRIM7_v2':DEEPDRIM7.parseOutput,
            'DEEPDRIM72_v2':DEEPDRIM7.parseOutput,
            'CNNC':CNNC.parseOutput,
            'INFERELATOR31':INFERELATOR3.parseOutput,
            'INFERELATOR32':INFERELATOR3.parseOutput,
            'INFERELATOR33':INFERELATOR3.parseOutput,
            'INFERELATOR34':INFERELATOR3.parseOutput,
            'INFERELATOR35':INFERELATOR3.parseOutput,
            'INFERELATOR36':INFERELATOR3.parseOutput,
            'INFERELATOR37':INFERELATOR3.parseOutput,
            'INFERELATOR38':INFERELATOR3.parseOutput,
            'INFERELATOR31_v2':INFERELATOR3.parseOutput,
            'INFERELATOR32_v2':INFERELATOR3.parseOutput,
            'INFERELATOR33_v2':INFERELATOR3.parseOutput,
            'INFERELATOR34_v2':INFERELATOR3.parseOutput,
            'INFERELATOR35_v2':INFERELATOR3.parseOutput,
            'INFERELATOR36_v2':INFERELATOR3.parseOutput,
            'INFERELATOR37_v2':INFERELATOR3.parseOutput,
            'INFERELATOR38_v2':INFERELATOR3.parseOutput,
            'INFERELATOR39_v2':INFERELATOR3.parseOutput,
            'INFERELATOR310_v2':INFERELATOR3.parseOutput}






class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    def __init__(self,
                params):
        self.name = params['name']
        self.inputDir = params['inputDir']
        self.fullInputDir = params['fullInputDir']
        self.params = params['params']
        self.exprData = params['exprData']
        self.cellData = params['cellData']
        self.trueEdges = params['trueEdges']
        self.singularityImage = params['singularityImage']
        self.singularityOverlay = params.get('singularityOverlay','')
        self.singularityGPUFlag = params.get('singularityGPUFlag','')
        
    def generateInputs(self):
        InputMapper[self.name](self)
        
        
    def run(self):
        AlgorithmMapper[self.name](self)


    def parseOutput(self):
        OutputParser[self.name](self)

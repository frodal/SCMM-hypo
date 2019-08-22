#import numpy as np
#import sys
import os
import glob
import enum
from pathlib import Path
##----------------------------------------------------------------------
## 
##----------------------------------------------------------------------
def read_input(filename):
    fp = open(filename,'r')
    lines = fp.readlines()
    fp.close()
    return lines

def write_job(lines,string):
    fp = open('jobaba','w')
    for line in lines:
        fp.write(line.replace('<<k>>',string))
    fp.close()
    return
##----------------------------------------------------------------------
## Abaqus solver Enum
##----------------------------------------------------------------------
@enum.unique
class AbaqusSolver(enum.Enum):
    explicit = enum.auto()
    implicit = enum.auto()
##----------------------------------------------------------------------
## Material class
##----------------------------------------------------------------------
class Material:

    nProps = 17
    nStatev = 28

    def __init__(self,name,props,density):
        self.name = name
        self.props = props
        self.density = density
        if len(props)!=self.nProps:
            print('Warning!')
            print('Wrong number of material properties given!')
            print('Please supply '+str(self.nProps)+' material properties')
##----------------------------------------------------------------------
## Test class
##----------------------------------------------------------------------
class Test:
    def __init__(self,name='Test'):
        self.name = name

class AbaqusTest(Test):

    def __init__(self,name,inputfile,solver,material):
        super().__init__(name)
        
        self.inputfile = Path(__file__).parent.joinpath(inputfile)
        if not os.path.isfile(self.inputfile):
            print('The inputfile provided could not be found!')
            print(self.inputfile)

        self.solver = solver
        if self.solver != AbaqusSolver.explicit and self.solver != AbaqusSolver.implicit :
            print('Error!')
            print('Unknown solver!')

        self.material = material
        if not isinstance(material,Material):
            print('Error!')
            print('Unknown material!')

        jobabaFile = Path(__file__).parent.joinpath('Abaqus/Snurre-jobaba')
        self.jobaba = jobabaFile.read_text()

class FortranTest(Test):
    def __init__(self,name):
        super().__init__(name)
##----------------------------------------------------------------------
## 
##----------------------------------------------------------------------
#input_file = 'SLM_input.k'
#specimen   = 'SMOOTH'
#
def main():
    material = Material('Cube',[106430.,  60350.,  28210.,    0.01,   0.005, 46.7301,     1.4,      1.,
                                     0.,      0.,      0.,      2., 411.256, 104.029, 1.35459,      0.,
                                     1.],2.7e-9)
    test = AbaqusTest('SimpleShear','Abaqus/SimpleShear/SimpleShear-Explicit.inp',AbaqusSolver.explicit,material)
    print(test.name)
    # joblines = read_input('jobSLM')
    # #
    # for file in glob.glob('F_*'+specimen+'*.csv'):
    #     filename = file.split('.csv')[0]
    #     element  = filename.split('_')[2]
    #     jobname  = 'jobSLM_'+specimen+'_'+element
    #     gradfile = 'F_'+specimen+'-1_'+element+'.csv'
    #     if 'SLM_'+specimen+'_'+element+'.out' not in glob.glob('*.out'):
    #         write_job(joblines,specimen+'_'+element)
    #         os.system('qsub'+' '+jobname+' '+input_file+' '+gradfile) # qsub jobSLM_DIABOLO_1900 SLM_6082_O.k F_DIABOLO-1_6082.25_O_1900.csv

if __name__=='__main__':
    main()
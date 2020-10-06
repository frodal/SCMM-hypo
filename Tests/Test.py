##----------------------------------------------------------------------
## Import libraries
##----------------------------------------------------------------------
import os
import enum
from pathlib import Path
import shutil
import argparse
import glob
import pandas as pd
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
##----------------------------------------------------------------------
## Context manager class
##----------------------------------------------------------------------
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
##----------------------------------------------------------------------
## Print colored text class
##----------------------------------------------------------------------
class printColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
##----------------------------------------------------------------------
## Abaqus solver Enum
##----------------------------------------------------------------------
@enum.unique
class AbaqusSolver(enum.Enum):
    # Generates the value to the auto method
    def _generate_next_value_(name, start, count, last_values):
        return name
    # Enums
    Explicit = enum.auto()
    Implicit = enum.auto()
##----------------------------------------------------------------------
## Material class
##----------------------------------------------------------------------
class EulerAngles:
    # Constructor
    def __init__(self,phi1,PHI,phi2):
        self.Set(phi1,PHI,phi2)
    # Sets the Euler angles
    def Set(self,phi1,PHI,phi2):
        self.phi1 = phi1
        self.PHI  = PHI
        self.phi2 = phi2
##----------------------------------------------------------------------
## Material class
##----------------------------------------------------------------------
class Material:
    # static class members
    nProps = 21
    nStatev = 30
    nDelete = 30

    # Constructor
    def __init__(self,name,abaqusMaterialName,density,props):
        self.name = name
        self.abaqusMaterialName = abaqusMaterialName
        self.density = density
        self.props = props
        assert(len(props)==Material.nProps),(
            'Wrong number of material properties given in material: '+self.name+'\n'+
            'Please supply '+str(Material.nProps)+' material properties')
        self.eulerAngles = EulerAngles(self.props[8],self.props[9],self.props[10])
    
    # Sets the Euler angles
    def SetEulerAngles(self,phi1,PHI,phi2):
        self.eulerAngles.Set(phi1,PHI,phi2)
        self.props[8] = phi1
        self.props[9] = PHI
        self.props[10]= phi2
    
    # Writes an Abaqus material card for the current material
    def WriteMaterailCard(self,writePath):
        materialCardPath = Path(writePath).joinpath('material.inp')
        with open(materialCardPath,'w') as fp:
            fp.write('*Material, name={}\n'.format(self.abaqusMaterialName))
            fp.write('*Density\n')
            fp.write('{:8},\n'.format(self.density))
            fp.write('*Depvar, delete={}\n'.format(self.nDelete))
            fp.write('{:8},\n'.format(self.nStatev))
            fp.write('*User Material, constants={}\n'.format(self.nProps))
            for i in range(self.nProps):
                if (i+1)%8 == 0:
                    fp.write('{:8}'.format(self.props[i]))
                    if (i+1)!=self.nProps:
                        fp.write('\n')
                else:
                    fp.write('{:8},'.format(self.props[i]))
##----------------------------------------------------------------------
## Test class
##----------------------------------------------------------------------
class Test:
    # Constructor
    def __init__(self,name='Test'):
        self.name = name
##----------------------------------------------------------------------
class AbaqusTest(Test):
    # Constructor
    def __init__(self,name,inputfile,solver,postScript,material,ncpu):
        # Calls the base (super) class's constructor
        super().__init__(name)
        # Current directory
        pythonPath = Path(__file__).parent
        # Saves the path to the input file provided
        self.inputfile = pythonPath.joinpath(inputfile)
        assert (self.inputfile.exists()),(
            'The inputfile provided could not be found: '+str(self.inputfile))
        # Saves the solver to use in Abaqus
        self.solver = solver
        assert (self.solver == AbaqusSolver.Explicit or self.solver == AbaqusSolver.Implicit),(
            'Unknown solver: '+str(self.solver))
        # Saves the path to the post-processing script provided as an absolute path
        self.postScript = pythonPath.joinpath(postScript).absolute()
        assert (self.postScript.exists()),(
            'The post-processing script provided could not be found: '+str(self.postScript))
        # Saves the material to be used
        self.material = material
        assert (isinstance(self.material,Material)),(
            'Unknown material: '+str(self.material))
        # Saves the number of cores to be used
        self.ncpu = ncpu
        assert (isinstance(self.ncpu,int)and(self.ncpu>0)),(
            'ncpu should be a positive integer: '+str(self.ncpu))
        # Sets up the Abaqus folder path, the test working directory path, and reference data path
        self.abaqusPath = pythonPath.joinpath('Abaqus')
        self.testPath = self.abaqusPath.joinpath('WorkingDirectory').joinpath(self.name).joinpath(
            self.solver.name).joinpath(self.material.name)
        self.referencePath = pythonPath.joinpath('ReferenceData').joinpath(
            self.name).joinpath(self.solver.name).joinpath(self.material.name)
        # Current input file path
        self.currentInputFile = self.testPath.joinpath(self.name+'.inp')
        # Initialize
        self.passed = False
        self.residual = np.inf

    # Sets up the Abaqus job files
    def _SetupAbaqusJob(self):
        # Creates the test working directory
        if not self.testPath.exists():
            self.testPath.mkdir(parents=True)
        # Copies the Abaqus input file to the test working directory
        shutil.copy(self.inputfile,self.currentInputFile)
        # Writes a material card to the test working directory
        self.material.WriteMaterailCard(self.testPath)
        # Copies the Abaqus environment file to the test working directory
        shutil.copy(self.abaqusPath.joinpath('abaqus_v6.env'),
                    self.testPath.joinpath('abaqus_v6.env'))

    # Runs the test on the CASA cluster. Note that the script should be run from the cluster
    def _RunOnCasa(self):
        # Sets up the test
        self._SetupAbaqusJob()
        # Reads the template jobaba file and writes a jobaba file to the test working directory
        jobabaName = 'jobaba'
        jobabaTemplateFile = self.abaqusPath.joinpath('CASA-jobaba')
        if self.ncpu>24:
            self.ncpu = 24
        jobabaContent = jobabaTemplateFile.read_text().replace('<<jobName>>',self.name).replace('<<ncpu>>',str(self.ncpu))
        jobabaFile = self.testPath.joinpath(jobabaName)
        jobabaFile.write_text(jobabaContent)
        # Submit the job to the queue on the cluster
        with cd(self.testPath):
            os.system('sbatch '+jobabaName)

    # Runs the test on the current PC
    def _RunHere(self,interactiveOff=False):
        # Sets up the test
        self._SetupAbaqusJob()
        # Run the abaqus solver
        with cd(self.testPath):
            if interactiveOff:
                os.system('abaqus double job='+str(self.name)+' cpus='+str(self.ncpu))
            else:
                os.system('abaqus double job='+str(self.name)+' cpus='+str(self.ncpu)+' interactive')
    
    # Runs the test
    def Run(self,onCluster=False,interactiveOff=False):
        if onCluster:
            self._RunOnCasa()
        else:
            self._RunHere(interactiveOff)
    
    # Post-process the test
    def Process(self,shouldPlot=False):
        self.passed = False
        self.residual = np.inf
        with cd(self.testPath):
            # Call abaqus python with the post-processing script
            os.system('abaqus python '+str(self.postScript))
        # Read test result
        try:
            testData = pd.read_csv(self.testPath.joinpath('Result.csv')).T.to_numpy()
        except:
            return self.passed
        # Read reference data
        try:
            referenceData = pd.read_csv(self.referencePath.joinpath('Result.csv')).T.to_numpy()
        except:
            assert (False),'Could not read the test reference data'
            return self.passed
        # Assume that the test result contains x and y data
        xRef  = referenceData[0]
        yRef  = referenceData[1]
        nRef  = len(xRef)
        # Check length of data
        if len(testData[0]) != nRef:
            return self.passed
        # Creates interpolation functions as to evaluate the difference at the same x-values 
        fTest = interpolate.interp1d(testData[0],testData[1],bounds_error=False,fill_value=(testData[1][0],testData[1][-1]))
        fRef  = interpolate.interp1d(xRef,yRef)
        # 
        x     = np.linspace(xRef.min(),xRef.max(),nRef)
        yTest = fTest(x)
        y     = fRef(x)
        # Calculates the residual of the test
        self.residual = np.sum(np.power((y-yTest)/(y.max()-y.min()),2))/nRef
        self.passed = self.residual<0.0001
        # Plot results
        if shouldPlot:
            plt.figure()
            plt.plot(xRef, yRef)
            plt.plot(testData[0],testData[1])
        return self.passed
##----------------------------------------------------------------------
class FortranTest(Test):
    # Constructor
    def __init__(self,name):
        # Calls the base (super) class's constructor
        super().__init__(name)
    
    # Runs the test
    def Run(self):
        pass

    # Post-process the test
    def Process(self):
        pass
    # TODO: Implement tests to test parts of the Fortran code
##----------------------------------------------------------------------
## Clean working directories
##----------------------------------------------------------------------
def Clean():
    workingDir = Path(__file__).parent.joinpath('Abaqus').joinpath('WorkingDirectory')
    if workingDir.exists():
        try:
            shutil.rmtree(workingDir)
            print('Working directory removed: '+str(workingDir))
        except:
            print('Could not remove the working directory: '+str(workingDir)+'\n'+'Try again later!')
##----------------------------------------------------------------------
## Run tests
##----------------------------------------------------------------------
def RunTests(tests,onCluster=False,interactiveOff=False):
    for test in tests:
        test.Run(onCluster,interactiveOff)
##----------------------------------------------------------------------
## Post-process tests
##----------------------------------------------------------------------
def PostProcess(tests,shouldPlot=False):
    for test in tests:
        test.Process(shouldPlot)
        name = test.name+'-'+test.solver.name+'-'+test.material.name
        if test.passed:
            print('{}PASSED{} test {:40} residual = {:e}'.format(printColors.OKGREEN,printColors.ENDC,name,test.residual))
        else:
            print('{}FAILED{} test {:40} residual = {:e}'.format(printColors.FAIL,printColors.ENDC,name,test.residual))
        if shouldPlot:
            plt.show()
##----------------------------------------------------------------------
## Creates the SimpleShear tests
##----------------------------------------------------------------------
def CreateSimpleShearTests():
    # Material density
    density = 2.7e-9
    # Euler angles to be used
    eulerAngles = [ EulerAngles(     0.0,     0.0,     0.0),
                    EulerAngles(    45.0,     0.0,     0.0),
                    EulerAngles(129.2315,114.0948,333.4349),
                    EulerAngles(206.5651,114.0948,50.76848)]
    
    # Creates Kalidindi materials
    kalidindiMaterialNames = ['000','4500','inverted','other']
    kalidindiMaterials = []
    for materialName,eAngles in zip(kalidindiMaterialNames,eulerAngles):
        kalidindiMaterials.append(
            Material(materialName,'Material-1',density,
            [    106430.,      60350.,       28210., 0.01,   0.005, 46.7301,     1.4, 1.,
            eAngles.phi1, eAngles.PHI, eAngles.phi2,   2., 411.256, 104.029, 1.35459, 0.,
                      1.,         0.0,          1.0,  0.0,     0.0]))
    
    # Creates Voce hardening materials
    voceMaterialNames = ['000-Voce','4500-Voce','inverted-Voce','other-Voce']
    voceMaterials = []
    for materialName,eAngles in zip(voceMaterialNames,eulerAngles):
        voceMaterials.append(
            Material(materialName,'Material-1',density,
            [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    1.,
            eAngles.phi1, eAngles.PHI, eAngles.phi2,   1., 20.48,   18.07, 157.3, 39.11,
                      1.,         0.0,          1.0,  0.0,   0.0]))
    
    # Creates Voce hardening materials with RT-damage
    voceRTMaterialNames = ['000-Voce-RT','4500-Voce-RT','inverted-Voce-RT','other-Voce-RT']
    voceRTMaterials = []
    for materialName,eAngles in zip(voceRTMaterialNames,eulerAngles):
        voceRTMaterials.append(
            Material(materialName,'Material-1',density,
            [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    1.,
            eAngles.phi1, eAngles.PHI, eAngles.phi2,   1., 20.48,   18.07, 157.3, 39.11,
                      1.,        0.01,          0.6,  1.5,   1.0]))
    
    # Add different tests to be run
    tests = []
    # Add SimpleShear tests using Abaqus/Explicit and the Kalidindi materials
    for material in kalidindiMaterials:
        tests.append(AbaqusTest('SimpleShear','Abaqus/SimpleShearTest/SimpleShear-Explicit.inp',
                    AbaqusSolver.Explicit,'Abaqus/SimpleShearTest/SimpleShearExtract.py',
                    material,1))
    # Add SimpleShear tests using Abaqus/Explicit and the Voce hardening materials
    for material in voceMaterials:
        tests.append(AbaqusTest('SimpleShear','Abaqus/SimpleShearTest/SimpleShear-Explicit.inp',
                    AbaqusSolver.Explicit,'Abaqus/SimpleShearTest/SimpleShearExtract.py',
                    material,1))
    # Add SimpleShear tests using Abaqus/Standard and the Kalidindi materials
    for material in kalidindiMaterials:
        tests.append(AbaqusTest('SimpleShear','Abaqus/SimpleShearTest/SimpleShear-Implicit.inp',
                    AbaqusSolver.Implicit,'Abaqus/SimpleShearTest/SimpleShearExtract.py',
                    material,1))
    # Add SimpleShear tests using Abaqus/Explicit and the Voce hardening materials with RT-damage
    for material in voceRTMaterials:
        tests.append(AbaqusTest('SimpleShear','Abaqus/SimpleShearTest/SimpleShear-Explicit.inp',
                    AbaqusSolver.Explicit,'Abaqus/SimpleShearTest/SimpleShearExtract.py',
                    material,1))
    
    return tests
##----------------------------------------------------------------------
## Creates the uniaxial tension tests
##----------------------------------------------------------------------
def CreateUniaxialTensionTests():
    # Material density
    density = 2.7e-9
    # Euler angles to be used
    eulerAngles = [EulerAngles(  0.0,  0.0,  0.0),
                   EulerAngles( 45.0,  0.0,  0.0),
                   EulerAngles(  0.0, 45.0,  0.0),
                   EulerAngles( 35.0, 45.0,  0.0)]
    
    # Creates Voce hardening materials
    MaterialNames = ['000-Voce','4500-Voce','0450-Voce','35450-Voce']
    Materials = []
    for materialName,eAngles in zip(MaterialNames,eulerAngles):
        Materials.append(
            Material(materialName,'AL',density,
            [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    1.,
            eAngles.phi1, eAngles.PHI, eAngles.phi2,   1., 20.48,   18.07, 157.3, 39.11,
                      2.,         0.0,          1.0,  0.0,   0.0]))
    
    # Creates Voce hardening materials
    rtMaterialNames = ['000-Voce-RT','4500-Voce-RT','0450-Voce-RT','35450-Voce-RT']
    rtMaterials = []
    for materialName,eAngles in zip(rtMaterialNames,eulerAngles):
        rtMaterials.append(
            Material(materialName,'AL',density,
            [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    1.,
            eAngles.phi1, eAngles.PHI, eAngles.phi2,   1., 20.48,   18.07, 157.3, 39.11,
                      2.,        0.01,          0.6,  1.5,   1.0]))
    
    # Add different tests to be run
    tests = []
    # Add Uniaxial tension tests using Abaqus/Explicit and the Voce hardening materials
    for material in Materials:
        tests.append(AbaqusTest('UniaxialTension','Abaqus/UniaxialTensionTest/UniaxialTension-Explicit.inp',
                    AbaqusSolver.Explicit,'Abaqus/UniaxialTensionTest/UniaxialTensionExtract.py',
                    material,1))
    # Add Uniaxial tension tests using Abaqus/Standard and the Voce hardening materials
    for material in Materials:
        tests.append(AbaqusTest('UniaxialTension','Abaqus/UniaxialTensionTest/UniaxialTension-Implicit.inp',
                    AbaqusSolver.Implicit,'Abaqus/UniaxialTensionTest/UniaxialTensionExtract.py',
                    material,1))
    # Add Uniaxial tension tests using Abaqus/Explicit and the Voce hardening materials with RT-damage
    for material in rtMaterials:
        tests.append(AbaqusTest('UniaxialTension','Abaqus/UniaxialTensionTest/UniaxialTension-Explicit.inp',
                    AbaqusSolver.Explicit,'Abaqus/UniaxialTensionTest/UniaxialTensionExtract.py',
                    material,1))
    
    return tests
##----------------------------------------------------------------------
## Creates the polycrystal tests
##----------------------------------------------------------------------
def CreatePolycrystalTests():
    # Material density
    density = 2.7e-9
    
    # Creates Voce hardening materials
    materialName = 'Voce'
    material = Material(materialName,'AL',density,
        [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    2.,
                 0.0,         0.0,          0.0,   1., 20.48,   18.07, 157.3, 39.11,
                  2.,         0.0,          1.0,  0.0,   0.0])
    
    materialRTName = 'Voce-RT'
    materialRT = Material(materialRTName,'AL',density,
        [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    2.,
                 0.0,         0.0,          0.0,   1., 20.48,   18.07, 157.3, 39.11,
                  2.,        0.01,          0.6,  1.5,   1.0])
    
    # Add different tests to be run
    tests = []
    # Add polycrystal tests using Abaqus/Explicit and the Voce hardening materials
    tests.append(AbaqusTest('Polycrystal','Abaqus/PolycrystalTest/PolycrystalUniaxialTension-Explicit.inp',
                AbaqusSolver.Explicit,'Abaqus/PolycrystalTest/PolycrystalExtract.py',
                material,1))
    # Add polycrystal tests using Abaqus/Standard and the Voce hardening materials
    tests.append(AbaqusTest('Polycrystal','Abaqus/PolycrystalTest/PolycrystalUniaxialTension-Implicit.inp',
                AbaqusSolver.Implicit,'Abaqus/PolycrystalTest/PolycrystalExtract.py',
                material,8))
    # Add polycrystal tests using Abaqus/Explicit and the Voce hardening materials with RT-damage
    tests.append(AbaqusTest('Polycrystal','Abaqus/PolycrystalTest/PolycrystalUniaxialTension-Explicit.inp',
                AbaqusSolver.Explicit,'Abaqus/PolycrystalTest/PolycrystalExtract.py',
                materialRT,1))
    
    return tests
##----------------------------------------------------------------------
## Creates the plane strain tests
##----------------------------------------------------------------------
def CreatePlaneStrainTests():
    # Material density
    density = 2.7e-9
    # Euler angles to be used
    eulerAngles = [EulerAngles(  0.0,  0.0,  0.0),
                   EulerAngles( 45.0,  0.0,  0.0),
                   EulerAngles(  0.0, 45.0,  0.0),
                   EulerAngles( 35.0, 45.0,  0.0)]
    
    # Creates Voce hardening materials
    MaterialNames = ['000-Voce','4500-Voce','0450-Voce','35450-Voce']
    Materials = []
    for materialName,eAngles in zip(MaterialNames,eulerAngles):
        Materials.append(
            Material(materialName,'CP',density,
            [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    1.,
            eAngles.phi1, eAngles.PHI, eAngles.phi2,   1., 20.48,   18.07, 157.3, 39.11,
                      2.]))
    
    # Add different tests to be run
    tests = []
    # Add plane strain tests using Abaqus/Explicit and the Voce hardening materials
    for material in Materials:
        tests.append(AbaqusTest('PlaneStrain','Abaqus/PlaneStrainTest/PlaneStrain-Explicit.inp',
                    AbaqusSolver.Explicit,'Abaqus/PlaneStrainTest/PlaneStrainExtract.py',
                    material,1))
    # Add plane strain tests using Abaqus/Standard and the Voce hardening materials
    for material in Materials:
        tests.append(AbaqusTest('PlaneStrain','Abaqus/PlaneStrainTest/PlaneStrain-Implicit.inp',
                    AbaqusSolver.Implicit,'Abaqus/PlaneStrainTest/PlaneStrainExtract.py',
                    material,1))
    
    return tests
##----------------------------------------------------------------------
## Creates the axisymmetric tests
##----------------------------------------------------------------------
def CreateAxisymmetricTests():
    # Material density
    density = 2.7e-9
    # Euler angles to be used
    eulerAngles = [EulerAngles(  0.0,  0.0,  0.0),
                   EulerAngles( 45.0,  0.0,  0.0),
                   EulerAngles(  0.0, 45.0,  0.0),
                   EulerAngles( 35.0, 45.0,  0.0)]
    
    # Creates Voce hardening materials
    MaterialNames = ['000-Voce','4500-Voce','0450-Voce','35450-Voce']
    Materials = []
    for materialName,eAngles in zip(MaterialNames,eulerAngles):
        Materials.append(
            Material(materialName,'CP',density,
            [    106430.,      60350.,       28210., 0.01, 0.005, 46.7301,   1.4,    1.,
            eAngles.phi1, eAngles.PHI, eAngles.phi2,   1., 20.48,   18.07, 157.3, 39.11,
                      2.]))
    
    # Add different tests to be run
    tests = []
    # Add axisymmetric tests using Abaqus/Explicit and the Voce hardening materials
    for material in Materials:
        tests.append(AbaqusTest('Axisymmetric','Abaqus/AxisymmetricTest/Axisymmetric-Explicit.inp',
                    AbaqusSolver.Explicit,'Abaqus/AxisymmetricTest/AxisymmetricExtract.py',
                    material,1))
    # Add axisymmetric tests using Abaqus/Standard and the Voce hardening materials
    for material in Materials:
        tests.append(AbaqusTest('Axisymmetric','Abaqus/AxisymmetricTest/Axisymmetric-Implicit.inp',
                    AbaqusSolver.Implicit,'Abaqus/AxisymmetricTest/AxisymmetricExtract.py',
                    material,1))
    
    return tests
##----------------------------------------------------------------------
## Main
##----------------------------------------------------------------------
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run tests for the SCMM-hypo subroutine.')
    parser.add_argument('action',type=str,choices=['run','clean','post'],
                        help='Choose what to do.'+
                        ' "run": For running the tests.'+
                        ' "clean": For cleaning the test working directories.'+
                        ' "post": For post-processing the results.')
    parser.add_argument('--casa',default=False,const=True,action='store_const',
                        help='Add this flag when running the tests on the CASA cluster.')
    parser.add_argument('--interactive_off',default=False,const=True,action='store_const',
                        help='Add this flag to turn off interactive mode of the Abaqus analyses.')
    parser.add_argument('--plot',default=False,const=True,action='store_const',
                        help='Add this flag to plot the reference data and the test data during post-processing.')
    args = parser.parse_args()
    onCluster = args.casa
    action = args.action
    shouldPlot = args.plot
    interactiveOff = args.interactive_off
    
    # Creates the tests
    tests = []
    tests += CreateSimpleShearTests()
    tests += CreateUniaxialTensionTests()
    tests += CreatePlaneStrainTests()
    tests += CreateAxisymmetricTests()
    tests += CreatePolycrystalTests()

    # Do stuff
    if action=='run':
        RunTests(tests,onCluster,interactiveOff)
    elif action=='clean':
        Clean()
    elif action=='post':
        PostProcess(tests,shouldPlot)
##----------------------------------------------------------------------
## Entry point
##----------------------------------------------------------------------
if __name__=='__main__':
    main()
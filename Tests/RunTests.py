##----------------------------------------------------------------------
## Import libraries
##----------------------------------------------------------------------
import os
import enum
from pathlib import Path
import shutil
##----------------------------------------------------------------------
## Abaqus solver Enum
##----------------------------------------------------------------------
@enum.unique
class AbaqusSolver(enum.Enum):
    # Generates the value to the auto method
    def _generate_next_value_(name, start, count, last_values):
        return name
    #Enums
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
    nProps = 17
    nStatev = 28
    # Constructor
    def __init__(self,name,density,props):
        self.name = name
        self.density = density
        self.props = props
        assert(len(props)==self.nProps),(
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
        abaqusMaterialName = 'Material-1'
        materialCardPath = Path(writePath).joinpath('material.inp')
        with open(materialCardPath,'w') as fp:
            fp.write('*Material, name={}\n'.format(abaqusMaterialName))
            fp.write('*Density\n')
            fp.write('{:8},\n'.format(self.density))
            fp.write('*Depvar\n')
            fp.write('{:8},\n'.format(self.nStatev))
            fp.write('*User Material, constants={}\n'.format(self.nProps))
            for i in range(self.nProps):
                if (i+1)%8 == 0:
                    fp.write('{:8}'.format(self.props[i]))
                    if i!=self.nProps:
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
    def __init__(self,name,inputfile,solver,material):
        # Calls the base (super) class's constructor
        super().__init__(name)
        # Saves the path to the input file provided as an absolute path
        self.inputfile = Path(__file__).parent.joinpath(inputfile)
        assert (self.inputfile.exists()),(
            'The inputfile provided could not be found: '+str(self.inputfile))
        # Saves the solver to use in Abaqus
        self.solver = solver
        assert (self.solver == AbaqusSolver.Explicit or self.solver == AbaqusSolver.Implicit),(
            'Unknown solver: '+str(self.solver))
        # Saves the material to be used
        self.material = material
        assert (isinstance(material,Material)),(
            'Unknown material: '+str(self.material))
        # Sets up the Abaqus folder path and the test working directory path
        self.abaqusPath = Path(__file__).parent.joinpath('Abaqus')
        self.testPath = self.abaqusPath.joinpath(self.name).joinpath(
            self.solver.name).joinpath(self.material.name)

    # Sets up the Abaqus job files
    def SetupAbaqusJob(self):
        if not self.testPath.exists():
            self.testPath.mkdir(parents=True)
        # Copies the Abaqus input file to the test working directory
        shutil.copy(self.inputfile,self.testPath.joinpath(self.name+'.inp'))
        #
        # Writes a material card to the test working directory
        self.material.WriteMaterailCard(self.testPath)

    # Runs the test on Snurre. Note that the script should be run from Snurre
    def RunOnSnurre(self):
        # Sets up the test
        self.SetupAbaqusJob()
        # Reads the template jobaba file and writes a jobaba file to the test working directory
        jobabaTemplateFile = self.abaqusPath.joinpath('Snurre-jobaba')
        jobabaContent = jobabaTemplateFile.read_text().replace('<<jobName>>',self.name)
        jobabaName = 'jobaba'
        jobabaFile = self.testPath.joinpath(jobabaName)
        jobabaFile.write_text(jobabaContent)
        # Submit the job to the queue on Snurre
        os.system('qsub '+str(jobabaFile))
        # print('qsub '+str(jobabaFile))


    # Runs the test on the current PC
    def Run(self):
        # TODO: Implement a method(this) to run the test on the current PC
        # Sets up the test
        self.SetupAbaqusJob()
##----------------------------------------------------------------------
class FortranTest(Test):
    # Constructor
    def __init__(self,name):
        # Calls the base (super) class's constructor
        super().__init__(name)
    # TODO: Implement tests to test parts of the Fortran code
##----------------------------------------------------------------------
## Main
##----------------------------------------------------------------------
def main():
    # Material density
    density = 2.7e-9
    # Euler angles to be used
    eulerAngles = [EulerAngles(0.0,0.0,0.0),
                    EulerAngles(45.0,0.0,0.0),
                    EulerAngles(129.2315,114.0948,333.4349),
                    EulerAngles(206.5651,114.0948,50.76848)]
    
    # Creates Kalidindi materials
    kalidindiMaterialNames = ['000','4500','inverted','other']
    kalidindiMaterials = []
    for i in range(len(kalidindiMaterialNames)):
        kalidindiMaterials.append(
            Material(kalidindiMaterialNames[i],density,
            [106430., 60350., 28210., 0.01, 0.005, 46.7301, 1.4, 1.,
            eulerAngles[i].phi1,eulerAngles[i].PHI,eulerAngles[i].phi2, 2., 411.256, 104.029, 1.35459, 0.,
                      1.]))
    
    # Creates Voce hardening materials
    voceMaterialNames = ['000 - Voce','4500 - Voce','inverted - Voce','other - Voce']
    voceMaterials = []
    for i in range(len(voceMaterialNames)):
        voceMaterials.append(
            Material(voceMaterialNames[i],density,
            [106430.,  60350., 28210., 0.01, 0.005, 46.7301, 1.4, 1.,
            eulerAngles[i].phi1,eulerAngles[i].PHI,eulerAngles[i].phi2, 1., 20.48, 18.07, 157.3, 39.11,
                      1.]))
    
    # Add different tests to be run
    tests = []
    # Add SimpleShear tests using Abaqus/Explicit and the kalidindi materials
    for material in kalidindiMaterials:
        tests.append(AbaqusTest('SimpleShear','Abaqus/SimpleShear-Explicit.inp',
                        AbaqusSolver.Explicit,material))
    # Add SimpleShear tests using Abaqus/Explicit and the Voce hardening materials
    for material in voceMaterials:
        tests.append(AbaqusTest('SimpleShear','Abaqus/SimpleShear-Explicit.inp',
                        AbaqusSolver.Explicit,material))
    # Add SimpleShear tests using Abaqus/Standard and the kalidindi materials
    for material in kalidindiMaterials:
        tests.append(AbaqusTest('SimpleShear','Abaqus/SimpleShear-Implicit.inp',
                        AbaqusSolver.Implicit,material))
    
    # Run Tests
    for test in tests:
        test.RunOnSnurre()
##----------------------------------------------------------------------
## Entry point
##----------------------------------------------------------------------
if __name__=='__main__':
    main()
##----------------------------------------------------------------------
## Import libraries
##----------------------------------------------------------------------
import odbAccess
from sys import argv, exit
import numpy as np


##----------------------------------------------------------------------
## Load node data
##----------------------------------------------------------------------
def LoadNodeOutput(step, string):
    Sets = step.historyRegions.keys()
    S_temp = []
    N = 0.0
    for Set in Sets:
        if ('Node' in Set) and (string in step.historyRegions[Set].historyOutputs.keys()):
            S_temp = np.array(step.historyRegions[Set].historyOutputs[string].data)
            N += 1.0
            if (N == 1.0):
                S = S_temp[:, 1]
            else:
                pass

    return [np.nan_to_num(S), N]


##----------------------------------------------------------------------
## main
##----------------------------------------------------------------------
def main():
    # Open ODB-FILE
    odbName = 'Axisymmetric.odb'
    odb = odbAccess.openOdb(path=odbName)
    step = odb.steps['Load']

    # Load node data
    [RF, N1] = LoadNodeOutput(step, 'RF2')
    [U, N2] = LoadNodeOutput(step, 'U2')

    # Close odb
    odb.close()

    # Check if the number of nodes extracted from is 1
    if N1 != 1.0 or N2 != 1.0:
        print 'Error in: '+__file__
        exit(1)

    # Write force displacement data to the result file
    with open('Result.csv', 'w') as fil:
        fil.write('%20s , %20s\n' % ('Displacement', 'Force'))
        for disp, force in zip(U, RF):
            fil.write('%20.8e , %20.8f\n' % (disp, force))


##----------------------------------------------------------------------
## Entry point
##----------------------------------------------------------------------
if __name__ == '__main__':
    main()
    exit(0)

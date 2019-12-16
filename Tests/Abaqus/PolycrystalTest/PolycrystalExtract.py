##----------------------------------------------------------------------
## Import libraries
##----------------------------------------------------------------------
import odbAccess
from sys import argv,exit
from math import *
##----------------------------------------------------------------------
## main
##----------------------------------------------------------------------
def main():
    # Open ODB-file
    odbName = 'Polycrystal.odb'
    odb = odbAccess.openOdb(path=odbName)
    step = odb.steps['Step-1']

    Sets = step.historyRegions.keys()

    # Force in y-direction
    RF1 = []
    U1  = []
    U2  = []
    U3  = []

    # Extract data from ODB
    for set in Sets:
        if ('Node' in set) and ('RF1' in step.historyRegions[set].historyOutputs.keys()):
            RF1.append(step.historyRegions[set].historyOutputs['RF1'].data)
        if ('Node' in set) and ('U1' in step.historyRegions[set].historyOutputs.keys()):
            U1 = step.historyRegions[set].historyOutputs['U1'].data
        if ('Node' in set) and ('U2' in step.historyRegions[set].historyOutputs.keys()):
            U2 = step.historyRegions[set].historyOutputs['U2'].data
        if ('Node' in set) and ('U3' in step.historyRegions[set].historyOutputs.keys()):
            U3 = step.historyRegions[set].historyOutputs['U3'].data

    RF = []
    for i in range(len(RF1)):
        A = []
        for j in range(len(RF1[i])):
            A.append(RF1[i][j][1])
        RF.append(A)

    F = [sum([enAvListene[i] for enAvListene in RF]) for i in range(len(RF[0]))]

    # Write log. strain and true stress to the result file
    with open('Result.csv', 'w') as fil:
        fil.write('%20s %20s\n' % ('log. strain', 'true stress'))
        for i in range(len(U1)):
            u = U1[i][1]
            u_w = U2[i][1]
            u_t = U3[i][1]
            epsilon = log(1.0 + u)
            sigma = -F[i]/((1.0+u_t)*(1.0+u_w))
            fil.write('%20.8f %20.8f\n' % (epsilon, sigma))
    # Close odb
    odb.close()
##----------------------------------------------------------------------
## Entry point
##----------------------------------------------------------------------
if __name__=='__main__':
    main()
    exit(0)

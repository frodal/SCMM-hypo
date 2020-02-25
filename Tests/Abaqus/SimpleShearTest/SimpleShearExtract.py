##----------------------------------------------------------------------
## Import libraries
##----------------------------------------------------------------------
import odbAccess
from sys import argv,exit
import numpy as np
##----------------------------------------------------------------------
## Load element data
##----------------------------------------------------------------------
def LoadElement(step,str):
	Sets = step.historyRegions.keys()
	S_temp=[]
	N=0.0
	for set in Sets:
		if ('Element' in set) and (str in step.historyRegions[set].historyOutputs.keys()):
			S_temp=np.array(step.historyRegions[set].historyOutputs[str].data)
			N+=1.0
			if (N==1.0):
				S=S_temp[:,1]
			else:
				S+=S_temp[:,1]

	return [np.nan_to_num(S),N]
##----------------------------------------------------------------------
## main
##----------------------------------------------------------------------
def main():
    # Open ODB-FILE 
    odbName = 'SimpleShear.odb'
    odb = odbAccess.openOdb(path=odbName)
    step = odb.steps['Step-1']

    # Load SDVs from element
    [SDV26,N] = LoadElement(step,'SDV26')
    [SDV27,N] = LoadElement(step,'SDV27')

    # Close odb
    odb.close()

    # Check if the number of integration points extracted from is 1
    if N!=1.0:
        print 'Error in: '+__file__
        exit(1)

    # Write equivalent von Mises plastic strain and stress to the result file
    with open('Result.csv','w') as fil:
        fil.write('%20s , %20s\n' % ('VM eq. pl. strain','VM eq. stress'))
        for eps,sig in zip(SDV27,SDV26):
	        fil.write('%20.8e , %20.8f\n' % (eps,sig))

##----------------------------------------------------------------------
## Entry point
##----------------------------------------------------------------------
if __name__=='__main__':
    main()
    exit(0)
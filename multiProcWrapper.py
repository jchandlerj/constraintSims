########################################################
#NAME: multi_sim_wrapper.py
#CREATED: XXX
#MODIFIED: XXX
#
#DESCRIPTION: Program is a wrapper that calls n instances of the core simulation, 
# where n is equal to: (number of cores on the machine) - 1.  
#
#NOTES:
#	1/ Tested in Windows 7 64-bit environments.
########################################################

########################################################
#STEP 0: Import python modules, functions, and global variables
########################################################

import multiprocessing
from constraintSim import constraintSim
#import sys, os	

########################################################
#WRAPPER LOOP
#DESC: Small loop that initiates a number of processes equal to the number of cores
########################################################

num_cores = multiprocessing.cpu_count() - 1
print(num_cores)

if __name__ == '__main__':
 jobs = []
 for i in range(num_cores):
  #print('allocating jobs to cores')
  p = multiprocessing.Process(target=constraintSim, args=(i,))
  jobs.append(p)
  p.start()


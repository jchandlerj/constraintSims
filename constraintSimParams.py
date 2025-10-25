###############################
#Name: constraintSimParams.py
#Created by: XXX
#Created: XXX
#Desc: Program generates parameters for Barabasi Albert (BA), Erdos Renyi (ER),
 #Watts Strogatz (SW), Holme Kim (HK), or Herrera Zufiria (HZ) random network.
#Used by: constraintSim.py
#Notes:
 #1/ All configurable parameters exist in step 0
 #2/ numSimNets refers to the number of simulations per core. Would ideally be number
  #in total (using math.floor to ensure no overlapping sim IDs).
 #3/ Based on Strat. Sci. reviewer comments, simulating small networks, as these tend to correspond to the examples people offer of where 
 #      constriant is perceived as working correctly. Adjusted min/max NetSize to simulate small nets.
###############################

import random

###############################
#STEP 0: Configurable parameters
###############################

numSimNets = 300 #number per core
minAvgDegree = 2
minDensity = 0.0
maxDensity = 1.0
#Changing range from [30, 400] to [5, 29], to explore small networks. 
minNetSize = 5# 30
maxNetSize = 30 #400
HZwalkLen = random.randint(1, 15) #HZ (2011, 4) used 7
min1StepP = 0.0 #used to generate cc in HZ model
max1StepP = 1.0 #used to generate cc in HZ model
minRewireP = 0.0 #used in assign P in SW model
maxRewireP = 0.5 #used in assign P in SW model
pTF = random.uniform(0.0, 0.30) #HK paper seems to have used 0.15. 
tsExponent = 2.0 #exponent applied to tie strengths in relevant tsMethod routines
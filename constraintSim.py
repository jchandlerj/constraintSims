###############################
#Name: constraintSim.py
#Created by: XXX
#Created: XXX
#Desc: Program repeatedly generates a random network using one of 5 network generators.
 #Program then computes Burt's constraint, and decomposes that into underlying indicators
 #of network properties. Program uses built-in Python functions to compute some network
 #properties, and writes network, edge, and ego-level data out to SQL database.
#Depends on: 
 #constraintSimParams.py (cSP)
 #graphGen.py
 #constraintDecomp.py
#Used by: multiProcWrapper.py (cD)
#Notes:
 #1/ Configurable parameters exist in constraintSimParams.py
 #2/ 'instance_num' is an integer associated with the core on which a given instance 
  #of constraintSim is running. At any time, (num_cores - 1) instances can be running.
  #instance_num ensures that each sim_id is unique.
 #3 Limited random graph type to HZ and HK, since won't be using others.
 #4 Can avoid saving edgelist data to save time.
###############################

###############################
#STEP 0: Import modules/functions
###############################

import constraintSimParams as cSP
import imp
import math
import random
import networkx as nx
import numpy as np
import scipy as sp
import decimal
import datetime
import pyodbc
import pandas as pd

from graphGen import graphGen
from constraintDecomp import constraintDecomp as cD

#delete after testing
#instance_num = 0

###############################
#Define constraintSim as while loop
###############################

def constraintSim(instance_num = 0):

 ###############################
 #STEP 0: Establish ODBC connection and ensure DB tables exist
 ###############################


 #Specify connection
 cnxn = pyodbc.connect("DSN=ConstraintSim")

 c = cnxn.cursor()

 #Create tables in which to store graph, ego, and edge output
 c.execute ('''if not exists (select 1 from INFORMATION_SCHEMA.TABLES where TABLE_NAME ='graphSum') 
  create table graphSum (sim_id int, runTime real, netType nvarchar(16), rewireP real, netSize int, 
  densTarget real, densActual real, avgCC real, transitivity real, walkLen int,
  cc real, linkAdd int, avgDegree real, pTF real, tsMethod nvarchar(16), tsExp real, symmetric int,

  Ci_DD real, Ci_sizeEffect real, Ci_varEffect real, Ci_ID real, Ci_TB real,
  Ci_QS real, Ci_OQD real, Ci_CQD real, Ci_betweenness real, Ci_clustering real,
  Ci_degree real, DD_sizeEffect real, DD_varEffect real, DD_ID real, DD_TB real,
  DD_QS real, DD_OQD real, DD_CQD real, DD_betweenness real, DD_clustering real,
  DD_degree real, sizeEffect_varEffect real, sizeEffect_ID real, sizeEffect_TB real,
  sizeEffect_QS real, sizeEffect_OQD real, sizeEffect_CQD real, sizeEffect_betweenness real,
  sizeEffect_clustering real, sizeEffect_degree real, varEffect_ID real, varEffect_TB real,
  varEffect_QS real, varEffect_OQD real, varEffect_CQD real, varEffect_betweenness real,
  varEffect_clustering real, varEffect_degree real, ID_TB real, ID_QS real, 
  ID_OQD real, ID_CQD real, ID_betweenness real, ID_clustering real, ID_degree real,
  TB_QS real, TB_OQD real, TB_CQD real, TB_betweenness real, TB_clustering real,
  TB_degree real, QS_OQD real, QS_CQD real, QS_betweenness real, QS_clustering real,
  QS_degree real, OQD_CQD real, OQD_betweenness real, OQD_clustering real, 
  OQD_degree real, CQD_betweenness real, CQD_clustering real, CQD_degree real,
  betweenness_clustering real, betweenness_degree real, clustering_degree real,
  betweenness_C_net_size real, betweenness_C_net_var real, betweenness_C_net_DD real,
  clustering_C_net_size real, clustering_C_net_var real, clustering_C_net_DD real

  )''')

 c.execute ('''if not exists (select 1 from INFORMATION_SCHEMA.TABLES where TABLE_NAME ='ego_time') 
  create table ego_time (sim_id int, time_id int, ego_id int, concentration real NULL, 
  output real, input real, degCent int, Ci real, DD real, varTS real, sqAvgTS real, TD real, 
  ID real, CQD real, OQD real, betweenness real, clustering real, sizeEffect real, varEffect real)''')

 c.execute ('''if not exists (select 1 from INFORMATION_SCHEMA.TABLES where TABLE_NAME ='edgelist_time') 
  create table edgelist_time (sim_id int, time_id int, ego_id int, alter_id int, tieStrength real, pij real,
  frequency int, aggIndirect real)''')

 cnxn.commit()

 sim_id = instance_num*cSP.numSimNets
 while sim_id < cSP.numSimNets*(instance_num + 1):

  startTime = datetime.datetime.now()

  ###############################
  #STEP 1: Generate random network from randomized parameters
  ###############################

  #generate random network parameters
  netSize = random.randint(cSP.minNetSize, cSP.maxNetSize)
  netDensity = random.uniform(cSP.minDensity, cSP.maxDensity)
  rewireP = random.uniform(cSP.minRewireP, cSP.maxRewireP)
  #Force density such that average degree >= 2
  if netDensity*(netSize-1) < cSP.minAvgDegree:
   netDensity = int(cSP.minAvgDegree)/(netSize-1)  
  avgDegree = (netSize - 1)*netDensity
  linkAdd = max(1, int(avgDegree/2)) #parameter for the BA, HK, and HZ generators 
   #~1.3 seems to be best denominator, but HZ use 2. This produces too few ties 
   #because random walk with length 2 can return to original node, those "marking" 
   #twice, and only producing 1 new tie. 
  walkLen = cSP.HZwalkLen
  cc = random.uniform(cSP.min1StepP, cSP.max1StepP)

  #randomly choose network type, tie strength method, and symmetry status
  netType = random.choice(['ER', 'BA', 'SW', 'HK', 'HZ'])
  #netType = random.choice(['HK', 'HZ'])
  tsMethod = random.choice(['equal', 'freq', 'freqExp', 'rand', 'randExp', 'revRandExp'])
  symmetric = random.choice([0, 1])

  print('initialized graph params')

  #call function to generate network
  net = graphGen([netType, rewireP, netSize, netDensity, walkLen, cc, 
   linkAdd, avgDegree, cSP.pTF, tsMethod, cSP.tsExponent, symmetric])

  print('created graph')    
    
  ###############################
  #STEP 2: Decompose constraint for random net and compute correlations on returned parameters.
  ###############################
  
  #call function to process net and decompose constraint
  net = cD(net)

  print('decomposed graph')
  ###############################
  #STEP 3: Write data out to DB
  ###############################

  #Write ego data
  for i in net:
   cur_values = (sim_id, 1, i, net.node[i]['conc'], net.node[i]['output'], net.node[i]['input'], 
    net.node[i]['degree'], net.node[i]['Ci'], net.node[i]['DD'], net.node[i]['varTS'], net.node[i]['sqAvgTS'],
    net.node[i]['TB'], net.node[i]['ID'], net.node[i]['CC'], net.node[i]['IR'], net.node[i]['betweenness'], 
    net.node[i]['clustering'], net.node[i]['sizeEffect'], net.node[i]['varEffect']) 
   #print (cur_values)
   c.execute('INSERT INTO ego_time VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', cur_values)
   cnxn.commit()

  #Write edgelist data
  #for i in net.node:  
  # for j in net.edge[i]:
  #  cur_values = (sim_id, 1, i, j, net.edge[i][j]['weight'], net.edge[i][j]['pij'], net.edge[i][j]['freq'], 
  #   net.edge[i][j]['aggIndirect']) 
  #  #print (cur_values)
  #  c.execute('INSERT INTO edgelist_time VALUES (?, ?, ?, ?, ?, ?, ?, ?)', cur_values)
  #  cnxn.commit()

  #approximate total processing time for this graph
  endTime = datetime.datetime.now()
  et = endTime - startTime

  #Write Graph data
  cur_values = (sim_id, et.total_seconds(), net.graph['netType'], net.graph['rewireP'], net.graph['netSize'], 
   net.graph['netDensity'], nx.density(net), net.graph['avgCC'], net.graph['transitivity'], 
   net.graph['walkLen'], net.graph['cc'], net.graph['linkAdd'], 
   net.graph['avgDegree'], net.graph['pTF'], net.graph['tsMethod'], net.graph['tsExp'], net.graph['symmetric'],
   net.graph['Ci_DD'], net.graph['Ci_sizeEffect'], net.graph['Ci_varEffect'], net.graph['Ci_ID'], net.graph['Ci_TB'],
   net.graph['Ci_QS'], net.graph['Ci_OQD'], net.graph['Ci_CQD'], net.graph['Ci_betweenness'], net.graph['Ci_clustering'],
   net.graph['Ci_degree'], net.graph['DD_sizeEffect'], net.graph['DD_varEffect'], net.graph['DD_ID'], net.graph['DD_TB'],
   net.graph['DD_QS'], net.graph['DD_OQD'], net.graph['DD_CQD'], net.graph['DD_betweenness'], net.graph['DD_clustering'],
   net.graph['DD_degree'], net.graph['sizeEffect_varEffect'], net.graph['sizeEffect_ID'], net.graph['sizeEffect_TB'],
   net.graph['sizeEffect_QS'], net.graph['sizeEffect_OQD'], net.graph['sizeEffect_CQD'], net.graph['sizeEffect_betweenness'],
   net.graph['sizeEffect_clustering'], net.graph['sizeEffect_degree'], net.graph['varEffect_ID'], net.graph['varEffect_TB'],
   net.graph['varEffect_QS'], net.graph['varEffect_OQD'], net.graph['varEffect_CQD'], net.graph['varEffect_betweenness'],
   net.graph['varEffect_clustering'], net.graph['varEffect_degree'], net.graph['ID_TB'], net.graph['ID_QS'], 
   net.graph['ID_OQD'], net.graph['ID_CQD'], net.graph['ID_betweenness'], net.graph['ID_clustering'], net.graph['ID_degree'],
   net.graph['TB_QS'], net.graph['TB_OQD'], net.graph['TB_CQD'], net.graph['TB_betweenness'], net.graph['TB_clustering'],
   net.graph['TB_degree'], net.graph['QS_OQD'], net.graph['QS_CQD'], net.graph['QS_betweenness'], net.graph['QS_clustering'],
   net.graph['QS_degree'], net.graph['OQD_CQD'], net.graph['OQD_betweenness'], net.graph['OQD_clustering'], 
   net.graph['OQD_degree'], net.graph['CQD_betweenness'], net.graph['CQD_clustering'], net.graph['CQD_degree'],
   net.graph['betweenness_clustering'], net.graph['betweenness_degree'], net.graph['clustering_degree'],
   net.graph['betweenness_C_net_size'], net.graph['betweenness_C_net_var'], net.graph['betweenness_C_net_DD'],
   net.graph['clustering_C_net_size'], net.graph['clustering_C_net_var'], net.graph['clustering_C_net_DD'],)
 
  #print (cur_values)

  c.execute('INSERT INTO graphSum VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', cur_values)

  cnxn.commit()
  
  print('inserted data into DB')

  ############################## 
  # round(decimal.Decimal(nx.average_clustering(net)),2), 
  # "trans is: ", round(decimal.Decimal(nx.transitivity(net)),2), 
  # "tsM is: ", tsMethod)

  #increment to instance's next network sim
  sim_id += 1

  #Generate new simulation parameters
  imp.reload(cSP)



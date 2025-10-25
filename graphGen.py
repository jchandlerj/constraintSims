###############################
#Name: graphGen.py
#Created by: XXX
#Created: XXX
#Desc: Program generates a single random graph from one of 5 methods: Erdos Renyi (ER),
 #Watt Strogatz (SW), Barabasi Allen (BA), Holme Kim (HK), and Herrera Zufiria (HZ). HZ
 #graphs are not a built-in Python function, but instead generated via a custom function.
 #After generating a random graph, graphGen assigns tie strengths. The resulting graph
 #is returned.
#Depends on: 
 #constraintSimParams.py
 #graphGenHerreraZufiria.py
#Used by: constraintSim.py
#Notes:
 #1/ Configurable parameters exist in constraintSimParams.py
 #2/ SW networks may not be connected, so using connected_small_world. This repeatedly executes
  #SW until a connected graph occurs or a max number of iterations passes. Thus, it may still
  #yield a disconnected graph. If an isolate occurs, and if the isolate is chosen in the random
  #walk, this can yield an error.
 #3/ Added an "if" clause to random walk to avoid (2)
 #4/ May want to make random walk a separate function, and only call if using 'freq'-dependent
  #tie strength method.
 #5/ Added concentration attribute when computing tie strength (because already in that loop)
###############################

###############################
#Import modules/functions
###############################

import math
import random
import networkx as nx
import decimal

from graphGenHerreraZufiria import graphGenHerreraZufiria as genHZ

###############################
#Define graphGen to generate 1 random graph using constraintSim parameters
###############################

def graphGen(input):

 netType = input[0]
 rewireP = input[1]
 netSize = input[2]
 netDensity = input[3]
 walkLen = input[4]
 cc = input[5]
 linkAdd = input[6]
 avgDegree = input[7]
 pTF = input[8]
 tsMethod = input[9]
 tsExp = input[10]
 symmetric = input[11]
  
 #Generate random network
 if netType == 'ER':
  net = nx.erdos_renyi_graph(netSize, netDensity)
 elif netType == 'BA':
  net = nx.barabasi_albert_graph(netSize, int(linkAdd), seed=None)
 #used connected SW to avoid isolates in subsequent frequency random walks. To
  #reduce number of cases where isolates can occur, increased tries from default (100)
  #to 1000.
 elif netType == 'SW':
  net = nx.connected_watts_strogatz_graph(netSize, int(avgDegree), rewireP, 1000)
 elif netType == 'HK':
  net = nx.powerlaw_cluster_graph(netSize, int(linkAdd), pTF)
 else:
  net = genHZ([netSize, netDensity, walkLen, cc, linkAdd])

 #store the graph properties
 net.graph['netType'] = input[0]
 net.graph['rewireP'] = input[1]
 net.graph['netSize'] = input[2]
 net.graph['netDensity'] = input[3]
 net.graph['walkLen'] = input[4]
 net.graph['cc'] = input[5]
 net.graph['linkAdd'] = input[6]
 net.graph['avgDegree'] = input[7]
 net.graph['pTF'] = input[8]
 net.graph['tsMethod'] = input[9]
 net.graph['tsExp'] = input[10]
 net.graph['symmetric'] = input[11]

 #For non-HZ methods, create freq attribute and use random walks to populate (only necessary 
  #to influence tie strength). HZ networks already have this attribute.
 if netType != 'HZ':
  #add freq attribute
  '''for i in net.node:
   for j in net.edge[i]:
    net.edge[i][j]['freq'] = 1'''
  for i in net:
   for j in net[i]:
    net[i][j]['freq'] = 1
  walks = 0
  while walks < netSize:

   # choose random starting point
   ego = random.choice(list(net.nodes()))
   l = 0

   #conduct random walk from random starting point. Because isolates can occur
    #from SW rewiring, only execute random walk if length of ego's alters > 0
   #if any(net.neighbors(ego)):
   #deprecated this old line
   if len(net[ego]) > 0:
    while l < walkLen:
     alter = random.choice(list(net[ego]))
     #store fact that edge has been traversed
     net[ego][alter]['freq'] = net[ego][alter]['freq']+1
     ego = alter
     l += 1
 
   #increment population counter
   walks += 1

 #Calculate some graph properties.
 net.graph['avgCC'] = float(decimal.Decimal(nx.average_clustering(net)))
 net.graph['transitivity'] = float(decimal.Decimal(nx.transitivity(net)))

 #add node-level attributes: betweenness and clustering coefficient. Have to do this before
  #setting graph to diGraph.
 cl = nx.clustering(net)
 bet = nx.betweenness_centrality(net)
 nx.set_node_attributes(net, name = 'clustering', values = cl)
 nx.set_node_attributes(net, name = 'betweenness', values = bet)

 #set to diGraph if net is supposed to be asymmetric. This will allow asymmetric tie weights
  #in the subsequent step.
 if symmetric == 0:
  net = nx.DiGraph(net)

 #create and populate tie strength (weight) attribute. Also use this loop to create a 
  #'concentration' node-level attribute
 #nx.set_node_attributes(net, 1.0, name='conc')
 for i in net:
  net.nodes[i]['conc'] = 1.0
  for j in net[i]:
   if tsMethod == 'equal':
    net.edges[i, j]['weight'] = float(1)
   elif tsMethod == 'freq':
    net.edges[i, j]['weight'] = float(net.edges[i, j]['freq'])
   elif tsMethod == 'freqExp':
    net.edges[i, j]['weight'] = math.pow(float(net.edges[i, j]['freq']), tsExp)  
   elif tsMethod == 'rand':
    net.edges[i, j]['weight'] = random.uniform(0, 1)
   elif tsMethod == 'randExp':
    net.edges[i, j]['weight'] = math.pow(random.uniform(0, 1), tsExp)
   elif tsMethod == 'revRandExp':
    net.edges[i, j]['weight'] = math.pow(1-random.uniform(0, 1), tsExp)

 #return graph
 return(net)
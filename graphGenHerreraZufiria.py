###############################
#Name: graphGenHerreraZufiria.py
#Created by: XXX
#Created: XXX
#Notes:
# 1/ Implements Herrera and Zufiria, 2011 (HF)
# 2/ The 2-step move in STEP 3 often returns to the origination point. This may be OK,
#  but seems like it may reduce avg(degree) below target. observed avg(degree) and 
#  density ~50% desired. Could change routine to avoid doubling back, but would need 
#  to ensure no infinite loops (random walk goes to end of branch and cannot return).
# 3/ At edge creation, set edge traffic ('freq') = 1.
###############################

###############################
#STEP 0: Import necessary libraries
###############################

import networkx as nx
import random

###############################
#STEP 1: Create variables
###############################

def graphGenHerreraZufiria(input):

 #create configurable variables
 netSize = input[0] #random.randint(40, 500)
 netDensity = input[1] #random.uniform(.1, .5)
 walkLen = input[2] #7
 cc = input[3] #random.uniform(0, 1)
 linkAdd = input[4]

 #ensure a sufficient starting size (HZ, 2011, 4) 
 if linkAdd > 10:
  initSize = linkAdd
 else:
  initSize = 10 

 ###############################
 #STEP 2: Populate base graph
 ###############################

 #create base ring graph of size initSize (HZ, 2011, 4)
 G = nx.cycle_graph(initSize)

 #add node attributes
 for i in G:
  #create lifespan in case we eventually want node death
  G.nodes[i]['lifespan'] = 0
  G.nodes[i]['marked'] = 0
  #implement P(triad formation) based on the proporition cc (HZ, 2011, 3)
  if random.uniform(0, 1) < cc:
   G.nodes[i]['pTF'] = 1.0
  else:
   G.nodes[i]['pTF'] = 0.0

 #add information about edge utilization (number of times edge is traversed
  #during network generation). Set initial traffic = 1. Otherwise, if tie is 
  #never traversed, end up with existing ties with potential weight = 0
 for i in G:
  for j in G[i]:
   G.edges[i, j]['freq'] = 1

 ###############################
 #STEP 3: Run graph generator, adding new nodes and edges
 ###############################

 pop = initSize
 while pop < netSize:

  #print ("NEW NODE ADDITION ROUTINE")

  # choose random starting point
  ego = random.choice(list(G.nodes()))
  l = 0

  #conduct random walk from random starting point
  while l < walkLen:
   alter = random.choice(list(G[ego]))
   #alter = random.choice(list(G[ego].keys()))
   #alter = random.choice(list(G.neighbors(ego))) 
   #store fact that edge has been traversed
   G.edges[ego, alter]['freq'] = G.edges[ego, alter]['freq']+1
   ego = alter
   l += 1

  #mark arrival node for subsequent tie creation
  G.nodes[ego]['marked'] = 1
  #print ("random walk start node is ", ego)

  #random walk from arrival node to mark remaining alters
  m = 0
  while m < linkAdd-1:
   #determine whether 1 or 2 step path
   if random.uniform(0, 1) < G.nodes[ego]['pTF']:
    #move 1 step
    #alter = random.choice(list(G.edge[ego]))
    alter = random.choice(list(G[ego]))
    G.edges[ego, alter]['freq'] = G.edges[ego, alter]['freq']+1
    #print ("   one step search; next node is ", alter)
   else:
    #move 2 steps
    #alter = random.choice(list(G.edge[ego]))
    alter = random.choice(list(G[ego]))
    G.edges[ego, alter]['freq'] = G.edges[ego, alter]['freq']+1
    #print ("   two step search; first ego is ", alter)
    ego = alter
    #alter = random.choice(list(G.edge[ego]))
    alter = random.choice(list(G[ego]))
    G.edges[ego, alter]['freq'] = G.edges[ego, alter]['freq']+1
    #print ("   two step search; second ego is ", alter)
   ego = alter
   G.nodes[ego]['marked'] = 1
   #increment to next link
   m += 1 

  #add node, and add edges with marked nodes
  if random.uniform(0, 1) < cc:
   p = 1.0
  else:
   p = 0.0
  G.add_node(pop, freq = 0, marked = 0, pTF = p)
  ego = pop
  #print ("adding node ", ego) 

  for i in G:
   if G.nodes[i]['marked'] == 1:
    G.add_edge(ego, i, freq = 1)
    #print ("adding an edge with ", ego, "as ego and alter equal ", i) 

  #reset all nodes' 'marked' values
  for i in G:
   G.nodes[i]['marked'] = 0
 
  #increment population counter
  pop += 1

 #print (" running genHZ. cc = ", cc)
 #print (nx.density(G))
 return(G)
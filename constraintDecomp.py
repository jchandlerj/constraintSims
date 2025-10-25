###############################
#Name: constraintDecomp.py
#Created by: XXX
#Created: XXX
#Desc: Program takes a Python networkx graph as input, ensures that the graph is a DiGraph,
# uses edge characteristics to calculate interim values, and ultimately isolates constraint's
# underlying properties based on Johnson/Sasson's formal decomposition. Function returns the 
# original Python graph object with constraint, constraint-level properties as node-level 
# attributes, and pairwise correlations of constraint's components as graph-level attributes. 
#Depends on: 
#Used by: constraintSim.py
#Notes:
# 1/ For robustness across applications, should add a step to delete self-referencing edges. Or
#  does DiGraph command already do this?
# 2/ Calculations of community closure and indirect redundancy (formerly "triadic redundancy") 
#  avoid duplicate structures (e.g., recomputing (k, q) when already calculated (q, k)) by 
#  only assessing k > q. Tested this with strings, and it should also work. But worth checking. 
# 3/ Calculating CC using: CC = Ci - TB - ID - DD - IR. Commented out code to compute CC 
#  independently, but can use that for verification if later necessary. 
# 4/ Clean up commented-out code
# 5/ Standardize names, perhaps as: direct dependence = (size dependence + hetero. dependence),
#  indirect dependence, triadic dependence, open quadriad dependence, closed quadriad
#  dependence.
# 6/ This does not isolate the Cov(Oj, pij) term, though it could by subtracting size and 
#  variance from DD.
###############################

###############################
#STEP 0: Import modules/functions
###############################

import math
import networkx as nx
import decimal
import numpy as np
import pandas as pd

###############################
#Define constraintDecomp as while loop over nodes in input network object
###############################

def constraintDecomp(net):

 ###############################
 #STEP 1: Prepare graph
 #Ensure graph is DiGraph, add empty edges, and add some node-level attributes
 ###############################

 #force digraph to make sure that we can add any empty edges
 net = nx.DiGraph(net)
 #add empty edges (possible if input graph was based on an asymmetric edgelist): if there's 
  #an i->j edge but no j->i edge, add the j->i edge with weight = 0. 
 for i in net:
  for j in net[i]:
   if i not in net[j]:
    net.add_edge(j, i, weight = 0.0)  

 #add node volume attributes
 for i in net:
  net.nodes[i]['output'] = 0
  net.nodes[i]['input'] = 0

  for j in net[i]:
   #accumulate output
   net.nodes[i]['output'] += net.edges[i, j]['weight']

   #accumulate input
   if net.has_edge(j, i):
    net.nodes[i]['input'] += net.edges[j, i]['weight']

 ###############################
 #STEP 2: Compute all direct elements
 # Add pij values to graph, and compute constraint's dyadic elements 
 #  This happens before Step 3 to ensure that each piq is available when needed. Ordering 
 #  could solve this, but I find this easier to read anyway. Also compute Blau elements 
 #  here.
 ###############################

 #populate pij values. create lists of weights and weights^2, and use these list to compute
  #Blau elements
 for i in net.node:

  #create attribute to hold Blau
  net.node[i]['DD'] = float(0)

  #specify two lists: one for tie weights and one for squared tie weights
  ts = []
  #tsSq = []

  #calculate total activity
  totActivity = net.nodes[i]['input'] + net.nodes[i]['output']

  #update pij and append pij, pij^2 to lists
  for j in net[i]:
   net.edges[i, j]['pij'] = (net.edges[i, j]['weight'] + net.edges[j, i]['weight'])/totActivity
   ts.append(net.edges[i, j]['pij'])
   #tsSq.append(math.pow(net.edge[i][j]['weight'], 2))

   #update Blau
   net.nodes[i]['DD'] += math.pow(net.edges[i, j]['pij'], 2)*net.nodes[j]['conc']

  #update Blau elements, with exceptions for isolates
  #net.node[i]['degree'] = len(ts)
  net.nodes[i]['degree'] = len(net[i])  
  net.nodes[i]['varTS'] = float(0)
  net.nodes[i]['sqAvgTS'] = float(0)
  if net.nodes[i]['degree'] > 1:
   #net.node[i]['sumSqTS'] = sum(tsSq)
   #net.node[i]['sqAvgTS'] = math.pow(np.mean(ts),2)
   #net.node[i]['degree'] = len(ts)
   if net.graph['tsMethod'] == 'equal':
    net.nodes[i]['varTS'] = float(0)
   else:
    net.nodes[i]['varTS'] = np.var(ts)
   #print(np.var(net.edge[i]))

 #Once all nodes updated, add Blau decomposition terms as node attributes
 degree = nx.get_node_attributes(net, 'degree')
 varTS = nx.get_node_attributes(net, 'varTS')

 varEffect = {k : v*degree[k] for k, v in varTS.items() if k in degree}
 #This creates noise, because sizeEffect = when degree = 0. That's consistent with 
  #setting constraint = 1 for isolates, but it would probably be better just to 
  #exclude isolates (i.e., delete nodes from graph where degree == 0).
 sizeEffect = {k: 1/max(v,1) for k, v in degree.items()}

 nx.set_node_attributes(net, varEffect, name = 'varEffect')
 nx.set_node_attributes(net, sizeEffect, name = 'sizeEffect')

 ###############################
 #STEP 3: Compute indirect constraint components.
 # Do this in a series of 3 nested loops, first calculate indirect (sum(piq*pij)), then 
 #  find shared alters to compute triadic stuff, then use additional nested loop to compute
 #  quadratic stuff. 
 ###############################

 for i in net:

  #Create some node attributes. Force constraint = 1 if ego is an isolate
  net.nodes[i]['Ci'] = float(0)
  if net.nodes[i]['degree'] == 0:
   net.nodes[i]['Ci'] = float(1)
  net.nodes[i]['TB'] = float(0)
  #net.node[i]['term3'] = float(0)
  net.nodes[i]['ID'] = float(0)
  net.nodes[i]['QS'] = float(0)
  net.nodes[i]['IR'] = float(0)
  net.nodes[i]['CC'] = float(0)
  #net.node[i]['degree'] = 0

  ###############################
  #3A: Loop for each (i, j) pair
  ###############################

  for j in net[i]:

   #specify edge attributes
   net.edges[i, j]['aggIndirect'] = 0.0
   #net.edge[i][j]['pjqSqConcq'] = 0.0
   #net.edge[i][j]['commClosure'] = 0.0
   #net.edge[i][j]['indRedun'] = 0.0

   ###############################
   #3AI: Evaluate each (i, j, q) triple
   # Find intersection of ego/alter alters, and iterate, calculating sum(piq*pqj). Also
   #  execute nested loop over quadriads.
   ###############################

   #find shared alters. since we've already verified that i->j exists if j-> exists,
    #we're sure to get successors and predecessors for both i and j.
   #convert lists of i's and j's alters to two sets
   ai = set(net.successors(i))
   aj = set(net.successors(j))
   #find intersection and convert to list -- iteration should be faster than over set
   sharedAlters = list(ai&aj)

   for q in sharedAlters:
    net.edges[i, j]['aggIndirect'] += net.edges[i, q]['pij']*net.edges[q, j]['pij']
    net.nodes[i]['ID'] += pow(net.edges[i, j]['pij']*net.edges[j, q]['pij'], 2)*net.nodes[q]['conc']

    ###############################
    #3AIa: Evaluate (i, j, q, k) quads
    # Find communities (CC) and edge-sharing triads (IR)
    ###############################

    #Quadratic stuff
    #find all of q's alters
    aq = set(net.successors(q))
    
    #comment all of this out. instead rely on: QS = CC + IR, and QS = term3 - ID
    ##for mutual alters of j, q, and k; compute community closure
    #comm = list(set(sharedAlters)&aq)
    #comm = [h for h in comm if h > q]
    #for k in comm:
     ##update inner part of dyadic community closure
     #net.node[i]['CC'] += 2*net.node[j]['conc']*net.edge[i][q]['pij']*net.edge[i][k]['pij']*net.edge[q][j]['pij']*net.edge[k][j]['pij']

    #for alters k connected to i and j but not q, compute indirect redundancy
    comp = set(sharedAlters) - aq
    #comp.remove(q)
    comp = list(comp)
    comp = [h for h in comp if h > q]
    for k in comp:
     #update inner part of dyadic indirect redundancy
     net.nodes[i]['IR'] += 2*net.nodes[j]['conc']*net.edges[i, q]['pij']*net.edges[i, k]['pij']*net.edges[q, j]['pij']*net.edges[k, j]['pij']

   ###############################
   #3AII: Increment node[i] indirect-related values based on edge[i][j] 
   ###############################

   net.nodes[i]['Ci'] += math.pow(net.edges[i, j]['pij'] + net.edges[i, j]['aggIndirect'], 2)*net.nodes[j]['conc']
   net.nodes[i]['TB'] += 2*net.edges[i, j]['pij']*net.edges[i, j]['aggIndirect']*net.nodes[j]['conc']
   #net.node[i]['degree'] += 1

   #net.node[i]['term3'] += math.pow(net.edge[i][j]['aggIndirect'], 2)*net.node[j]['conc']
   #net.node[i]['ID'] += net.edge[i][j]['pjqSqConcq']
   #net.node[i]['IR'] += net.edge[i][j]['indRedun']
   #net.node[i]['CC'] += net.edge[i][j]['commClosure']

  ###############################
  #3B: Solve for CC using existing values
  ###############################
 
  #net.node[i]['QS'] = net.node[i]['term3'] - net.node[i]['ID']
  net.nodes[i]['CC'] = net.nodes[i]['Ci'] - (net.nodes[i]['DD'] + net.nodes[i]['TB']  + net.nodes[i]['ID'] + net.nodes[i]['IR'])

 ###############################
 #STEP 4: Generate node attribute correlation matrix, and add correlations
 # as graph attributes
 ###############################

 # Now, using the Ci values, create some other node attributes
 for i in net:
  net.nodes[i]['C_net_size'] = net.nodes[i]['Ci'] - net.nodes[i]['sizeEffect']
  net.nodes[i]['C_net_var'] = net.nodes[i]['Ci'] - net.nodes[i]['varEffect']
  net.nodes[i]['C_net_DD'] = net.nodes[i]['Ci'] - net.nodes[i]['DD']

 #generate and transpose a dataframe with all node attributes
 #df=pd.DataFrame(net.node)
 df=pd.DataFrame.from_dict(dict(net.nodes(data=True)), orient='index')
 dfT=pd.DataFrame.transpose(df)

 #generate the correlation matrix of node attributes. why was this previously transposed?
 #corrs = dfT.corr()
 corrs = df.corr()
 #replace NaN correlations with zeroes
 corrs.fillna(0.0, inplace = True)

 #add correlations as graph attributes
 net.graph['Ci_DD'] = corrs.at["Ci", "DD"]
 net.graph['Ci_sizeEffect'] = corrs.at["Ci", "sizeEffect"]
 net.graph['Ci_varEffect'] = corrs.at["Ci", "varEffect"]
 net.graph['Ci_ID'] = corrs.at["Ci", "ID"]
 net.graph['Ci_TB'] = corrs.at["Ci", "TB"]
 net.graph['Ci_QS'] = corrs.at["Ci", "QS"]
 net.graph['Ci_OQD'] = corrs.at["Ci", "IR"]
 net.graph['Ci_CQD'] = corrs.at["Ci", "CC"]
 net.graph['Ci_betweenness'] = corrs.at["Ci", "betweenness"]
 net.graph['Ci_clustering'] = corrs.at["Ci", "clustering"]
 net.graph['Ci_degree'] = corrs.at["Ci", "degree"]

 net.graph['DD_sizeEffect'] = corrs.at["DD", "sizeEffect"]
 net.graph['DD_varEffect'] = corrs.at["DD", "varEffect"]
 net.graph['DD_ID'] = corrs.at["DD", "ID"]
 net.graph['DD_TB'] = corrs.at["DD", "TB"]
 net.graph['DD_QS'] = corrs.at["DD", "QS"]
 net.graph['DD_OQD'] = corrs.at["DD", "IR"]
 net.graph['DD_CQD'] = corrs.at["DD", "CC"]
 net.graph['DD_betweenness'] = corrs.at["DD", "betweenness"]
 net.graph['DD_clustering'] = corrs.at["DD", "clustering"]
 net.graph['DD_degree'] = corrs.at["DD", "degree"]

 net.graph['sizeEffect_varEffect'] = corrs.at["sizeEffect", "varEffect"]
 net.graph['sizeEffect_ID'] = corrs.at["sizeEffect", "ID"]
 net.graph['sizeEffect_TB'] = corrs.at["sizeEffect", "TB"]
 net.graph['sizeEffect_QS'] = corrs.at["sizeEffect", "QS"]
 net.graph['sizeEffect_OQD'] = corrs.at["sizeEffect", "IR"]
 net.graph['sizeEffect_CQD'] = corrs.at["sizeEffect", "CC"]
 net.graph['sizeEffect_betweenness'] = corrs.at["sizeEffect", "betweenness"]
 net.graph['sizeEffect_clustering'] = corrs.at["sizeEffect", "clustering"]
 net.graph['sizeEffect_degree'] = corrs.at["sizeEffect", "degree"]

 net.graph['varEffect_ID'] = corrs.at["varEffect", "ID"]
 net.graph['varEffect_TB'] = corrs.at["varEffect", "TB"]
 net.graph['varEffect_QS'] = corrs.at["varEffect", "QS"]
 net.graph['varEffect_OQD'] = corrs.at["varEffect", "IR"]
 net.graph['varEffect_CQD'] = corrs.at["varEffect", "CC"]
 net.graph['varEffect_betweenness'] = corrs.at["varEffect", "betweenness"]
 net.graph['varEffect_clustering'] = corrs.at["varEffect", "clustering"]
 net.graph['varEffect_degree'] = corrs.at["varEffect", "degree"]

 net.graph['ID_TB'] = corrs.at["ID", "TB"]
 net.graph['ID_QS'] = corrs.at["ID", "QS"]
 net.graph['ID_OQD'] = corrs.at["ID", "IR"]
 net.graph['ID_CQD'] = corrs.at["ID", "CC"]
 net.graph['ID_betweenness'] = corrs.at["ID", "betweenness"]
 net.graph['ID_clustering'] = corrs.at["ID", "clustering"]
 net.graph['ID_degree'] = corrs.at["ID", "degree"]

 net.graph['TB_QS'] = corrs.at["TB", "QS"]
 net.graph['TB_OQD'] = corrs.at["TB", "IR"]
 net.graph['TB_CQD'] = corrs.at["TB", "CC"]
 net.graph['TB_betweenness'] = corrs.at["TB", "betweenness"]
 net.graph['TB_clustering'] = corrs.at["TB", "clustering"]
 net.graph['TB_degree'] = corrs.at["TB", "degree"]

 net.graph['QS_OQD'] = corrs.at["QS", "IR"]
 net.graph['QS_CQD'] = corrs.at["QS", "CC"]
 net.graph['QS_betweenness'] = corrs.at["QS", "betweenness"]
 net.graph['QS_clustering'] = corrs.at["QS", "clustering"]
 net.graph['QS_degree'] = corrs.at["QS", "degree"]

 net.graph['OQD_CQD'] = corrs.at["IR", "CC"]
 net.graph['OQD_betweenness'] = corrs.at["IR", "betweenness"]
 net.graph['OQD_clustering'] = corrs.at["IR", "clustering"]
 net.graph['OQD_degree'] = corrs.at["IR", "degree"]

 net.graph['CQD_betweenness'] = corrs.at["CC", "betweenness"]
 net.graph['CQD_clustering'] = corrs.at["CC", "clustering"]
 net.graph['CQD_degree'] = corrs.at["CC", "degree"]

 net.graph['betweenness_clustering'] = corrs.at["betweenness", "clustering"]
 net.graph['betweenness_degree'] = corrs.at["betweenness", "degree"]
 net.graph['betweenness_C_net_size'] = corrs.at["betweenness", "C_net_size"]
 net.graph['betweenness_C_net_var'] = corrs.at["betweenness", "C_net_var"]
 net.graph['betweenness_C_net_DD'] = corrs.at["betweenness", "C_net_DD"]
	
 net.graph['clustering_degree'] = corrs.at["clustering", "degree"]
 net.graph['clustering_C_net_size'] = corrs.at["clustering", "C_net_size"]
 net.graph['clustering_C_net_var'] = corrs.at["clustering", "C_net_var"]
 net.graph['clustering_C_net_DD'] = corrs.at["clustering", "C_net_DD"]

 ###############################
 #STEP 5: Return updated graph
 ###############################

 #return updated graph
 return(net)
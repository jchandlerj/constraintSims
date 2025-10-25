# constraintSims
This repo constains all of the .py scripts used in Johnson and Sasson (2025) to study constraint's empirical meaning via synthetic networks. 

#Name: readme.txt
#Desc: Contains info about network constraint simulation programs

These programs generate random graphs of various types, decompose Burt’s (1992) constraint index for each graph, compute metrics across 
the graphs, and store output to a DB using an ODBC connection. 

mutliProcWrapper.py is a wrapper around constraintSim.py, which allows multiple simulations to run simultaneously, horizontally scaling across the host machine’s available cores - 1.

constraintSim.py repeatedly generates random graphs, decomposes those graph’s into constraint’s underlying terms, computes graph-level metrics, and saves ego, edge, and graph-level data to a DB via ODBC. This program relies on:
	-constraintSimParams.py: A program that specifies graph parameterizations
	-graphGen.py: A program that generates a random graph using one of five generators: an Erdos Renyi generator, a small-world
		generator, a scale-free generator, the Holme-Kim SW/SF generator, or the Herrera-Zufiria SW/SF generator. The latter is 		implemented in:
			-graphGenHerreraZufiria.py: a program that implements the HZ (2011) algorithm


The programs rely heavily on Python’s networkx 2.2 library. networkx has changed syntax significantly over this project’s life, so earlier versions of programs will typically not work with newer versions of other programs (e.g., the current version of constraintSim.py will fail if attempting to call an older version of constraintDecomp.py). 

Separately, the authors can make available a generalized version of constraintDecomp. This ingests an (optionally) weighted edgelist, several processing parameters, and (optionally) node-level concentration data. For the provided graph, it returns node-level constraint and constraint’s constituent terms. This program serves as the back-end for a web app with a simple UI.

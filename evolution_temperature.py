# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 15:30:29 2023

@author: Tom Guérinel

Ce programme a pour but de lire les données créées par evolution_mcmc.c
Il sépare les différentes matrices stockées dans le même fichier mat_T.dat
puis les affiches et enregitre l'image.
C'est également ce programme que j'ai utilisé pour tracer les grandeurs en
fonction de la température.
"""

import numpy as np
import matplotlib.pyplot as plt


plt.close('all')



dataaim = np.loadtxt("aimantation_T.dat")
dataT = np.loadtxt("mat_T.dat")
dataen = np.loadtxt("energie_T.dat")
datasus = np.loadtxt("susceptibilite_T.dat")

N = dataT.shape[1]

nb_mat = dataT.shape[0]//N

nd = str(10)

for i in range(nb_mat):
    plt.matshow(dataT[i*N:(i+1)*N,:])
    plt.savefig("evolution_temperature/evo"+nd+"/mat"+str(i))
    plt.close()



plt.figure("aimantation par site")
plt.plot(dataaim[:,0],abs(dataaim[:,1]))
plt.savefig("evolution_temperature/evo"+nd+"/aimantation_par_site")

plt.figure("énergie par site")
plt.plot(dataen[:,0],dataen[:,1])
plt.savefig("evolution_temperature/evo"+nd+"/energie_par_site")

plt.figure("capacité thermique")
capa = np.gradient(dataen[:,1], dataen[:,0])
plt.plot(dataen[:,0],capa)
plt.savefig("evolution_temperature/evo"+nd+"/capacite_thermique")

plt.figure("susceptibilité magnétique")
plt.plot(datasus[2:,0],datasus[2:,1])
plt.savefig("evolution_temperature/evo"+nd+"/susceptibilite_magnetique")
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:45:11 2023


@author: Tom Guérinel


Ce programme a pour but de lire les données créées par evolution_mcmc.c
Il sépare les différentes matrices stockées dans le même fichier matev2.dat
puis les affiches et enregitre l'image.
C'est également ce programme que j'ai utilisé pour tracer l'aimantation au cours
du temps mais cette partie a été surpprimée, elle fonctionnait de la même façon
que la tracé des grandeurs dans le fichier evolution_temperature.py.
"""


import numpy as np
import matplotlib.pyplot as plt


plt.close('all')


dataev = np.loadtxt("matev2.dat")

N = dataev.shape[1]

nb_mat = dataev.shape[0]//N


for i in range(nb_mat):
    plt.matshow(dataev[i*N:(i+1)*N,:])
    plt.savefig("evolution_metro/evo16/mat"+str(i))
    plt.close()



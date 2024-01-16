# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 18:16:40 2023

@author: Tom Guérinel

Ce programme a pour but de lire les données créées par hystresis.c
Il permet de créer les courbes sur lesquelles on visualise l'hystérésis et les
enregistre. Il permet également de tracer le champ coercitif en focntion de la 
température et de le fit.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so


def sign(x):
    if x < 0:
        return -1
    elif x > 0 :
        return 1
    else :
        return 0


plt.close('all')


dataaim_in = np.loadtxt("aimantation_h.dat")
dataaim2_in = np.loadtxt("aimantation_h2.dat")


print(dataaim_in.shape)

N = 65 #nombre de valeurs de champ magnétique extérieur, déterminé dans hysteresis.c

list_aim = []
list_aim2 = []

for i in range(26): #nombre de valeurs de température, déterminé dans hysteresis.c
    list_aim.append(dataaim_in[i*(N+1):(i+1)*(N+1)])
    list_aim2.append(dataaim2_in[i*N:(i+1)*N])




print(list_aim[0])


list_T = []
list_delta = []

for j in range(len(list_aim)): #trace l'aimantation au cours de l'augmentation
                                #et la diminution de h et enregistre chaque image
    
    entete = list_aim[j][0]
    dataaim = list_aim[j][1:]#données de l'augmentation
    dataaim2 = list_aim2[j]#données de la diminution
    
    
    T, nd = entete[0], str(int(entete[1]))

    for i in range(1,len(dataaim)):
        if sign(dataaim[i][1])!=sign(dataaim[i-1][1]):
            chgt_up = (dataaim[i][0]+dataaim[i-1][0])/2 #valeur du champ coercitif
            print(chgt_up)
    
    for i in range(1,len(dataaim2)):
        if sign(dataaim2[i][1])!=sign(dataaim2[i-1][1]):
            chgt_down = (dataaim2[i][0]+dataaim2[i-1][0])/2 #valeur du champ coercitif
            print(chgt_down)
    
    delta = abs(chgt_up-chgt_down)
    
    print(delta, end="\n \n")
    
    list_T.append(T)
    list_delta.append(delta)
    
    
    
    plt.figure("aimantation par site")
    plt.title(f"Aimantation par site en fonction du champ extérieur, T={T}, $\Delta = {delta}$")
    plt.xlabel('h')
    plt.ylabel('M')
    plt.plot(dataaim[:,0],dataaim[:,1], 'r', label = "augmentation")
    plt.plot(dataaim2[:,0],dataaim2[:,1], 'b', label = "diminution")
    plt.plot([chgt_down,chgt_down],[-1,1],"k--")
    plt.plot([chgt_up,chgt_up],[-1,1],"k--")
    plt.legend()
    plt.savefig("evolution_champ_ext/aimantation_par_site_evo"+nd)
    plt.close('all')

list_T = np.array(list_T)
list_delta = np.array(list_delta)

def f(x,a,T_0):
    return a*np.exp(-x/T_0)



#Figure champ coercitif en fonction de la température

plt.figure("hysteresis_T_H_C")
plt.plot(list_T, list_delta/2, 'b',label='Données')
plt.errorbar(list_T, list_delta/2, 0.0625, 0)

popt, pcov = so.curve_fit(f, list_T, list_delta/2) #fit exponentiel
print(popt)
a_H_C,T_0_H_C =popt[0],popt[1]



plt.title("Champ coercitif en fonction de la température")

plt.xlabel('T')
plt.ylabel('$H_C$')

plt.plot(list_T, f(list_T,a_H_C,T_0_H_C),'r', label=f'Fit : {a_H_C:.3f}exp(-T/{T_0_H_C:.3f})')

plt.legend()

print(pcov)



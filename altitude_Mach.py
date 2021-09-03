# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 13:53:56 2021

@author: geigerr
"""

import numpy as np
import math
import matplotlib.pyplot as plt

	
R=1716.55	#ft^2/(sec^2degR)
gamma=1.4
g0=32.17405	#ft/sec^2
RE=20926476	# Radius of the Earth, ft
D = [0,1,2,3,4,5,6,7]
K=[-.00356616,0,.00054864,.00153619,0,-.00109728,-.00219456,0]	#degR/ft
T=[518.67,389.97,389.97,411.57,487.17,487.17,454.17,325.17]	#degR
H=[0,36089.239,65616.798,104986.878,154199.475,170603.675,200131.234,259186.352]	#ft
P=[2116.22,472.675801650081,114.343050672041,18.1283133205764,2.31620845720195,1.23219156244977,0.38030066501701,0.0215739175227548]	#lbf/ft^2
RHO=[0.00237691267925741,0.000706115448911997,0.000170813471460564,2.56600341257735e-05,2.76975106424479e-06,1.47347009326248e-06,4.87168173794687e-07,3.86714900013768e-08]	#slug/ft^3

q_list = [800,1200,1600]
M_list = np.linspace(0.4,10,100)
M_plot_list = []
altitude_list = []
for y in q_list:
    M_plot_list = []
    altitude_list = []
    T_list = []
    T0_list = []
    P0_list = []
    P_list = []
    for x in M_list:
        M = x
        p_static = (y*2/gamma)*(1/M**2)
        #print(p_static)
        if(p_static <= 2116.22):
            if(p_static > 472.675801650081):
                #print('test1')
                index = 0
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                PonPi = p_static / P_i
                TonTi = PonPi**(-K_i*R/g0)
                temp = TonTi * T_i
                h = H_i + (T_i*(TonTi - 1))/K_i 
                altitude_list.append(h)
                T_list.append(temp)
            
            elif(p_static > 114.343050672041):
                #print('test2')
                index = 1
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                temp = T_i
                PonPi = p_static / P_i
                h = H_i + (T_i*R*np.log(PonPi))/(-1*g0) 
                altitude_list.append(h)
                T_list.append(temp)
            elif(p_static > 18.1283133205764):
                #print('test3')
                index = 2
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                PonPi = p_static / P_i
                TonTi = PonPi**(-K_i*R/g0)
                temp = TonTi * T_i
                h = H_i + (T_i*(TonTi - 1))/K_i 
                altitude_list.append(h)
                T_list.append(temp)
            elif(p_static > 2.31620845720195):
                #print('test4')
                index = 3
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                PonPi = p_static / P_i
                TonTi = PonPi**(-K_i*R/g0)
                temp = TonTi * T_i
                h = H_i + (T_i*(TonTi - 1))/K_i 
                altitude_list.append(h)
                T_list.append(temp)
            elif(p_static > 1.23219156244977):
                #print('test5')
                index = 4
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                temp = T_i
                PonPi = p_static / P_i
                h = H_i + (T_i*R*np.log(PonPi))/(-1*g0) 
                altitude_list.append(h)
                T_list.append(temp)
            elif(p_static > 0.38030066501701):
                #print('test6')
                index = 5
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                PonPi = p_static / P_i
                TonTi = PonPi**(-K_i*R/g0)
                temp = TonTi * T_i
                h = H_i + (T_i*(TonTi - 1))/K_i 
                altitude_list.append(h)
                T_list.append(temp)
            elif(p_static > 0.0215739175227548):
                #print('test7')
                index = 6
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                PonPi = p_static / P_i
                TonTi = PonPi**(-K_i*R/g0)
                temp = TonTi * T_i
                h = H_i + (T_i*(TonTi - 1))/K_i 
                altitude_list.append(h)
                T_list.append(temp)
            elif(p_static <= 0.0215739175227548):
                #print('test8')
                index = 7
                K_i = K[index]
                T_i = T[index]
                H_i = H[index]
                RHO_i = RHO[index]
                P_i = P[index]
                
                temp = T_i
                PonPi = p_static / P_i
                h = H_i + (T_i*R*np.log(PonPi))/(-1*g0) 
                altitude_list.append(h)
                T_list.append(temp)
            #print(index)
            T0 = temp*(1 + ((gamma-1)/2)*M**2)
            print(p_static)
            P0 = p_static*(1 + ((gamma-1)/2)*M**2)**(gamma/(gamma-1))
            T0_list.append(T0)
            P0_list.append(P0)
            M_plot_list.append(M)
            P_list.append(p_static)

    plt.plot(altitude_list,P0_list,label='q Trajectory:'+str(y))
plt.grid()
plt.legend()
plt.xlabel('Altitude (ft)')
plt.ylabel('Recovery Temperature (R)')  
plt.title('Altitude vs. Recovery Temperature')


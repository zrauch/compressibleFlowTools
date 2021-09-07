# Generated with SMOP  0.41
#from libsmop import *
# atmosphere.m

from numpy import *
import matplotlib
import matplotlib.pyplot as plt
# setunits defines unit system-dependent things like R, g, and Cp
from setunits import *

# define global constants for both atmosphere functions
global RE,D,K,T,H,P,RHO
RE=20926476

D=array(
[[1,    -0.00356616,    518.67, 0,              2116.22,                0.00237691267925741],
[2,     0,              389.97, 36089.239,      472.675801650081,       0.000706115448911997],
[3,     0.00054864,     389.97, 65616.798,      114.343050672041,       0.000170813471460564],
[4,     0.00153619,     411.57, 104986.878,     18.1283133205764,       2.56600341257735e-05],
[5,     0,              487.17, 154199.475,     2.31620845720195,       2.76975106424479e-06],
[6,     -0.00109728,    487.17, 170603.675,     1.23219156244977,       1.47347009326248e-06],
[7,     -0.00219456,    454.17, 200131.234,     0.38030066501701,       4.87168173794687e-07],
[8,     0,              325.17, 259186.352,     0.0215739175227548,     3.86714900013768e-08]])

K=D[:,1]
T=D[:,2]
H=D[:,3]
P=D[:,4]
RHO=D[:,5]


def atmosphere4(q,Hvector):
    # default set GeometricFlag to 1 because who the fuck is going to give geopotential altitude as an input???
    GeometricFlag = 1 
    #varargin = atmosphere4.varargin
    #nargin = atmosphere4.nargin
    # function [temp,press,rho,Hgeopvector] = atmosphere4(Hvector,GeometricFlag)
    # Standard Atmospheric data based on the 1976 NASA Standard Atmoshere.
    # Hvector is a vector of altitudes.
    # If Hvector is Geometric altitude set GeometricFlag=1.
    # If Hvector is Geopotential altitude set GeometricFlag=0.
    # Temp, press, and rho are temperature, pressure and density
    # output vectors the same size as Hgeomvector.
    # Output vector Hgeopvector is a vector of corresponding geopotential altitudes (ft).
    # This atmospheric model is good for altitudes up to 295,000 geopotential ft.
    # Ref: Intoduction to Flight Test Engineering by Donald T. Ward and Thomas W. Strganac
    # index   Lapse rate   Base Temp     Base Geopo Alt        Base Pressure            Base Density
    #   i     Ki(degR/ft)  Ti(degR)        Hi(ft)              P, lbf/ft^2           RHO, slug/ft^3

    M = zeros(size(Hvector))
    T0 = zeros(size(Hvector))
    temp = zeros(size(Hvector))
    p0 = zeros(size(Hvector))
    press = zeros(size(Hvector))
    rho = zeros(size(Hvector))
    Hgeopvector = zeros(size(Hvector))

    # Convert from geometric altitude to geopotental altitude, if necessary.
    if GeometricFlag:
        Hgeopvector=(dot(RE,Hvector)) / (RE + Hvector)
        #disp('Convert from geometric altitude to geopotential altitude in feet')
    else:
        Hgeopvector=copy(Hvector)

    ih=len(Hgeopvector)
    n1=argwhere(Hgeopvector <= H[1])
    n2=argwhere((Hgeopvector<=H[2]) & (Hgeopvector>H[1]))
    n3=argwhere((Hgeopvector<=H[3]) & (Hgeopvector>H[2]))
    n4=argwhere((Hgeopvector<=H[4]) & (Hgeopvector>H[3]))
    n5=argwhere((Hgeopvector<=H[5]) & (Hgeopvector>H[4]))
    n6=argwhere((Hgeopvector<=H[6]) & (Hgeopvector>H[5]))
    n7=argwhere((Hgeopvector<=H[7]) & (Hgeopvector>H[6]))
    n8=argwhere((Hgeopvector<=295000) & (Hgeopvector>H[7]))
    icorrect=len(n1) + len(n2) + len(n3) + len(n4) + len(n5) + len(n6) + len(n7) + len(n8)
    if icorrect < ih:
        disp('One or more altitutes is above the maximum for this atmospheric model')
        icorrect
        ih

    # Index 1, Troposphere, K1= -.00356616
    if len(n1) > 0:
        i=0
        h=Hgeopvector[n1]
        TonTi=1 + dot(K[i],(h - H[i])) / T[i]
        temp[n1]=dot(TonTi,T[i])
        PonPi=TonTi ** (- g0 / (dot(K[i],R)))
        press[n1]=dot(P[i],PonPi)
        RonRi=TonTi ** (- g0 / (dot(K[i],R)) - 1)
        rho[n1]=dot(RHO[i],RonRi)
        M[n1] = sqrt(2.0*q/(gamma*press[n1]))
        p0[n1] = press[n1]*(1 + (gamma-1)/2 * M[n1]**2)**(gamma/(gamma-1))
        T0[n1] = temp[n1]*(1 + (gamma-1)/2 * M[n1]**2)

    # Index 2,  K2= 0
    if len(n2) > 0:
        i=1
        h=Hgeopvector[n2]
        temp[n2]=T[i]
        PonPi=exp(dot(- g0,(h - H[i])) / (dot(T[i],R)))
        press[n2]=dot(P[i],PonPi)
        RonRi=copy(PonPi)
        rho[n2]=dot(RHO[i],RonRi)
        M[n2] = sqrt(2.0*q/(gamma*press[n2]))
        p0[n2] = press[n2]*(1 + (gamma-1)/2 * M[n2]**2)**(gamma/(gamma-1))
        T0[n2] = temp[n2]*(1 + (gamma-1)/2 * M[n2]**2)

    # Index 3,  K3= .00054864
    if len(n3) > 0:
        i=2
        h=Hgeopvector[n3]
        TonTi=1 + dot(K[i],(h - H[i])) / T[i]
        temp[n3]=dot(TonTi,T[i])
        PonPi=TonTi ** (- g0 / (dot(K[i],R)))
        press[n3]=dot(P[i],PonPi)
        RonRi=TonTi ** (- g0 / (dot(K[i],R)) - 1)
        rho[n3]=dot(RHO[i],RonRi)
        M[n3] = sqrt(2.0*q/(gamma*press[n3]))
        p0[n3] = press[n3]*(1 + (gamma-1)/2 * M[n3]**2)**(gamma/(gamma-1))
        T0[n3] = temp[n3]*(1 + (gamma-1)/2 * M[n3]**2)

    # Index 4,  K4= .00153619
    if len(n4) > 0:
        i=3
        h=Hgeopvector[n4]
        TonTi=1 + dot(K[i],(h - H[i])) / T[i]
        temp[n4]=dot(TonTi,T[i])
        PonPi=TonTi ** (- g0 / (dot(K[i],R)))
        press[n4]=dot(P[i],PonPi)
        RonRi=TonTi ** (- g0 / (dot(K[i],R)) - 1)
        rho[n4]=dot(RHO[i],RonRi)
        M[n4] = sqrt(2.0*q/(gamma*press[n4]))
        p0[n4] = press[n4]*(1 + (gamma-1)/2 * M[n4]**2)**(gamma/(gamma-1))
        T0[n4] = temp[n4]*(1 + (gamma-1)/2 * M[n4]**2)

    # Index 5,  K5= 0
    if len(n5) > 0:
        i=4
        h=Hgeopvector[n5]
        temp[n5]=T[i]
        PonPi=exp(dot(- g0,(h - H[i])) / (dot(T[i],R)))
        press[n5]=dot(P[i],PonPi)
        RonRi=copy(PonPi)
        rho[n5]=dot(RHO[i],RonRi)
        M[n5] = sqrt(2.0*q/(gamma*press[n5]))
        p0[n5] = press[n5]*(1 + (gamma-1)/2 * M[n5]**2)**(gamma/(gamma-1))
        T0[n5] = temp[n5]*(1 + (gamma-1)/2 * M[n5]**2)

    # Index 6,  K6= -.00109728
    if len(n6) > 0:
        i=5
        h=Hgeopvector[n6]
        TonTi=1 + dot(K[i],(h - H[i])) / T[i]
        temp[n6]=dot(TonTi,T[i])
        PonPi=TonTi ** (- g0 / (dot(K[i],R)))
        press[n6]=dot(P[i],PonPi)
        RonRi=TonTi ** (- g0 / (dot(K[i],R)) - 1)
        rho[n6]=dot(RHO[i],RonRi)
        M[n6] = sqrt(2.0*q/(gamma*press[n6]))
        p0[n6] = press[n6]*(1 + (gamma-1)/2 * M[n6]**2)**(gamma/(gamma-1))
        T0[n6] = temp[n6]*(1 + (gamma-1)/2 * M[n6]**2)

    # Index 7,  K7= -.00219456
    if len(n7) > 0:
        i=6
        h=Hgeopvector[n7]
        TonTi=1 + dot(K[i],(h - H[i])) / T[i]
        temp[n7]=dot(TonTi,T[i])
        PonPi=TonTi ** (- g0 / (dot(K[i],R)))
        press[n7]=dot(P[i],PonPi)
        RonRi=TonTi ** (- g0 / (dot(K[i],R)) - 1)
        rho[n7]=dot(RHO[i],RonRi)
        M[n7] = sqrt(2.0*q/(gamma*press[n7]))
        p0[n7] = press[n7]*(1 + (gamma-1)/2 * M[n7]**2)**(gamma/(gamma-1))
        T0[n7] = temp[n7]*(1 + (gamma-1)/2 * M[n7]**2)

    # Index 8,  K8= 0
    if len(n8) > 0:
        i=7
        h=Hgeopvector[n8]
        temp[n8]=T[i]
        PonPi=exp(dot(- g0,(h - H[i])) / (dot(T[i],R)))
        press[n8]=dot(P[i],PonPi)
        RonRi=copy(PonPi)
        rho[n8]=dot(RHO[i],RonRi)
        M[n8] = sqrt(2.0*q/(gamma*press[n8]))
        p0[n8] = press[n8]*(1 + (gamma-1)/2 * M[n8]**2)**(gamma/(gamma-1))
        T0[n8] = temp[n8]*(1 + (gamma-1)/2 * M[n8]**2)

    return [M,T0,temp,p0,press,rho,Hgeopvector]

# @author: geigerr
# revised to be formatted as a function by zrauch
def atmosphereMach(q,M_list):
    # initialize necessary lists for calculation and output of properties
    altitude_list = []
    T0_list = []; T_list = []
    P0_list = []; P_list = []
    for M in M_list:
        # atmosphereMach operates over a range of Mach numbers, not a range of altitudes like atmosphere4
        p_static = 2.0*q/(gamma*M**2)
        # next determine which range of h you are in based on Mach number and static p
        if(p_static > P[1]):
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
        
        elif(p_static > P[2]):
            index = 1
            K_i = K[index]
            T_i = T[index]
            H_i = H[index]
            RHO_i = RHO[index]
            P_i = P[index]
            
            temp = T_i
            PonPi = p_static / P_i
            h = H_i + (T_i*R*log(PonPi))/(-1*g0) 
            altitude_list.append(h)
            T_list.append(temp)
        elif(p_static > P[3]):
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
        elif(p_static > P[4]):
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
        elif(p_static > P[5]):
            index = 4
            K_i = K[index]
            T_i = T[index]
            H_i = H[index]
            RHO_i = RHO[index]
            P_i = P[index]
            
            temp = T_i
            PonPi = p_static / P_i
            h = H_i + (T_i*R*log(PonPi))/(-1*g0) 
            altitude_list.append(h)
            T_list.append(temp)
        elif(p_static > P[6]):
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
        elif(p_static > P[7]):
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
        elif(p_static <= P[7]):
            index = 7
            K_i = K[index]
            T_i = T[index]
            H_i = H[index]
            RHO_i = RHO[index]
            P_i = P[index]
            
            temp = T_i
            PonPi = p_static / P_i
            h = H_i + (T_i*R*log(PonPi))/(-1*g0) 
            altitude_list.append(h)
            T_list.append(temp)
        elif(p_static>P[0]):
            print("Error! Static pressure value is greater than static pressure at sea level!")
            print("Exiting now, please check your inputs to make sure the static pressure is physically possible!")
            exit(0)

        T0_list.append(temp*(1 + ((gamma-1)/2)*M**2))
        P0_list.append(p_static* (1 + ((gamma-1)/2)*M**2)**(gamma/(gamma-1)))
        P_list.append(p_static)

    return [T0_list,T_list,P0_list,P_list,altitude_list]

test = 0
if test:
    # input values necessary for test case of each function
    qlist = [800, 1200, 1600] # given values of dynamic pressure [lbf/ft^2]
    altrange = linspace(0,100000,101)
    M_list = linspace(0.4,10,1001)
    #atmosphere4 test case
    # plot output results for flight corridor over a range of altitudes
    fig,ax1 = plt.subplots()
    plt.title(r"Mach Number with Altitude for Fixed Dynamic Pressures")
    ax1.set_xlabel(r"Mach Number")
    ax1.set_ylabel(r"Altitude [ft]")
    for i,q in enumerate(qlist):
        [Mach,T0,Ts,P0,Ps,rho,hgeop] = atmosphere4(q,altrange)
        ax1.plot(Mach,altrange,label="q="+str(q)+"psf")
    ax1.grid()
    plt.legend()
    fig.tight_layout()
    plt.savefig('altitude_Mach.pdf')
    plt.close()

    #atmosphereMach test case
    fig,ax1 = plt.subplots()
    for i,q in enumerate(qlist):
        [T0,Ts,P0,Ps,h] = atmosphereMach(q,M_list)
        plt.plot(M_list,h,label='q Trajectory:'+str(q))

    fig,ax1 = plt.subplots()
    plt.title(r"Mach Number with Altitude for Fixed Dynamic Pressures")
    ax1.set_xlabel(r"Mach Number")
    ax1.set_ylabel(r"Altitude [ft]")
    for i,q in enumerate(qlist):
        [T0,Ts,P0,Ps,h] = atmosphereMach(q,M_list)
        plt.plot(M_list,h,label='q Trajectory:'+str(q))
    ax1.grid()
    plt.legend()
    fig.tight_layout()
    plt.savefig('altitude_Mach_2.pdf')
    plt.close()
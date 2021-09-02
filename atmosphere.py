# Generated with SMOP  0.41
#from libsmop import *
# atmosphere.m

from numpy import *

def atmosphere4(Hvector=None, GeometricFlag=None, *args, **kwargs):
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

    # define constants
    D=array(
    [[1,    -0.00356616,    518.67, 0,              2116.22,                0.00237691267925741],
    [2,     0,              389.97, 36089.239,      472.675801650081,       0.000706115448911997],
    [3,     0.00054864,     389.97, 65616.798,      114.343050672041,       0.000170813471460564],
    [4,     0.00153619,     411.57, 104986.878,     18.1283133205764,       2.56600341257735e-05],
    [5,     0,              487.17, 154199.475,     2.31620845720195,       2.76975106424479e-06],
    [6,     -0.00109728,    487.17, 170603.675,     1.23219156244977,       1.47347009326248e-06],
    [7,     -0.00219456,    454.17, 200131.234,     0.38030066501701,       4.87168173794687e-07],
    [8,     0,              325.17, 259186.352,     0.0215739175227548,     3.86714900013768e-08]])
    
    R=1716.55
    gamma=1.4
    g0=32.17405
    RE=20926476

    K=D[:,1]
    T=D[:,2]
    H=D[:,3]
    P=D[:,4]
    RHO=D[:,5]

    temp=zeros(size(Hvector))
    press=zeros(size(Hvector))
    rho=zeros(size(Hvector))
    Hgeopvector=zeros(size(Hvector))

    # Convert from geometric altitude to geopotental altitude, if necessary.
    if GeometricFlag:
        Hgeopvector=(dot(RE,Hvector)) / (RE + Hvector)
        disp('Convert from geometric altitude to geopotential altitude in feet')
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

    # Index 2,  K2= 0
    if len(n2) > 0:
        i=1
        h=Hgeopvector[n2]
        temp[n2]=T[i]
        PonPi=exp(dot(- g0,(h - H[i])) / (dot(T[i],R)))
        press[n2]=dot(P[i],PonPi)
        RonRi=copy(PonPi)
        rho[n2]=dot(RHO[i],RonRi)

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

    # Index 5,  K5= 0
    if len(n5) > 0:
        i=4
        h=Hgeopvector[n5]
        temp[n5]=T[i]
        PonPi=exp(dot(- g0,(h - H[i])) / (dot(T[i],R)))
        press[n5]=dot(P[i],PonPi)
        RonRi=copy(PonPi)
        rho[n5]=dot(RHO[i],RonRi)

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

    # Index 8,  K8= 0
    if len(n8) > 0:
        i=7
        h=Hgeopvector[n8]
        temp[n8]=T[i]
        PonPi=exp(dot(- g0,(h - H[i])) / (dot(T[i],R)))
        press[n8]=dot(P[i],PonPi)
        RonRi=copy(PonPi)
        rho[n8]=dot(RHO[i],RonRi)

    return [temp,press,rho,Hgeopvector]
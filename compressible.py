# Import libraries
import numpy as np
import scipy
import scipy.sparse as scysparse
import scipy.sparse.linalg
import re
import os
import math
from math import pi as PI

## variables
global gamma, R_SI,R_ENG

gamma = 1.4
R_SI = 287 # J/kg-K


def D(M,gamma): # mass flow function, see Mattingly and Boyer
	M_2 = pow(M,2)
	isen = (1 + ((gamma-1)/2)*M_2)
	power = (gamma+1)/(2*(gamma-1))
	D = M/pow(isen,power)

	if M == 0:
		dD_dM = float("inf")
	else:
		dD_dM = (D/M)*(1-M_2)/isen
	return D,dD_dM

def G(M,gamma): # thrust/momentum flow function, see Mattingly and Boyer
	M_2 = pow(M,2)
	isen = (1 + ((gamma-1)/2)*M_2)
	G = (1+gamma*M_2)/pow(isen,(gamma/(gamma-1)))
	dG_dM = gamma*G*M*(2/(1+gamma*M_2) - 1/(isen))
	return G, dG_dM

def N(M,gamma): # the no-named function?? see Mattingly and Boyer
	M_2 = pow(M,2)
	isen = (1 + ((gamma-1)/2)*M_2)
	N = D(M,gamma)[0]/G(M,gamma)[0]
	if M == 0:
		dN_dM = float("inf")
	else:
		dN_dM = N*(1/M + ((gamma-1)/2)*(M/isen) - (2*gamma*M/(1+gamma*M_2)))
	return N, dN_dM

def areaMachRelation(M_value,A_ratio,M_init): # derived from mass flow using isentropic relations, see Mattingly&Boyer or Zucrow (1976)
	DM = D(M_value,gamma)[0]
	M_guess = M_init
	eps = 1 # initialize with eps > delta
	while (eps > 1e-6):
		D_guess,dD_dM_guess = D(M_guess,gamma)
		M_new = M_guess + (A_ratio*DM-D_guess)/dD_dM_guess
		eps = np.abs(M_new-M_guess)
		M_guess = M_new
	return M_new

def newtonRaphsonMach(M1,A_ratio): # newton-raphson solver for M2 given M1 and area ratio of a nozzle
	eps = 1
	D1,unused = D(M1,gamma)
	M_guess = 0.5
	while eps > 1e-4:
		D_guess,dD_dM_guess = D(M_guess,gamma)
		M2 = M_guess + (D1*A_ratio-D_guess)/dD_dM_guess
		eps = np.absolute(M2-M_guess)
		M_guess = M2
	return M_guess

def newtonRaphsonMachSHOCK(N1,gamma): # newton-raphson solver for M2 after a shock, given N(M1)
	eps = 1
	M_guess = 0.75
	while eps > 1e-4:
		N_guess,dN_dM_guess = N(M_guess,gamma)
		M2 = M_guess + (N1-N_guess)/dN_dM_guess
		eps = np.absolute(M2-M_guess)
		M_guess = M2
	return M_guess

def newtonRaphsonSupersonic(M1,gamma): # newton-raphson solver for a supersonic M2 given a subsonic M1
	eps = 1
	D1,unused = D(M1,gamma)
	M_guess = 2.5
	while eps > 1e-4:
		D_guess,dD_dM_guess = D(M_guess,gamma)
		M2 = M_guess + (D1-D_guess)/dD_dM_guess
		eps = np.absolute(M2-M_guess)
		M_guess = M2
	return M_guess

def isentropicRelationCalculator(M,gamma,printout): # isentropic relations calculator, does not currently include sonic relations
	isen = 1 + (gamma-1)/2*M**2
	T0_T = isen
	p0_p = isen**(gamma/(gamma-1))
	rho0_rho = isen**(1/(gamma-1))

	if printout:
		print("Normal Shock Calculator Results")
		print("Give inputs M =",round(M,3),"gamma =",gamma)
		print("==================================")
		print("T0/T =",T0_T)
		print("p0/p =",p0_p)
		print("rho0/rho =",rho0_rho)
		print("\n")
	return [T0_T, p0_p, rho0_rho]

def normalShockCalculator(M,gamma,printout): # normal shock calculator, equations derived in Mattingly and Boyer
	M2 = np.sqrt(((gamma-1)*M**2 + 2)/(2*gamma*M**2 - (gamma-1)))
	p02_p01 = ( ((gamma+1)*M**2)/((gamma-1)*M**2 + 2) )**(gamma/(gamma-1)) * ( (gamma+1)/(2*gamma*M**2-(gamma-1)) )**(1/(gamma-1))
	T2_T1 = (2*gamma*M**2-(gamma-1))*((gamma-1)*M**2+2)/((gamma+1)**2*M**2)
	T02_T01 = 1
	p2_p1 = (2*gamma*M**2-(gamma-1))/(gamma+1)
	rho2_rho1 = ((gamma+1)*M**2)/((gamma-1)*M**2+2)

	if printout:
		print("Normal Shock Calculator Results")
		print("Give inputs M =",round(M,3),"gamma =",gamma)
		print("==================================")
		print("M2 =",M2)
		print("p02/p01 =",p02_p01)
		print("p2/p1 =",p2_p1)
		print("T02/T01 =",T02_T01)
		print("T2/T1 =",T2_T1)
		print("rho2/rho1 =",rho2_rho1,"\n")
	return [M2,p02_p01,p2_p1,T02_T01,T2_T1,rho2_rho1]


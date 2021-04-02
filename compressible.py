# Import libraries
import numpy as np
import scipy
import scipy.sparse as scysparse
import scipy.sparse.linalg
import re
import os
import math
from math import pi as PI
from scipy.optimize import fsolve

## variables
global gamma, R_SI,R_ENG

gamma = 1.4
R_SI = 287 # J/kg-K
R_ENG = 53.353 # ft-lbf/lbm-R -- needs multiplied by g for most applications

def massFlow(p0,T0,A,M):
	mdot = p0*A/np.sqrt(T0)*np.sqrt(gamma/R_SI)*D(M,gamma)[0]
	return mdot

def solveMfromMassFlow(mdot,p0,T0,A,regime):
	LHS = mdot*np.sqrt(T0)/(p0*A)*np.sqrt(R_SI/gamma)
	def equation(x):
		RHS = x / (1+(gamma-1)/2 * x**2)**((gamma+1)/(2*(gamma-1)))
		return RHS-LHS

	if regime == 'subsonic':
		M = fsolve(equation,0.5)
	elif regime == 'supersonic':
		M = fsolve(equation,2.5)
	return M[0]



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

def shockAreaIterator(As_At_guess,At,p0,pb):
	eps = 1 # initialize with eps > delta
	while eps > 5e-8:
		Mshock = areaMachRelation(1,1/As_At_guess,2.0)
		[Tratio,Pratio,rhoRatio] = isentropicRelationCalculator(Mshock,gamma,0)
		p1 = 1/Pratio*p0
		[M2,stagPratio,Pratio,stagTratio,Tratio,rhoRatio] = normalShockCalculator(Mshock,gamma,0)
		pb_guess = Pratio*p1
		if pb_guess>pb:
			As_At_guess+=1e-6
		elif pb_guess<pb:
			As_At_guess-=1e-6
		eps = np.abs(pb_guess-pb)/pb
	As = As_At_guess*At
	print("Shock Area =",As,"[m2]")
	return As


def lstarFromMach(M,f,D):
	term1 = -1/(gamma) - (gamma+1)/(2*gamma)*np.log(1.0/(1+(gamma-1)/2))
	term2 = -1/(gamma*M**2) - (gamma+1)/(2*gamma)*np.log(M**2/(1+(gamma-1)*M**2/2))
	flD = term1-term2
	lstar = flD*D/(4*f)
	return lstar

def MachFromLstar(f,Lstar,D,regime):
	def equation(x):
		term1 = -1/(gamma) - (gamma+1)/(2*gamma)*np.log(1/(1+(gamma-1)/2))
		term2 = -1/(gamma*x**2) - (gamma+1)/(2*gamma)*np.log(x**2/(1+(gamma-1)*x**2/2))
		return term1-term2-4*f*Lstar/D
	if regime == 'subsonic':
		M = fsolve(equation,0.5)
	elif regime == 'supersonic':
		M = fsolve(equation,2.5)

	return M[0]

def fannoFlowProperties(M):
	Trat = ((gamma+1)/2)/(1+(gamma-1)/2*M**2)
	prat = (1/M)*np.sqrt(Trat)
	vrat = M*np.sqrt(Trat)
	stagPrat = (1/M)*(Trat)**((gamma+1)/(2*(1-gamma)))
	rhorat = 1/vrat
	return [Trat,prat,vrat,stagPrat,rhorat]

def rayleighStarProperties(M):
	term = (1+gamma)/(1+gamma*M**2)
	pratio = term
	Tratio = M**2 * term**2
	vratio = term * M**2
	p0ratio = term * (1 + (gamma-1)/2 * M**2)**(gamma/(gamma-1))
	T0ratio = term**2 * M**2 *(1+(gamma-1)/2 * M**2)/(1+(gamma-1)/2)
	return [pratio,Tratio,vratio,p0ratio,T0ratio]

def rayleighPropertyChanges(M1,M2):
	term1 = 1 + gamma*M1**2
	term2 = 1 + gamma*M2**2
	term3 = 1 + (gamma-1)/2 * M1**2
	term4 = 1 + (gamma-1)/2 * M2**2

	pratio = term1/term2
	Tratio = (term1/term2)**2 * (M2/M1)**2
	rhoRatio = (term1/term2) * (M1/M2)**2
	p0ratio = (term1/term2) * (term4/term3)**(gamma/(gamma-1))
	T0ratio = (term1/term2)**2 * (M2/M1)**2 * (term4/term3)
	deltaS_Cp = np.log((M2/M1)**2 * (term1/term2)**((gamma-1)/gamma))
	return [pratio,Tratio,rhoRatio,p0ratio,T0ratio,deltaS_Cp]


def rayleighFlowSolver(value, ratio, regime):
	def equation(x):
		term1 = value
		simp = (1+gamma)/(1+gamma*x**2)
		if ratio == 'p':
			term2 = simp * (1 + (gamma-1)/2 * x**2)**(gamma/(gamma-1))
		elif ratio == 'T':
			term2 = x**2 * simp**2
		elif ratio == 'v':
			term2 = simp * x**2
		elif ratio == 'p0':
			term2 = simp * (1 + (gamma-1)/2 * x**2)**(gamma/(gamma-1))
		elif ratio == 'T0':
			term2 = simp**2 * x**2 *(1+(gamma-1)/2 * x**2)/(1+(gamma-1)/2)
		return term2-term1

	if regime == 'subsonic':
		M = fsolve(equation,0.5)
	elif regime == 'supersonic':
		M = fsolve(equation,2.5)

	return M[0]

def obliqueShockCalculator(theta, M, gamma, shockType):
    #betaRange = np.linspace(0.001,np.pi/2,1001) # brute force through all possible beta values
    #for b in betaRange: # loop through range of beta values until the RHS is larger than the LHS (requires computing/precision)
    def equation(x):
    	LHS = np.tan(theta) # LHS is constant
    	RHS = (2/np.tan(x))*(M**2*np.sin(x)**2 - 1)/(M**2*(gamma+np.cos(2*x)) + 2)
    	return RHS-LHS

    if shockType == 'Weak':
    	b = fsolve(equation,0.2)[0]
    elif shockType == 'Strong':
    	b = fsolve(equation,np.pi/2-0.2)[0]

    term = M**2*np.sin(b)**2
    M2 = np.sqrt((1/np.sin(b-theta)**2) * (1 + (gamma-1)/2*term)/(gamma*term - (gamma-1)/2)) # from NASA website
    T_ratio = (2*gamma*term - (gamma-1))*((gamma-1)*term + 2)/((gamma+1)**2*term) # from NASA website
    p_ratio = (2*gamma*term - (gamma-1))/(gamma+1) # from NASA website
    term1 = pow(((gamma+1)*term)/((gamma-1)*term + 2),(gamma/(gamma-1)))
    term2 = pow(((gamma+1)/(2*gamma*term - (gamma-1))),(1/(gamma-1)))
    stagp_ratio = term1*term2
    return b, M2, T_ratio, p_ratio, stagp_ratio

def obliqueShockCalculator2(beta, M, gamma):
    #betaRange = np.linspace(0.001,np.pi/2,1001) # brute force through all possible beta values
    #for b in betaRange: # loop through range of beta values until the RHS is larger than the LHS (requires computing/precision)
    def equation(x):
    	LHS = np.tan(x) # LHS is constant
    	RHS = (2/np.tan(beta))*(M**2*np.sin(beta)**2 - 1)/(M**2*(gamma+np.cos(2*beta)) + 2)
    	return RHS-LHS

    theta = fsolve(equation,0.2)[0]

    term = M**2*np.sin(beta)**2
    M2 = np.sqrt((1/np.sin(beta-theta)**2) * (1 + (gamma-1)/2*term)/(gamma*term - (gamma-1)/2)) # from NASA website
    T_ratio = (2*gamma*term - (gamma-1))*((gamma-1)*term + 2)/((gamma+1)**2*term) # from NASA website
    p_ratio = (2*gamma*term - (gamma-1))/(gamma+1) # from NASA website
    term1 = pow(((gamma+1)*term)/((gamma-1)*term + 2),(gamma/(gamma-1)))
    term2 = pow(((gamma+1)/(2*gamma*term - (gamma-1))),(1/(gamma-1)))
    stagp_ratio = term1*term2
    return theta, M2, T_ratio, p_ratio, stagp_ratio

def prandtlMeyerExpansion(M,gamma):
	nu = np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)*(M**2-1)/(gamma+1))) - np.arctan(np.sqrt(M**2-1))
	return nu

def solvePMexpansion(nu,gamma):
	M2guess = 1.5
	def eqn(x,*data):
		nu,gamma = data
		LHS = nu
		RHS = np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)*(x**2-1)/(gamma+1))) - np.arctan(np.sqrt(x**2-1))
		return RHS-LHS
	return fsolve(eqn,M2guess,args=(nu,gamma))[0]

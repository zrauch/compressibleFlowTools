import numpy as np
import scipy as sp
import sympy as sym
import math
from scipy.optimize import fsolve
from setunits import *

# generic RK4 stencil for numerical integration
def RK4(y,x,dx,function,parameters):
	k1 = function(x, y, parameters)
	k2 = function(x+0.5*dx, y+0.5*k1*dx, parameters)
	k3 = function(x+0.5*dx, y+0.5*k2*dx, parameters)
	k4 = function(x+dx, y+k3*dx, parameters)
	y_new = y + dx/6*(k1+2*k2+2*k3+k4)
	return y_new

def RK4_stagP(M,P0,x,dx,function,parameters):
	k1 = function(x, M, parameters)
	k2 = function(x+0.5*dx, M, parameters)
	k3 = function(x+0.5*dx, M, parameters)
	k4 = function(x+dx, M, parameters)
	P0_new = P0 + P0*dx/6*(k1+2*k2+2*k3+k4)
	return P0_new

# massFlow :: equation for mass flow rate in terms of Mach, Area, and stagnation conditions
# inputs :: stagnation pressure, p0			[Pa or psf]
# 			stagnation temperature,T0 		[K or R]
# 			area at point of evaluation, A 	[m2 or ft2]
# 			Mach number, M
# output :: mass flow rate, mdot 			[kg/s or lbm/s]
def massFlow(p0,T0,A,M):
	mdot = p0*A/np.sqrt(T0)*np.sqrt(gamma/R)*D(M,gamma)[0]
	return mdot

def solveAfromMassFlow(mdot,M,p0,T0,R,gamma):
	A = mdot*np.sqrt(T0*R/gamma)/(p0*M) * (1 + ((gamma-1)/2) * M**2)**((gamma+1)/(2*(gamma-1)))
	return A

# solveMfromMassFlow :: uses fsolve to determine the solutions of Mach number for which a 
# 						certain mass flow rate value can be attained at given conditions.
# inputs :: required mass flow rate, mdot [kg/s or lbm/s]
#			stagnation pressure, p0			[Pa or psf]
# 			stagnation temperature,T0 		[K or R]
# 			area at point of evaluation, A 	[m2 or ft2]
#			flow regime, regime				['subsonic' or 'supersonic']
# outputs :: Mach number, M	
def solveMfromMassFlow(mdot,p0,T0,A,regime):
	LHS = mdot*np.sqrt(T0)/(p0*A)*np.sqrt(R/gamma)
	def equation(x):
		RHS = x / (1+(gamma-1)/2 * x**2)**((gamma+1)/(2*(gamma-1)))
		return RHS-LHS
	if regime == 'subsonic':
		M = fsolve(equation,0.5)
	elif regime == 'supersonic':
		M = fsolve(equation,2.5)
	return M[0]

# D (mass flow function) :: Returns the mass flow and derivative for a given M and gamma
def D(M,gamma): # see Mattingly and Boyer
	M_2 = pow(M,2)
	isen = (1 + ((gamma-1)/2)*M_2)
	power = (gamma+1)/(2*(gamma-1))
	D = M/pow(isen,power)

	if M == 0:
		dD_dM = float("inf")
	else:
		dD_dM = (D/M)*(1-M_2)/isen
	return D,dD_dM

def Dfunction(M,gamma): # see Mattingly and Boyer
	M_2 = pow(M,2)
	isen = (1 + ((gamma-1)/2)*M_2)
	power = (gamma+1)/(2*(gamma-1))
	D = M/pow(isen,power)

	if M == 0:
		dD_dM = float("inf")
	else:
		dD_dM = (D/M)*(1-M_2)/isen
	return D,dD_dM

# G (momentum flow function) :: Returns the thrust/momentum flow and derivative for a given M and gamma
def G(M,gamma): # see Mattingly and Boyer
	M_2 = pow(M,2)
	isen = (1 + ((gamma-1)/2)*M_2)
	G = (1+gamma*M_2)/pow(isen,(gamma/(gamma-1)))
	dG_dM = gamma*G*M*(2/(1+gamma*M_2) - 1/(isen))
	return G, dG_dM

# N (stream flow function) :: Returns the stream flow and derivative for a given M and gamma
def N(M,gamma): # see Mattingly and Boyer
	M_2 = pow(M,2)
	isen = (1 + ((gamma-1)/2)*M_2)
	N = D(M,gamma)[0]/G(M,gamma)[0]
	if M == 0:
		dN_dM = float("inf")
	else:
		dN_dM = N*(1/M + ((gamma-1)/2)*(M/isen) - (2*gamma*M/(1+gamma*M_2)))
	return N, dN_dM

# areaMachRelation :: uses previously defined D function to iteratively solve for the value of M
#					  at a downstream location for an isentropic flow
# inputs ::	known Mach number, M_value
#			the area ratio of geometry, A_ratio
# 			initial Mach number guess, M_init
# outputs :: final (downstream) Mach number, M_new
# NOTE :: the return value of M_new depends on M_init. if M_init > 1, M_new will be > 1. If M_init < 1, M_new < 1.
def areaMachRelation(M_value,A_ratio,M_init,*param): # see Mattingly&Boyer or Zucrow (1976)
	gamma = param[0]
	DM = D(M_value,gamma)[0]
	M_guess = M_init
	eps = 1 # initialize with eps > delta
	while eps > 1e-8:
		D_guess,dD_dM_guess = D(M_guess,gamma)
		M_new = M_guess + (A_ratio*DM-D_guess)/dD_dM_guess
		eps = np.abs(M_new-M_guess)
		M_guess = M_new
	return M_new


def areaMachRelation2(M1, Mx):
	# NOTE: A_ratio is h1/hx
	A_ratio = (Mx/M1) * ((2 + (gamma-1)*M1**2)/(2 + (gamma-1)*Mx**2))**((gamma+1)/(2*(gamma-1)))
	return A_ratio

def nonIsentropicAreaMachRelation(M1,M2, T01_T02, P02_P01,gamma):
	# NOTE: A_ratio is h1/h2
	A_ratio = np.sqrt(T01_T02) * (P02_P01) *(M2/M1) * ((2 + (gamma-1)*M1**2)/(2 + (gamma-1)*M2**2))**((gamma+1)/(2*(gamma-1)))
	return A_ratio


# newtonRaphsonMach :: uses previously defined D function to iteratively solve for the value of Mach number
# 					   downstream in a subsonic flow knowing the geometry.
# inputs :: initial mach number, M1
# 			area ratio of geometry, A_ratio
# outputs :: final (downstream) Mach number, M_guess
# NOTE :: this function only works for subsonic relationships, areaMachRelation will work for subsonic or
# 		  supersonic, depending on the given M_init
def areaMachSubsonic(M1,A_ratio): # newton-raphson solver for M2 given M1 and area ratio of a nozzle
	eps = 1
	D1,unused = D(M1,gamma)
	M_guess = 0.5
	while eps > 1e-8:
		D_guess,dD_dM_guess = D(M_guess,gamma)
		M2 = M_guess + (D1*A_ratio-D_guess)/dD_dM_guess
		eps = np.absolute(M2-M_guess)
		M_guess = M2
	return M_guess

# newtonRaphsonShock :: uses previously defined N function to iteratively solve for the value of Mach number
# downstream in a flow, since N is constant across a shock.
# inputs :: upstream SUPERSONIC mach number, M1
# 			ratio of specific heats of working fluid, gamma
# outputs :: the subsonic Mach number which will be immediately on the downstream side of the shock, M_guess
# NOTE :: This function only works if N1 was evaluated for a Mach number > 1. Otherwise, it will fail
def newtonRaphsonShock(M1,gamma):
	eps = 1
	N1, unused = N(M1,gamma)
	M_guess = 0.75
	while eps > 1e-4:
		N_guess,dN_dM_guess = N(M_guess,gamma)
		M2 = M_guess + (N1-N_guess)/dN_dM_guess
		eps = np.absolute(M2-M_guess)
		M_guess = M2
	return M_guess

# newtonRaphsonSupersonic :: uses previously defined D function to iteratively solve for the value of Mach number
# downstream in a flow which develops from subsonic to supersonic.
# inputs :: upstream SUBSONIC Mach number, M1
# 			ratio of specific heats of working fluid, gamma
# outputs :: the supersonic Mach number after the flow choke, M_guess
# NOTE :: this function requires M1 < 1 to be correct!
def newtonRaphsonSupersonic(M1,gamma):
	eps = 1
	D1,unused = D(M1,gamma)
	M_guess = 2.5
	while eps > 1e-4:
		D_guess,dD_dM_guess = D(M_guess,gamma)
		M2 = M_guess + (D1-D_guess)/dD_dM_guess
		eps = np.absolute(M2-M_guess)
		M_guess = M2
	return M_guess

# isentropicRelationCalculator returns all isentropic ratios given Mach number and gamma for a flow.
# NOTE: printout is an optional feature which will print out the results to terminal.
def isentropicRelationCalculator(M,gamma,printout):
	isen = 1 + (gamma-1)/2*M**2
	T0_T = isen
	p0_p = isen**(gamma/(gamma-1))
	rho0_rho = isen**(1/(gamma-1))

	if printout:
		print("Isentropic Flow Calculator Results")
		print("Give inputs M =",round(M,3),"gamma =",gamma)
		print("==================================")
		print("T0/T =",T0_T)
		print("p0/p =",p0_p)
		print("rho0/rho =",rho0_rho)
		print("\n")
	return [T0_T, p0_p, rho0_rho]

# normalShockCalculator returns the values of all property ratios across a normal shock given the upstream Mach number and gamma
# NOTE: printout is an optional feature which will print out the results to terminal.
def normalShockCalculator(M,gamma,printout):
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

# shockAreaIterator :: returns the geometric area at which a normal shock forms in a flow channel given proper conditions
# inputs :: Initial guess of the shock location with respect to the throat aream, As_At_guess
#			throat area, At [m2 or ft2]
#			upstream stagnation pressure, p0 [Pa or psf]
# 			back pressure, pb [Pa or psf]
# outputs :: shock area, As [m2 or ft2]
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

########## FANNO FLOW FUNCTIONS #########
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

########## RAYLEIGH FLOW FUNCTIONS #########
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

######### OBLIQUE SHOCK CALCULATORS ##########
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

    r = (gamma+1)*term / ((gamma-1)*term + 2)
    thetaMax = np.arctan((r-1)/(2*r**0.5))

    return b, M2, T_ratio, p_ratio, stagp_ratio, thetaMax

def obliqueShockCalculator2(beta, M, gamma):
    #betaRange = np.linspace(0.001,np.pi/2,1001) # brute force through all possible beta values
    #for b in betaRange: # loop through range of beta values until the RHS is larger than the LHS (requires computing/precision)
    def equation(x):
    	LHS = np.tan(x) # LHS is constant
    	RHS = (2/np.tan(beta))*(M**2*np.sin(beta)**2 - 1)/(M**2*(gamma+np.cos(2*beta)) + 2)
    	return RHS-LHS

    theta = fsolve(equation,0.2)[0]

    term = M**2*np.sin(beta)**2
    M2 = np.sqrt((1/np.sin(beta-theta)**2) * (2 + (gamma-1)*term)/(2*gamma*term - (gamma-1))) # from NASA website
    T_ratio = (2*gamma*term - (gamma-1))*((gamma-1)*term + 2)/((gamma+1)**2*term) # from NASA website
    p_ratio = (2*gamma*term - (gamma-1))/(gamma+1) # from NASA website
    term1 = pow(((gamma+1)*term)/((gamma-1)*term + 2),(gamma/(gamma-1)))
    term2 = pow(((gamma+1)/(2*gamma*term - (gamma-1))),(1/(gamma-1)))
    stagp_ratio = term1*term2
    return theta, M2, T_ratio, p_ratio, stagp_ratio

########## PRANDTL-MEYER EXPANSION CALCULATORS ##########
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

########## GENERALIZED 1D COMPRESSIBLE (INVISCID) FLOW TOOLS ##########
def findSonicLocation(flowClass):
	def equation(x):
		dT0_dx = flowClass.energyEquation(x)
		A_x, dA_dx = flowClass.flowArea(x)
		D_x, dD_dx = flowClass.flowDiameter(x)
		mdot_x, dmdot_dx = flowClass.massflowrate(x)
		return (-dA_dx/A_x + flowClass.gamma*4*flowClass.f/(2*D_x) + (1+flowClass.gamma)*dT0_dx/(2*(flowClass.T01 + x*dT0_dx)) + (1+flowClass.gamma)*dmdot_dx/mdot_x)
	
	x_sonic = fsolve(equation,flowClass.xguess)[0]
	dT0_dx = flowClass.energyEquation(x_sonic)
	T0_x = flowClass.T01 + x_sonic*dT0_dx
	A_x, dA_dx = flowClass.flowArea(x_sonic)
	D_x, dD_dx = flowClass.flowDiameter(x_sonic)
	mdot_x, dmdot_dx = flowClass.massflowrate(x_sonic)
	G,dGdx = Gfunction(x_sonic,1,flowClass)
	b = ((flowClass.gamma+1)/8)*2*flowClass.gamma*(4*flowClass.f/D_x + dT0_dx/T0_x + 2*dmdot_dx/mdot_x)
	c = ((flowClass.gamma+1)/8)*dGdx

	def equation2(x):
		return (x**2 + b*x + c)

	dM_dx_subsonic = fsolve(equation2, -0.5)
	dM_dx_supersonic = fsolve(equation2, 0.5)
	sonicPoints = [dM_dx_subsonic[0],dM_dx_supersonic[0]]

	return x_sonic,sonicPoints

def Gfunction(x,M,flowClass):
	z = sym.symbols('z')
	T0, dT0_dx = flowClass.energyEquation(x,z)
	#T0, dT0_dx = flowClass.T0x,flowClass.dT0dx
	#print(T0,dT0_dx)
	A_x, dA_dx = flowClass.flowArea(x,z)
	D_x, dD_dx = flowClass.flowDiameter(x,z)
	mdot_x, dmdot_dx = flowClass.massflowrate(x,z)
	G_sym = (-dA_dx/A_x + flowClass.gamma*M**2*4*flowClass.f/(2*D_x) + (1+flowClass.gamma*M**2)*dT0_dx/(2*T0) + (1+flowClass.gamma*M**2)*dmdot_dx/mdot_x) 
	dG_dx_sym = -2*sym.diff(dA_dx/A_x) + flowClass.gamma*sym.diff(4*flowClass.f/D_x) + (1+flowClass.gamma)*sym.diff(dT0_dx/T0) + 2*(1+flowClass.gamma)*sym.diff(dmdot_dx/mdot_x)
	if isinstance(G_sym,sym.Basic):
		G_num = float(G_sym.subs(z,x))
		dGdx_num = float(dG_dx_sym.subs(z,x))
		return G_num,dGdx_num
	else:
		return G_sym,dG_dx_sym

def generalized1Dflow(x,M,flowClass):
	if np.abs(M-1)<5e-2: # singularity case, if M is near 1	
		if flowClass.MachStart<1:
			dM_dx = flowClass.dMdx_sp[1] # positive solution -- M should increase towards and after the throat if flow starts subsonic
		elif flowClass.MachStart>1:
			dM_dx = flowClass.dMdx_sp[0] # negative solution -- M should decrease towards and after the throat if flow starts supersonic
		else:
			print("you can't specify a starting Mach number of 1 you fucking idiot")
	else: # general case
		psi_M = 1 + (flowClass.gamma-1)/2 * M**2
		dM_dx = (M*psi_M/(1-M**2)) * Gfunction(x,M,flowClass)[0]
	return dM_dx

def Gfunction_numeric(x,M,flowClass):
	z = "garbo"
	#T0, dT0_dx = flowClass.energyEquation(x,z)
	T0, dT0_dx = flowClass.T0x,flowClass.dT0dx
	#print(T0,dT0_dx)
	A_x, dA_dx = flowClass.flowArea(x,z)
	D_x, dD_dx = flowClass.flowDiameter(x,z)
	mdot_x, dmdot_dx = flowClass.massflowrate(x,z)
	G_num = (-dA_dx/A_x + flowClass.gamma*M**2*4*flowClass.f/(2*D_x) + (1+flowClass.gamma*M**2)*dT0_dx/(2*T0) + (1+flowClass.gamma*M**2)*dmdot_dx/mdot_x) 
	# dGdx_num isn't actually used (???) so i'm ignoring it because the sym.diff is throwing errors
	#dGdx_num = -2*sym.diff(dA_dx/A_x) + flowClass.gamma*sym.diff(4*flowClass.f/D_x) + (1+flowClass.gamma)*sym.diff(dT0_dx/T0) + 2*(1+flowClass.gamma)*sym.diff(dmdot_dx/mdot_x)
	return G_num

def generalized1Dflow_numeric(x,M,flowClass):
	if np.any((M-1)<1e-8): # singularity case, if M is near 1
		if flowClass.MachStart<1:
			dM_dx = flowClass.dMdx_sp[1] # positive solution -- M should increase towards and after the throat if flow starts subsonic
		elif flowClass.MachStart>1:
			dM_dx = flowClass.dMdx_sp[0] # negative solution -- M should decrease towards and after the throat if flow starts supersonic
		else:
			print("you can't specify a starting Mach number of 1 you fucking idiot")
	else: # general case
		psi_M = 1 + (flowClass.gamma-1)/2 * M**2
		dM_dx = (M*psi_M/(1-M**2)) * Gfunction_numeric(x,M,flowClass)
	return dM_dx

def stagPressureRK4(x,M,flowClass):
	z = "garbo"
	T0, dT0_dx = flowClass.T0x,flowClass.dT0dx
	D_x, dD_dx = flowClass.flowDiameter(x,z)
	mdot_x, dmdot_dx = flowClass.massflowrate(x,z)
	dP0_P0 = (-flowClass.gamma*M**2/2)*flowClass.gamma*M**2*4*flowClass.f/(2*D_x) - (-flowClass.gamma*M**2/2)*(dT0_dx/T0) - flowClass.gamma*M**2*(dmdot_dx/mdot_x)
	print(dP0_P0)
	return dP0_P0

def stagTempRK4(x,T0,flowClass):
	return (flowClass.deltaH_B/flowClass.Cp - flowClass.T0x) * (flowClass.massflowrate(x,"garbo")[1]/flowClass.massflowrate(x,"garbo")[0])



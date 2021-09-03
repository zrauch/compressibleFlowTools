import numpy as np
from numpy.linalg import inv
import scipy as sp
from scipy.signal import unit_impulse
import sympy as sym
import math
import sys,os
from scipy.optimize import fsolve
import spatial_discretization

## variables
global gamma, R, Cp, R_SI,CP_SI,R_ENG,CP_ENG
R_SI = 287 # J/kg-K -- m^2/s^2-K\
CP_SI = 1004.5 # J/kg-K
R_ENG = 53.353 # ft-lbf/lbm-R -- needs multiplied by g for most applications
#R_ENG = 1717 # ft^2/s^2-R
CP_ENG = 0.3 # Btu/lbm-R

global g,convertBTU,StoHR,IN2toFT2,FtoR
g = 32.174 # [ft-lbm/lbf-s2]
convertBTU = 778.17 # [ft-lbf/BTU]
StoHR = 3600 # [s/hr]
IN2toFT2 = 144 # [144 in2/ft2]
FtoR = 459.67

inputfound = 0; count = 0; limit = 4
unitsys = input("Are you working in the Metric (M) or Imperial (I) unit system?\t")
options = ["Metric","metric","M","English","english","E","Imperial","imperial","I","Ass","ass"]
if (unitsys in options):
	inputfound = 1

while (inputfound != 1):
	if (count == 0):
		print(unitsys,"is not a valid input, please enter one of the following:")
		print(options)
	unitsys = input("Are you working in the Metric (M) or Imperial (I) unit system?\t")
	count+=1
	if (unitsys in options):
		inputfound = 1
	elif (count == limit):
		print("\nHEY! \nstop being a dumbass and enter one of the values above. Last chance.")
	elif (count > limit):
		print("\nget outta here you filthy animal")
		exit(0)

if (options.index(unitsys) < 3):
	R = R_SI
	Cp = CP_SI
else:
	R = R_ENG
	Cp = CP_ENG

gamma = 1.4

################################################################################################### 
# enginesystems.py is a library intended to contain all relevant python codes developed for the
# complete solving of air-breathing propulsion systems such as turbofans, turbojets, A/B turbojets,
# and sc/ramjets.
# 
# most included routines are functions of M, stagnation properties, combustion product temperature,
# and the necessary parameters for each unique system.
################################################################################################### 


# turbofan :: a function to solve the complete IDEAL turbofan engine system.
# Note that the equations derived for use in this function depend on some critical assumptions,
# namely that the engine is perfectly expanded (pe = p0), and the flow rate of fuel is negligible
# compared to the mass flow rate of air through the engine (i.e fuel flow rate terms are neglected)
# inputs:
# 1. Mach number, M
# 2. free stream static pressure, p [Pa or psf]
# 3. free stream static temperature, T [K or R]
# 4. total pressure ratio across compressor, pi_c
# 5. total pressure ratio across fan, pi_f
# 6. engine bypass ratio, B
# outputs:
# 1. specific thrust, ST [s]
# 2. specific fuel consumption, SFC [check these units]
def turbofan(M, p, T, pi_c, pi_f, B, T04, delta_HB):
	## IDEAL Turbofan Calculations
	a0 = np.sqrt(gamma*R*g*T)
	u0 = M*a0
	# 0 --> 2 : compression through ideal diffuser
	if M<1:
		p02 = p*(1+((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	elif M>=1 and M<3:
		p02 = p*(1 - 0.075*(M-1)**1.35)
	
	T02 = T00 = T*(1+((gamma-1)/2)*M**2)
	tau_r = T00/T

	# 2 --> 3: compression through ideal compressor
	p03 = pi_c*p02
	T03 = T02*(pi_c)**((gamma-1)/gamma)
	tau_c = T03/T02

	# 2 --> 3' : compression through fan
	p03_p = pi_f*p02
	T03_p = T02*(pi_f)**((gamma-1)/gamma)
	tau_F = T03_p/T02

	# 3 --> 4: isobaric combustion process
	tau_b = T04/T03
	tau_lamba = tau_r*tau_c*tau_b
	p04 = p03
	f = Cp*(T04-T03)/(delta_HB-Cp*T04)
	#f = Cp*T00/delta_HB * (tau_lamba - tau_r*tau_c)

	# 4 --> 5: expansion through ideal turbine
	T05 = T04*(1 - (T02/T04)*((tau_c-1) + B*(tau_F-1)))
	tau_t = T05/T04
	p05 = p04*tau_t**((gamma-1)/gamma)

	# 5 -- 9: expansion through ideal nozzle
	p09 = p05
	T09 = T05

	## 3' --> 9': isentropic nozzle flow
	p09_p = p03_p
	T09_p = T03_p

	## Turbofan
	u9 = u0*np.sqrt(tau_b*(tau_r*tau_c*tau_t - 1)/(tau_r-1))
	M9_p = np.sqrt((2/(gamma-1))*(tau_r*tau_F - 1))
	T9_p = T09_p/(1+(gamma-1)/2*M9_p**2)
	a9 = np.sqrt(gamma*R*g*T9_p)
	u9_p = M9_p*a9

	ST = ((a0*M/(1+B))*((u9/u0 - 1) + B*(u9_p/u0 - 1)))
	SFC = f/((1+B)*ST)*g#*StoHR
	return [ST, SFC]

# turbojet :: a function to solve the complete IDEAL turbojet engine system.
# Note that the equations derived for use in this function depend on some critical assumptions,
# namely that the engine is perfectly expanded (pe = p0), and the flow rate of fuel is negligible
# compared to the mass flow rate of air through the engine (i.e fuel flow rate terms are neglected)
# inputs:
# 1. Mach number, M
# 2. free stream static pressure, p [Pa or psf]
# 3. free stream static temperature, T [K or R]
# 4. total pressure ratio across compressor, pi_c
# outputs:
# 1. specific thrust, ST [s]
# 2. specific fuel consumption, SFC [check these units]
def turbojet(M,p,T,pi_c,delta_HB):
	## IDEAL Turbofan Calculations
	a0 = np.sqrt(gamma*R*g*T)
	u0 = M*a0
	# 0 --> 2 : compression through ideal diffuser
	p02 = p00 = p*(1+((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	T02 = T00 = T*(1+((gamma-1)/2)*M**2)
	tau_r = T00/T

	# 2 --> 3: compression through ideal compressor
	p03 = pi_c*p02
	T03 = T02*(pi_c)**((gamma-1)/gamma)
	tau_c = T03/T02

	# 3 --> 4: isobaric combustion process
	tau_b = T04/T03
	tau_lamba = tau_r*tau_c*tau_b
	p04 = p03
	f = Cp*(T04-T03)/(delta_HB-Cp*T04)
	#f = Cp*T00/delta_HB * (tau_lamba - tau_r*tau_c)

	# 4 --> 5: expansion through ideal turbine
	T05 = T04*(1 - (T02/T04)*(tau_c-1))
	tau_t = T05/T04
	p05 = p04*tau_t**((gamma-1)/gamma)

	# 7 --> 9: Nozzle
	T09 = T05
	p09 = p05

	## Ideal Afterburning Turbojet
	#f_total = f+f_ab
	M9 = np.sqrt(2/(gamma-1)*(tau_r*tau_c*tau_t-1))
	T9 = T09/(1+(gamma-1)/2*M9**2)
	a9 = np.sqrt(gamma*R*T9*g)
	u9 = M9*a9
	ST = a0*M*((1+f)*u9/u0 -1)/g
	SFC = f/ST*StoHR
	return [ST, SFC]

def afterburningTurbojet(M,p,T,pi_c,T07):
	## IDEAL Turbofan Calculations
	a0 = np.sqrt(gamma*R*g*T)
	u0 = M*a0
	# 0 --> 2 : compression through ideal diffuser
	p02 = p00 = p*(1+((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	T02 = T00 = T*(1+((gamma-1)/2)*M**2)
	tau_r = T00/T

	# 2 --> 3: compression through ideal compressor
	p03 = pi_c*p02
	T03 = T02*(pi_c)**((gamma-1)/gamma)
	tau_c = T03/T02

	# 3 --> 4: isobaric combustion process
	tau_b = T04/T03
	tau_lamba = tau_r*tau_c*tau_b
	p04 = p03
	T03 = T02*(pi_c)**((gamma-1)/gamma)
	f = Cp*(T04-T03)/(delta_HB-Cp*T04)
	#f = Cp*T00/delta_HB * (tau_lamba - tau_r*tau_c)

	# 4 --> 5: expansion through ideal turbine
	T05 = T04*(1 - (T02/T04)*(tau_c-1))
	tau_t = T05/T04
	p05 = p04*tau_t**((gamma-1)/gamma)

	# 5 -- 7: Afterburner
	p07 = p05
	tau_ab = T07/T05
	tau_lamba_ab = tau_r*tau_c*tau_b*tau_t*tau_ab
	f_ab = Cp*T/delta_HB*(tau_lamba_ab - tau_lamba*tau_t)
	f_total = Cp*T/delta_HB*(tau_lamba_ab - tau_r)
	## NOTE: T07 fixed for the afterburner case (assumed to be maximum)

	# 7 --> 9: Nozzle
	T09 = T07
	# p09 not needed

	## Ideal Afterburning Turbojet
	#f_total = f+f_ab
	M9 = np.sqrt(2/(gamma-1)*(tau_r*tau_c*tau_t-1))
	T9 = T09/(1+(gamma-1)/2*M9**2)
	a9 = np.sqrt(gamma*R*T9*g)
	u9 = M9*a9
	ST = a0*M*((1+f_total)*u9/u0 -1)/g
	SFC = f_total/ST*StoHR
	return [ST, SFC]

def afterburningTurbojetWithBleeds(M,p,T,tpr,T07,eps1,eps2,eps3,p0ratio):
	# 0 : Ambient Stagnation Conditions
	a0 = np.sqrt(gamma*R*T*g)
	u0 = M*a0
	isen0 = 1 + ((gamma-1)/2)*M**2
	p00 = p*isen0**(gamma/(gamma-1))
	T00 = T*isen0
	# 0 --> 2: Diffuser
	#p02=p00
	p02 = p0ratio*p00
	T02 = T00 # same as isentropic case since diffuser is adiabatic
	tau_r = T02/T

	# 2 --> 3: Compressor
	p03 = p02*tpr;
	T03 = T02*(1 + (tpr**((gamma-1)/gamma)-1)/n_c);
	tau_c = T03/T02

	# 3 --> 4: Combustor
	p04 = p03 # no change in stagnation pressure
	f_act = (1-eps1-eps2-eps3)*Cp*(T04-T03)/(n_b*delta_HB-Cp*T04)
	tau_b = T04/T03
	tau_lamba = tau_r*tau_c*tau_b

	# 4 --> 5: Turbine
	T02_p = (T03+T02)/2
	T05 = T04 - (T02_p-T02 + (1-eps3)*(T03-T02_p))/(1-eps2-eps3+f_act)
	tau_t = T05/T04
	p05 = p04 * (1 - ((1 - tau_t)/ n_t))**(gamma/(gamma-1))

	# 5 --> 7: Afterburner
	p07 = p05 # SUPER SKETCHY ASSUMPTION
	f_ab = ((1-eps3)*Cp*(T07-T05)+Cp*eps3*(T07-T02_p))/(n_b*delta_HB - Cp*T07)
	f_total = f_act+f_ab
	## NOTE: T07 fixed for the afterburner case (assumed to be maximum)

	# 7 --> 9: Nozzle
	T09 = T07
	p09 = p07

	p0_list = [p00,p02,p02,p03,p04,p05,p05,p07,p07,p09]
	T0_list = [T00,T02,T02,T03,T04,T05,T05,T07,T07,T09]

	# Specific Fuel Consumption Calculations
	u9 = np.sqrt(2*Cp*g*convertBTU*n_n*T07*(1-((p/p09)**((gamma-1)/gamma))))
	M9 = np.sqrt(2/(gamma-1)*((p09/p)**((gamma-1)/gamma) - 1))
	ST= ((1+f_total)*u9-u0)/g
	SFC = f_total/ST*StoHR

	return [ST,SFC,T0_list,p0_list,M9,f_total,f_ab]

def ramjet(M,p,T):
	a0 = np.sqrt(gamma*R*g*T)
	u0 = M*a0
	# 0 --> 2 : compression through ideal diffuser
	p02 = p00 = p*(1+((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	T02 = T00 = T*(1+((gamma-1)/2)*M**2)
	tau_r = T00/T

	# 2 --> 3: compression through ideal compressor
	p03 = p02 # SKETCHY ASSUMPTIONS AAAAA!!!
	T03 = T02
	tau_c = T03/T02

	# 3 --> 4: isobaric combustion process
	tau_b = T04/T03
	tau_lamba = tau_r*tau_c*tau_b
	p04 = p03
	f = Cp*(T04-T03)/(delta_HB-Cp*T04)
	#f = Cp*T00/delta_HB * (tau_lamba - tau_r*tau_c)

	# 4 --> 5: expansion through ideal turbine
	T05 = T04*(1 - (T02/T04)*(tau_c-1))
	tau_t = T05/T04
	p05 = p04*tau_t**((gamma-1)/gamma)

	# 7 --> 9: Nozzle
	T09 = T05
	p09 = p05

	#f_total = f+f_ab
	M9 = np.sqrt(2/(gamma-1)*(tau_r*tau_c*tau_t-1))
	T9 = T09/(1+(gamma-1)/2*M9**2)
	a9 = np.sqrt(gamma*R*T9*g)
	u9 = M9*a9
	ST = a0*M*((1+f)*u9/u0 -1)/g
	SFC = f/ST*StoHR
	return [ST, SFC]

def nozzleExitArea(M9,T09,p09,mdot):
	T9 = T09/(1+(gamma-1)/2*M9**2)
	a9 = np.sqrt(gamma*R*g*T9)
	u9 = a9*M9
	rho9 = p/(R*g*T9)
	A9 = mdot/(rho9*u9)/g # calculate nozzle exit area based on M9
	return A9

def nozzleThroatArea(M8,T08,p08,mdot):
	T8 = T08/(1+(gamma-1)/2*M8**2)
	a8 = np.sqrt(gamma*R*g*T8)
	u8 = a8*M8
	p8 = p08/(1+(gamma-1)/2*M8**2)**(gamma/(gamma-1))
	rho8 = p8/(R*g*T8)
	A8 = mdot/(rho8*u8)/g # fix throat area for mdot_SLS
	return A8
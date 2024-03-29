import numpy as np
from numpy.linalg import inv
import scipy as sp
from scipy.signal import unit_impulse
import sympy as sym
import math
import sys,os
from scipy.optimize import fsolve
import spatial_discretization

# setunits defines unit system-dependent things like R, g, and Cp
from setunits import *

################################################################################################### 
# enginesystems.py is a library intended to contain all relevant python codes developed for the
# complete solving of air-breathing propulsion systems such as turbofans, turbojets, A/B turbojets,
# and sc/ramjets.
# 
# most included routines are functions of M, stagnation properties, combustion product temperature,
# and the necessary parameters for each unique system.
################################################################################################### 
def milstd_5008_B(M):
	if M<=1:
		p02_p0 = 1
	elif M>1 and M<5:
		p02_p0 = 1 - 0.075*(M-1)**1.35
	elif M>=5:
		p02_p0 = 800/(M**4 + 935)
	else:
		print("M cannot be negative!")
	return p02_p0

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
def turbofan(M,h, p, T, pi_c, pi_f, B, T04, delta_HB):
	## IDEAL Turbofan Calculations
	u0 = M*np.sqrt(gamma*R*T)
	# 0 --> 2 : compression through ideal diffuser
	p0 = p * (1 + ((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	T02 = T00 = T*(1+((gamma-1)/2)*M**2)
	p02 = p0*milstd_5008_B(M)
	# 2 --> 3: compression through ideal compressor
	p03 = pi_c*p02
	T03 = T02*(pi_c)**((gamma-1)/gamma)
	# 2 --> 3' : compression through fan
	p03_p = pi_f*p02
	T03_p = T02*(pi_f)**((gamma-1)/gamma)
	# 3 --> 4: isobaric combustion process
	p04 = p03
	f = Cp*(T04-T03)/(delta_HB-Cp*T04)
	# 4 --> 5: expansion through ideal turbine
	T05 = T04 - T03 + T02 - B*(T03_p - T02) # assuming constant cp across engine
	p05 = p04*(T05/T04)**(gamma/(gamma-1))
	# 5 -- 9: expansion through ideal nozzle
	p09 = p05
	T09 = T05
	## 3' --> 9': isentropic nozzle flow
	p09_p = p03_p
	T09_p = T03_p
	u9 = np.sqrt(convertBTU*g* 2*Cp*T09*(1 - (p/p09)**((gamma-1)/gamma) ))
	u9_p = np.sqrt(convertBTU*g* 2*Cp*T03_p*(1 - (p/p03_p)**((gamma-1)/gamma)))

	specificThrust = ((u9 - u0) + B*(u9_p - u0))/(1+B)
	Isp = (((1+f)*u9-u0) + B*(u9_p-u0))/(f*g)
	sfc = 1/Isp
	return [specificThrust,sfc,Isp,f]

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
def turbojet(M, p, T, pi_c, T04, delta_HB):
	## IDEAL turbojet Calculations
	a0 = np.sqrt(gamma*R*T)
	u0 = M*a0
	# 0 --> 2 : compression through ideal diffuser
	p0 = p*(1+((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	T02 = T00 = T*(1+((gamma-1)/2)*M**2)
	p02 = p0*milstd_5008_B(M)
	# 2 --> 3: compression through ideal compressor
	p03 = pi_c * p02
	T03 = T02 * pi_c**((gamma-1)/gamma)
	# 3 --> 4: isobaric combustion process
	p04 = p03
	f = Cp*(T04-T03)/delta_HB
	# 4 --> 5: expansion through ideal turbine
	T05 = T04 + (T02-T03)/(1+f)
	p05 = p04 * (T05/T04)**(gamma/(gamma-1))
	# 7 --> 9: Nozzle
	T09 = T05
	p09 = p05
	## Ideal Turbojet
	u9 = np.sqrt(convertBTU*g* 2*Cp*T09*(1 - (p/p09)**((gamma-1)/gamma)))
	ST = ((1+f)*u9 - u0)
	Isp = ST/(f*g)
	SFC = 1/Isp
	return [ST,SFC,Isp,f]

def afterburningTurbojet(M, p, T, pi_c, T04, T07, delta_HB):
	## IDEAL Afterburning turbojet Calculations
	a0 = np.sqrt(gamma*R*T)
	u0 = M*a0
	# 0 --> 2 : compression through ideal diffuser
	p0 = p*(1+((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	T02 = T00 = T*(1+((gamma-1)/2)*M**2)
	p02 = p0*milstd_5008_B(M)

	# 2 --> 3: compression through ideal compressor
	p03 = pi_c*p02
	T03 = T02*(pi_c)**((gamma-1)/gamma)

	# 3 --> 4: isobaric combustion process
	p04 = p03
	T03 = T02*(pi_c)**((gamma-1)/gamma)
	f = Cp*(T04-T03)/delta_HB

	# 4 --> 5: expansion through ideal turbine
	T05 = T04 + (T02-T03)/(1+f)
	p05 = p04*(T05/T04)**(gamma/(gamma-1))

	# 5 -- 7: Afterburner
	p07 = p05
	f_ab = (1+f)*Cp*(T07-T05)/delta_HB
	f_total = f + f_ab

	# 7 --> 9: Nozzle
	T09 = T07
	p09 = p07

	## Ideal Afterburning Turbojet
	u9 = np.sqrt(convertBTU*g * 2*Cp*T09*(1 - (p/p09)**((gamma-1)/gamma)))
	ST = ((1+f_total)*u9 - u0)
	SFC = f_total/ST
	Isp = ST/(f_total*g)
	return [ST,SFC,Isp,f]

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

def ramjet(M, p, T, T04, delta_HB):
	# ideal ramjet calculations
	a0 = np.sqrt(gamma*R*T)
	u0 = M*a0
	p0 = p*(1+((gamma-1)/2)*M**2)**(gamma/(gamma-1))
	T00 = T * (1+((gamma-1)/2)*M**2)
	# 0 --> 2: adiabatic mil-std inlet
	T02 = T00
	p02 = p0*milstd_5008_B(M)
	
	# 2 --> 3: isentropic ram compression
	#p04 = p02*(T04/T02)**((gamma-1)/gamma)
	#print(p04,p02)

	# 3 --> 4: isobaric compression
	f = Cp*(T04-T02)/delta_HB

	# 4 --> 9: perfectly expanded exhaust conditions
	u9 = np.sqrt(convertBTU*g * 2*Cp*T04*(1 - (p/p02)**((gamma-1)/gamma)))
	st = ((1+f)*u9 - u0)
	Isp = st/(f*g)
	sfc = 1/Isp
	return [st,sfc,Isp,f]

def nozzleExitArea(M9,T09,p09,mdot):
	T9 = T09/(1+(gamma-1)/2*M9**2)
	a9 = np.sqrt(gamma*R*T9)
	u9 = a9*M9
	rho9 = p/(R*g*T9)
	A9 = mdot/(rho9*u9)/g # calculate nozzle exit area based on M9
	return A9

def nozzleThroatArea(M8,T08,p08,mdot):
	T8 = T08/(1+(gamma-1)/2*M8**2)
	a8 = np.sqrt(gamma*R*T8)
	u8 = a8*M8
	p8 = p08/(1+(gamma-1)/2*M8**2)**(gamma/(gamma-1))
	rho8 = p8/(R*g*T8)
	A8 = mdot/(rho8*u8)/g # fix throat area for mdot_SLS
	return A8
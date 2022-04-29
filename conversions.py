import numpy as np
##### DEFINE CONVERSION CONSTANTS #####
R_U = 8.31446261815324  # [J/mol-K]
G_ENG = 32.174          # [ft/s2]
G_ENG_IN = G_ENG*12     # [in/s2]


# TIME CONVERSIONS
MIN_TO_S = 60           # 1 min =  60 s

# ANGULAR CONVERSIONS
REV_TO_RAD = 2*np.pi    # 1 rev = 2*pi rad

# DISTANCE CONVERSIONS
IN_TO_M = 0.0254        # 1 in = 0.0254 m
FT_TO_IN = 12           # 1 ft = 12 in
M_TO_FT = 3.28084       # 1 m = 3.28084 ft

# VOLUME CONVERSIONS
M3_TO_L = 1e-3          # 1 m3 = 1000 L

# MASS CONVERSIONS
KG_TO_LBM = 2.20462     # 1 kg = 2.20462 lbm

# DENSITY CONVERSIONS
LBMIN3_TO_LBMFT3 = 1728  # 1 lbm/in3 = lbm/ft3
LBMIN3_TO_KGM3 = 27679.9 # 1 lbm/in3 = 27679.9 kg/m3

# PRESSURE CONVERSIONS
ATM_TO_PSIA = 14.6959   # 1 atm = 14.6959 psia
BAR_TO_PA = 100000      # 1 bar = 100000 Pa
PSIA_TO_PA = 6894.76    # 1 psi = 6894.76 Pa
PSIA_TO_KPA = 6.89476   # 1 psi = 6.89476 kPa

# TEMPERATURE CONVERSIONS
K_TO_R = 1.8            # 1 K = 1.8 R

def F_TO_R(temp):
    temp += 459.67
    return temp

def R_TO_F(temp):
    temp -= 459.67
    return temp

# ENERGY CONVERSIONS
JKG_TO_FTLBFLBM = 0.334553                          # 1 J/kg = 0.334553 ft-lbf/lbm
KJKGK_TO_FTLBFLBMR = 185.862535173871               # 1 kJ/kg-K = 185.862535173871 ft-lbf/lbm-R
JKGK_TO_FTLBFLBMR = 0.185862535177479               # 1 J/kg-K = 0.185862535177479 ft-lbf/lbm-R
IN2S2_TO_FTLBFLBM = 1/(G_ENG*FT_TO_IN*FT_TO_IN)     # 1 in2/s2 = 1/(G_ENG*FT_TO_IN*FT_TO_IN) 

# POWER CONVERSIONS
HP_TO_FTLBFS = 550 # 1 hp = 550 ft-lbf/s

# FLOWRATE CONVERSIONS
M3S_TO_GPM = 15850.3    # 1 m3/s = 15850.3 gpm
IN3S_TO_GPM = 0.25974 # 1 in3/s = 0.25974 gpm
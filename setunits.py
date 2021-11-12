## variables
global gamma,Cp,R,g,g0,Pr,R_SI,CP_SI,R_ENG,CP_ENG
global convertBTU,StoHR,IN2toFT2
options = ["Metric","metric","M","English","english","E","Imperial","imperial","I","Ass","ass"]
hardcode = 1
# SI constants
R_SI =      287 # J/kg-K -- m^2/s^2-K
CP_SI =     1004.5 # J/kg-K
G_SI =      9.81 # m/s2
# English constants
R_ENG =     1716.55 # ft^2/s^2-R
CP_ENG =    0.3 # Btu/lbm-R
G_ENG =     32.17405 # ft-lbm/lbf-s2

# conversion factors
convertBTU =    778.169 # [ft-lbf/BTU]
StoHR =         3600 # [s/hr]
IN2toFT2 =      144 # [144 in2/ft2]

def FtoR(tempF):
    return tempF + 459.67

def CtoK(tempC):
    return tempC + 273.15


try:
    g
except NameError:
    if (hardcode): # this case is for when I want to run hardcoded values
        unitsys = "I"
        gamma = 1.4
        Cp = 0.24 # Btu/lbm-R
        R = convertBTU*G_ENG*Cp*(1 - (1/gamma)) # ft2/s2-R
        Pr = 0.72 # assume a constant Prandtl number for now
        g = g0 = G_ENG
    else:
        inputfound = 0
        count = 0; limit = 4
        unitsys = input("Are you working in the Metric (M) or Imperial (I) unit system?\t")

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
            g = g0 = G_SI
        else:
            R = R_ENG
            Cp = CP_ENG
            g = g0 = G_ENG
        Pr = 0.72 # assume a constant Prandtl number for now
else:
    print('setunits has already been called!')
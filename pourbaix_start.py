import matplotlib.pyplot as plt
import numpy as np

# Constants
constants = {
'R': 8.31451, # [J/(K*mol)],
'T': 298.16,  # [K],
'n': 2,       # Number of electrons transferred
'F': 96485    # [C/mol]
}

## THERMODYNAMIC DATA ##

# Gibbs free energy of formation of species at STP [REVISED POURBAIX DIAGRAM]
constants_deltaG_formation = {
'deltaG_Zn': 0,                     # Zn(s)
'deltaG_Zn^2+': -147.203*10**3,     # Zn^2+(aq)
'deltaG_Zn(OH)^+': -333.20*10**3,   # Zn(OH)^+(aq)
'deltaG_Zn(OH)2': -525.02*10**3,    # Zn(OH)2(s)  ----  
'deltaG_Zn(OH)2-eps':-555.82*10**3,  # Solid form
'deltaG_Zn(OH)3^-': -696.52*10**3,  # Zn(OH)3^-(aq)
'deltaG_Zn(OH)4^2-': -860.59*10**3, # Zn(OH)4^2-(aq)
'deltaG_ZnO': -320.479*10**3,       # ZnO
'deltaG_H^+': 0,                    # H^+(aq)
'deltaG_H2O': -237.1*10**3,         # H2O
'deltaG_OH^-':-157.2*10**3,         # OH^-(aq)
}

# Entropy of formation of species at STP [REBVISED POURBAIX DIAGRAM]
constants_S_formation = {
'S_Zn': 41.63,          # Zn(s)
'S_Zn^2+': -109.8,      # Zn^2+(aq)
'S_Zn(OH)^+': -24,      # Zn(OH)^+(aq)
'S_Zn(OH)2': 43,        # Zn(OH)2(s)
'S_Zn(OH)2-eps': 77,
'S_Zn(OH)3^-': 40,      # Zn(OH)3^-(aq)
'S_Zn(OH)4^2-': 15,     # Zn(OH)4^2-(aq)
'S_ZnO': 43.65,         # ZnO
'S_H^+': 0,             # H^+(aq)
'S_H2O': 70,  # H2O
'S_OH^-':-11,  # OH^-(aq)
}

# Molar heat capacity of species at STP [REVISED POURBAIX DIAGRAM]
constants_Cp = {
'Cp_Zn': [21.334, 11.648*10**(-3), 0.054*10**(6)],
'Cp_ZnO':[45.338, 7.289*10**(-3), -0.573*10**(6)],
'Cp_Zn(OH)2-eps':74.27,
'Cp_Zn^2+':-25.8,
'Cp_Zn(OH)^+': 10,
'Cp_Zn(OH)2': 70,
'Cp_Zn(OH)3^-': 94,
'Cp_Zn(OH)4^2-': -284,
'Cp_H^+': 0,
'Cp_H2O': 75,
'Cp_OH^-': -149
}



equilibrium_constants_2 = {'Zn(OH)': 10**4.4, 'ZnOH2_aq': 10**11.3, 'Zn(OH)3': 10**14.14, 'Zn(OH)4': 10**17.66,\
                            'Zn(NH3)': 10**2.37, 'Zn(NH3)2': 10**4.81, 'Zn(NH3)3': 10**7.31, 'Zn(NH3)4': 10**9.46,\
                            'Zn(NH3)(OH)': 10**9.23, 'Zn(NH3)2(OH)': 10**10.80, 'Zn(NH3)3(OH)': 10**12,  'Zn(NH3)(OH)2_aq': 10**13, 'Zn(NH3)2(OH)2_aq': 10**13.6, 'Zn(NH3)(OH)3': 10**14.50,\
                            'ZnOH2_sat': 10**(-14.82), 'ZnO': 10**(-15.96), 'ZnCO3': 10**(-10), \
                            'H2CO3': 10**(6.33), 'HCO3': 10**(9.56), 'pCO2': 10**(-1.55),\
                            'NH3': 10**(-9.246)}  # Initialize with a non-zero value

## Calculating E0 from thermodynamic data

# Gibbs free energy for reaction st STP
constants_deltaG = {
'deltaG_I': constants_deltaG_formation['deltaG_Zn'] - constants_deltaG_formation['deltaG_Zn^2+'],                                                          # Zn^2+ + 2e^- <--> Zn(s)
'deltaG_II': constants_deltaG_formation['deltaG_Zn'] + constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)^+']) + constants_deltaG_formation['deltaG_H^+'],      # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
'deltaG_III': constants_deltaG_formation['deltaG_Zn'] + 2*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)2-eps']+ 2*constants_deltaG_formation['deltaG_H^+']),   # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
'deltaG_IV': constants_deltaG_formation['deltaG_Zn'] + 3*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)3^-'] + 3*constants_deltaG_formation['deltaG_H^+']), # Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + 2H2O
'deltaG_V': constants_deltaG_formation['deltaG_Zn'] + 4*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)4^2-'] + 4*constants_deltaG_formation['deltaG_H^+']), # Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
'deltaG_VI': constants_deltaG_formation['deltaG_Zn'] + constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_ZnO']) + constants_deltaG_formation['deltaG_H^+'],           # ZnO + 2e^- + 2H^+ <--> Zn + H2O
'deltaG_VIII': constants_deltaG_formation['deltaG_Zn(OH)2-eps'] - (constants_deltaG_formation['deltaG_Zn^2+'] + 2*constants_deltaG_formation['deltaG_OH^-']),  # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
'deltaG_IX': constants_deltaG_formation['deltaG_Zn(OH)3^-'] - (constants_deltaG_formation['deltaG_Zn(OH)2-eps'] + constants_deltaG_formation['deltaG_OH^-']),  # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
'deltaG_X': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_Zn(OH)3^-'] + constants_deltaG_formation['deltaG_OH^-']),# [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
'deltaG_XI': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_Zn(OH)2-eps'] + 2*constants_deltaG_formation['deltaG_OH^-'])# [J/mol] - Zn(OH)2 +2OH^- <--> Zn(OH)4^2-
}

# Entropy of reaction at STP
constants_deltaS = {
'deltaS_I': constants_S_formation['S_Zn'] - constants_S_formation['S_Zn^2+'],                                                          # Zn^2+ + 2e^- <--> Zn(s)
'deltaS_II': constants_S_formation['S_Zn'] + constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)^+']) + constants_S_formation['S_H^+'],      # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
'deltaS_III': constants_S_formation['S_Zn'] + 2*constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)2']+ 2*constants_S_formation['S_H^+']),   # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
'deltaS_IV': constants_S_formation['S_Zn'] + 3*constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)3^-'] + 3*constants_S_formation['S_H^+']), # Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + 2H2O
'deltaS_V': constants_S_formation['S_Zn'] + 4*constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)4^2-'] + 4*constants_S_formation['S_H^+']), # Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
'deltaS_VI': constants_S_formation['S_Zn'] + constants_S_formation['S_H2O'] - (constants_S_formation['S_ZnO']) + constants_S_formation['S_H^+'],           # ZnO + 2e^- + 2H^+ <--> Zn + H2O
'deltaS_VIII': constants_S_formation['S_Zn(OH)2'] - (constants_S_formation['S_Zn^2+'] + 2*constants_S_formation['S_OH^-']),  # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
'deltaS_IX': constants_S_formation['S_Zn(OH)3^-'] - (constants_S_formation['S_Zn(OH)2'] + constants_S_formation['S_OH^-']),  # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
'deltaS_X': constants_S_formation['S_Zn(OH)4^2-'] - (constants_S_formation['S_Zn(OH)3^-'] + constants_S_formation['S_OH^-']),# [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
'deltaS_XI': constants_S_formation['S_Zn(OH)4^2-'] - (constants_S_formation['S_Zn(OH)2'] + 2*constants_S_formation['S_OH^-'])# [J/mol] - Zn(OH)2 +2OH^- <--> Zn(OH)4^2-
}

# Enthalpy for reaction as STP
constants_deltaH = {
'deltaH_I': constants_deltaG['deltaG_I'] + constants['T']*constants_deltaS['deltaS_I'],         # Zn^2+ + 2e^- <--> Zn(s)
'deltaH_II': constants_deltaG['deltaG_II'] + constants['T']*constants_deltaS['deltaS_II'],      # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
'deltaH_III': constants_deltaG['deltaG_III'] + constants['T']*constants_deltaS['deltaS_III'],   # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
'deltaH_IV': constants_deltaG['deltaG_IV'] + constants['T']*constants_deltaS['deltaS_IV'],      # Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + 2H2O
'deltaH_V': constants_deltaG['deltaG_V'] + constants['T']*constants_deltaS['deltaS_V'],         # Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
'deltaH_VI': constants_deltaG['deltaG_VI'] + constants['T']*constants_deltaS['deltaS_VI'],      # ZnO + 2e^- + 2H^+ <--> Zn + H2O
'deltaH_VIII': constants_deltaG['deltaG_VIII'] - constants['T']*constants_deltaS['deltaS_VIII'],# [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
'deltaH_IX': constants_deltaG['deltaG_IX'] - constants['T']*constants_deltaS['deltaS_IX'],      # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
'deltaH_X': constants_deltaG['deltaG_X'] - constants['T']*constants_deltaS['deltaS_X'],         # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
'deltaH_XI': constants_deltaG['deltaG_XI'] - constants['T']*constants_deltaS['deltaS_XI']       # [J/mol] - Zn(OH)2 +2OH^- <--> Zn(OH)4^2-
}

# Molar heat capacity for reaction
constants_deltaCp = {
'deltaCp_I': constants_Cp['Cp_Zn'][0] - constants_Cp['Cp_Zn^2+'],                                                           # Zn^2+ + 2e^- <--> Zn(s)                   ----- There are two more coefficients in Cp for Zn that must be considered when taking the integral over temperature
'deltaCp_II': constants_Cp['Cp_Zn'][0] + constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)^+']) + constants_Cp['Cp_H^+'],      # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O     ----- There are two more coefficients in Cp for Zn that must be considered when taking the integral over temperature
'deltaCp_III': constants_Cp['Cp_Zn'][0] + 2*constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)2']+ 2*constants_Cp['Cp_H^+']),   # Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O     ----- There are two more coefficients in Cp for Zn that must be considered when taking the integral over temperature
'deltaCp_IV': constants_Cp['Cp_Zn'][0] + 3*constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)3^-'] + 3*constants_Cp['Cp_H^+']), # Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + 2H2O   ----- There are two more coefficients in Cp for Zn that must be considered when taking the integral over temperature
'deltaCp_V': constants_Cp['Cp_Zn'][0] + 4*constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)4^2-'] + 4*constants_Cp['Cp_H^+']), # Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O----- There are two more coefficients in Cp for Zn that must be considered when taking the integral over temperature
'deltaCp_VI': constants_Cp['Cp_Zn'][0] + constants_Cp['Cp_H2O'] - (constants_Cp['Cp_ZnO'][0]) + constants_Cp['Cp_H^+'],           # ZnO + 2e^- + 2H^+ <--> Zn + H2O           ----- There are two more coefficients in Cp for Zn that must be considered when taking the integral over temperature
'deltaCp_VIII': constants_Cp['Cp_Zn(OH)2'] - (constants_Cp['Cp_Zn^2+'] + 2*constants_Cp['Cp_OH^-']),  # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
'deltaCp_IX': constants_Cp['Cp_Zn(OH)3^-'] - (constants_Cp['Cp_Zn(OH)2'] + constants_Cp['Cp_OH^-']),  # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
'deltaCp_X': constants_Cp['Cp_Zn(OH)4^2-'] - (constants_Cp['Cp_Zn(OH)3^-'] + constants_Cp['Cp_OH^-']),# [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
'deltaCp_XI': constants_Cp['Cp_Zn(OH)4^2-'] - (constants_Cp['Cp_Zn(OH)2'] + 2*constants_Cp['Cp_OH^-'])# [J/mol] - Zn(OH)2 +2OH^- <--> Zn(OH)4^2-    
}

## Calculation of Standard reduction potentials for reactions
# Zn^2+ + 2e^- <--> Zn(s)
EI_0 = -(constants_deltaG_formation['deltaG_Zn'] - constants_deltaG_formation['deltaG_Zn^2+'])/(constants['n']*constants['F'])

# Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
EII_0 = -(constants_deltaG_formation['deltaG_Zn'] + constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)^+']) + constants_deltaG_formation['deltaG_H^+'])/(constants['n']*constants['F'])

# Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + 2H2O
EIII_0 = -(constants_deltaG_formation['deltaG_Zn'] + 2*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)2-eps']+ 2*constants_deltaG_formation['deltaG_H^+']))/(constants['n']*constants['F'])

# Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
EIV_0 = -(constants_deltaG_formation['deltaG_Zn'] + 3*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)3^-'] + 3*constants_deltaG_formation['deltaG_H^+']))/(constants['n']*constants['F'])

# Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
EV_0 = -(constants_deltaG_formation['deltaG_Zn'] + 4*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)4^2-'] + 4*constants_deltaG_formation['deltaG_H^+']))/(constants['n']*constants['F'])

# ZnO + 2e^- + 2H^+ <--> Zn + H2O
EVI_0 = -(constants_deltaG_formation['deltaG_Zn'] + constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_ZnO']) + constants_deltaG_formation['deltaG_H^+'])/(constants['n']*constants['F'])

## E0 constants
constants_E0 = {
'EI_0': EI_0,      # [V] - Zn^2+(aq) --> Zn(s)
'EII_0': EII_0,    # [V] - Zn(OH)^+ --> Zn(s) - Not used
'EIII_0': EIII_0,  # [V] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
'EIV_0': EIV_0,    # [V] - Zn(OH)3^- --> Zn(s)
'EV_0': EV_0,      # [V] - Zn(OH)4^2- --> Zn(s)
'EVI_0': EVI_0,    # [V] - ZnO --> Zn(s) - Not used
'EHER_0': 0,       # [V] - HER
'EOER_0': 1.229    # [V] - OER
}

# The pZn where the equilibrium Zn(OH)2 <--> Zn(OH)3^-1 and Zn(OH)3^-1 <--> Zn(OH)4^2- becomes the same and the domain for Zn(OH)3^- vanishes
pZn_threshold = (constants_deltaG['deltaG_IX']-constants_deltaG['deltaG_X'])/(constants['R']*constants['T']*np.log(10))

## Plotting

pH = np.arange(0, 16, 0.01)     # pH range
EHER = constants_E0['EHER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH   # HER
EOER = constants_E0['EOER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH   # OER
#slope = -2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)
#angle = np.arctan(slope)

# Making the figure
plt.figure()
plt.plot(pH,EHER, '--')     # Line for the HER
plt.plot(pH,EOER, '--')     # Line for the OER

# Concentration of Zn ions
pZn_values = np.array([4, 6])           # Values we iterate through for the activity of dissolved Zn - pZn = -log(c_Zn) --- Maximum value is 6 and minimum value is 0
pZn_values_dict = {'pZn': pZn_values}   # Dictionary with values
linestyle = ['-', '--']

x_Zn2_ZnOH2 = []
x_ZnOH2_ZnOH4_2 = []
x_ZnOH2_ZnOH3_1 = []
x_ZnOH3_1_ZnOH4_2 = []

for i in range(len(pZn_values)):
    pZn_string = f'pZn =  {int(pZn_values_dict['pZn'][i])}' # The string printed for the different pZn values in the plot. 

    # The amount of Zn in solution
    constants_p = {             # Setting all values to be the same because they are the dominating species within their domain
    'pZn': pZn_values[i],       # pZn = -log(Zn^2+)
    'pZnOH': pZn_values[i],     # pZnOH = -log(Zn(OH)^+)
    'pZnOH2': pZn_values[i],    # pZnOH2 = -log(Zn(OH)2)
    'pZnOH3': pZn_values[i],    # pZnOH3 = -log(Zn(OH)3^-)
    'pZnOH4': pZn_values[i]     # pZnOH4 = -log(Zn(OH)4^-2)
    }

    ## E calculations

    # Zn^2+ --> Zn(s)
    EI = np.ones(len(pH))*constants_E0['EI_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZn'])
    
    # Zn(OH)^+ --> Zn(s)
    EII = constants_E0['EII_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH']+pH)

    # Zn(OH)2 --> Zn(s)
    EIII = constants_E0['EIII_0'] - 2*constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(0*constants_p['pZnOH2']/2+pH)

    # Zn(OH)3^- --> Zn(s)
    EIV = constants_E0['EIV_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH3']+3*pH)

    # Zn(OH)4^2- --> Zn(s)
    EV = constants_E0['EV_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH4']+4*pH)

    # ZnO --> Zn(s)
    EVI = constants_E0['EVI_0'] - 2*constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(pH)

    ## pH calculations

    # Zn^2+ --> Zn(OH)2
    pHVIII = 14 + constants_deltaG['deltaG_VIII']/(2*constants['R']*constants['T']*np.log(10)) + constants_p['pZn']/2 - 0*constants_p['pZnOH2']/2
    x_Zn2_ZnOH2.append(pHVIII)

    # Zn(OH)2 --> Zn(OH)3^-
    pHIX = 14 + constants_deltaG['deltaG_IX']/(constants['R']*constants['T']*np.log(10)) + 0*constants_p['pZnOH2'] - constants_p['pZnOH3']
    x_ZnOH2_ZnOH3_1.append(pHIX)
    # Zn(OH)3^- --> Zn(OH)4^2-
    pHX = 14 + constants_deltaG['deltaG_X']/(constants['R']*constants['T']*np.log(10)) + constants_p['pZnOH3'] - constants_p['pZnOH4']
    x_ZnOH3_1_ZnOH4_2.append(pHX)
    # Zn(OH)2 --> Zn(OH)4^2-
    pHXI = 14 + constants_deltaG['deltaG_XI']/(2*constants['R']*constants['T']*np.log(10)) - (1/2)*constants_p['pZnOH4']
    x_ZnOH2_ZnOH4_2.append(pHXI)

    # Finding intersections for the different lines
    index_I_III = np.where(EI <= EIII)[0][-1] if np.any(EI <= EIII) else None       # pH between Zn^2+ and Zn(OH)2 equilibrium
    index_III_IV = np.where(EIII <= EIV)[0][-1] if np.any(EIII <= EIV) else None    # pH between Zn(OH)2 and Zn(OH)3^- equilibrium
    index_IV_V = np.where(EIV <= EV)[0][-1] if np.any(EIV <= EV) else None          # pH between Zn(OH)3^-1 and Zn(OH)4^2- equilibrium
    index_III_V = np.where(EIII <= EV)[0][-1] if np.any(EIII <=EV) else None        # pH where Zn(OH)2 and Zn(OH)4^-2 equilibrium

    if any(pZn_values > pZn_threshold): # CXhecking if any of the elements pZn values are bigger than tyhe threshold


        # Plotting the lines of the Pourbaix diagram
        plt.plot(pH[:index_I_III], EI[:index_I_III], 'k', label='Zn$^{2+}$ - Zn', linestyle=linestyle[i])                               # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        plt.vlines(pHVIII, EI[index_I_III], 1.5, 'k', label='Zn$^{2+}$ - Zn(OH)_{2}', linestyle=linestyle[i])                           # Chemical        - Equilibrium between Zn^2+ and Zn(OH)2
        
        

        # Adding pZn text --> Should be adjusted later
        plt.text(pH[index_I_III+10], 1, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        
        
        if pZn_values[i] >= pZn_threshold:
            plt.vlines(pHX, EV[index_IV_V], 1.5, 'm', label='Zn(OH)$_{3}^-$ - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[i])                   # Chemical        - Equilibrium between Zn(OH)3^- and Zn(OH)4^2-
            plt.vlines(pHIX, EIV[index_III_IV], 1.5, 'r', label='Zn(OH)$_{2}$ - Zn(OH)$_{3}^-$', linestyle=linestyle[i])                    # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)3^-
            plt.plot(pH[index_III_IV:index_IV_V], EIV[index_III_IV:index_IV_V], 'r', label='Zn(OH)$_{3}^{-}$ - Zn', linestyle=linestyle[i]) # Electrochemical - Equilibrium between Zn(OH)3^- and Zn
            plt.text(pH[index_III_IV+10], 1, pZn_string, color='r', rotation = 90)  # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)3^- -- must be fixed
            plt.text(pH[index_IV_V+10], 1, pZn_string, color='m', rotation = 90)    # Adding pZn value to equilibrium between Zn(OH)3^-1 and Zn(OH)4^2 -- must be fixed
            plt.plot(pH[index_IV_V:], EV[index_IV_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[i])                         # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn   
            plt.plot(pH[index_I_III: index_III_IV], EIII[index_I_III: index_III_IV], 'c', label='Zn(OH)$_{2}$ - Zn')                          # Electrochemical - Equilibrium between Zn(OH)2 and Zn
            
        elif pZn_values[i]<pZn_threshold:
            plt.plot(pH[index_III_V:], EV[index_III_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[i])                       # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn
            plt.vlines(pHXI, EV[index_III_V], 1.5, 'm', label='Zn(OH)$_{2}$ - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[i])                   # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)4^2-
            plt.text(pH[index_III_V+10], 1, pZn_string, color='m', rotation = 90)    # Adding pZn value to equilibrium between Zn(OH)3^-1 and Zn(OH)4^2 -- must be fixed
            plt.plot(pH[index_I_III: index_III_V], EIII[index_I_III: index_III_V], 'c', label='Zn(OH)$_{2}$ - Zn')                          # Electrochemical - Equilibrium between Zn(OH)2 and Zn

        if i == len(pZn_values)-1:
            # Adding text to different domains --> Should be adjusted later
            plt.text(5, -1.25, 'Zn(s)', fontsize=12, color='black')                                 # Zn domain
            plt.text(np.min(x_Zn2_ZnOH2)-3, 0.45, 'Zn$^{2+}$(aq)', fontsize=12, color='black')                           # Zn^2+ domain
            plt.text(2, -0.32, 'HER', fontsize=12, color='black', rotation=-10)                     # HER line
            plt.text(2, 0.9, 'OER', fontsize=12, color='black', rotation=-10)                       # OER line
            plt.text(np.max([np.max(x_ZnOH3_1_ZnOH4_2), np.max(x_ZnOH2_ZnOH4_2)])+0.5, -0.65, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain
            plt.text(np.min(x_ZnOH2_ZnOH3_1)+0.2, -0.4, 'Zn(OH)$^{-}_{3}$(aq)', fontsize=12, color='black', rotation=90)     # Zn(OH)3^- domain
            plt.text(np.max(x_Zn2_ZnOH2)+1, -0.25, 'Zn(OH)$_2$(s)', fontsize=12, color='black', rotation=90)              # Zn(OH)2 domain
        
    else:
        # Plotting the lines of the Pourbaix diagram
        plt.plot(pH[:index_I_III], EI[:index_I_III], 'k', label='Zn$^{2+}$ - Zn', linestyle=linestyle[i])                               # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        plt.vlines(pHVIII, EI[index_I_III], 1.5, 'k', label='Zn$^{2+}$ - Zn(OH)_{2}', linestyle=linestyle[i])                           # Chemical        - Equilibrium between Zn^2+ and Zn(OH)2
        plt.plot(pH[index_I_III: index_III_V], EIII[index_I_III: index_III_V], 'c', label='Zn(OH)$_{2}$ - Zn')# Electrochemical - Equilibrium between Zn(OH)2 and Zn
        plt.plot(pH[index_III_V:], EV[index_III_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[i])                       # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn
        plt.vlines(pHXI, EV[index_III_V], 1.5, 'm', label='Zn(OH)$_{2}$ - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[i])                   # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)4^2-
        
        # Adding pZn text --> Should be adjusted later
        plt.text(pH[index_I_III+10], 1, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        plt.text(pH[index_III_V+10], 1, pZn_string, color='m', rotation = 90)   # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)4^2 -- must be fixed

        if i == len(pZn_values)-1:
            # Adding text to different domains --> Should be adjusted later
            plt.text(5, -1.25, 'Zn(s)', fontsize=12, color='black')                                 # Zn domain
            plt.text(np.min(x_Zn2_ZnOH2)-3, 0.45, 'Zn$^{2+}$(aq)', fontsize=12, color='black')      # Zn^2+ domain
            plt.text(2, -0.32, 'HER', fontsize=12, color='black', rotation=-10)                     # HER line
            plt.text(2, 0.9, 'OER', fontsize=12, color='black', rotation=-10)                       # OER line
            plt.text(np.max(x_ZnOH2_ZnOH4_2)+0.5, -0.65, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain
            plt.text(np.max(x_Zn2_ZnOH2)+1, -0.25, 'Zn(OH)$_2$(s)', fontsize=12, color='black', rotation=0)              # Zn(OH)2 domain 

plt.xlabel('pH  /  []')
plt.ylabel('Potential - E  /  V')
plt.xlim(xmin=0, xmax=max(pH))  # Set the x-axis range
plt.ylim(ymin=-1.5, ymax=1.5)
plt.show()
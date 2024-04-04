import matplotlib.pyplot as plt
import numpy as np

# Constants
constants = {
'R': 8.31451, # [J/(K*mol)],
'T': 298.16,  # [K],
'n': 2,       # Number of electrons transferred
'F': 96485    # [C/mol]
}

# Gibbs free energy of formation of species at STP [REVISED POURBAIX DIAGRAM]
constants_deltaG_formation = {
'deltaG_Zn': 0,                     # Zn(s)
'deltaG_Zn^2+': -147.203*10**3,     # Zn^2+(aq)
'deltaG_Zn(OH)^+': -333.20*10**3,   # Zn(OH)^+(aq)
'deltaG_Zn(OH)2': -555.82*10**3,    # Zn(OH)2(s)
'deltaG_Zn(OH)3^-': -696.52*10**3,  # Zn(OH)3^-(aq)
'deltaG_Zn(OH)4^2-': -860.59*10**3, # Zn(OH)4^2-(aq)
'deltaG_ZnO': -320.479*10**3,       # ZnO
'deltaG_H^+': 0,                    # H^+(aq)
'deltaG_H2O': -237.1*10**3,         # H2O
'deltaG_OH^-':-157.2*10**3,         # OH^-(aq)
}

## Calculating E0 from thermodynamic data

# Zn^2+ + 2e^- <--> Zn(s)
EI_0 = -(constants_deltaG_formation['deltaG_Zn'] - constants_deltaG_formation['deltaG_Zn^2+'])/(constants['n']*constants['F'])

# Zn(OH)^+ + 2e- + H^+ <--> Zn(s) + H2O
EII_0 = -(constants_deltaG_formation['deltaG_Zn'] + constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)^+']) + constants_deltaG_formation['deltaG_H^+'])/(constants['n']*constants['F'])

# Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + 2H2O
EIII_0 = -(constants_deltaG_formation['deltaG_Zn'] + 2*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)2']+2*constants_deltaG_formation['deltaG_H^+']))/(constants['n']*constants['F'])

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

# deltaG for reactions
constants_deltaG = {
'deltaG_VIII': constants_deltaG_formation['deltaG_Zn(OH)2'] - (constants_deltaG_formation['deltaG_Zn^2+'] + 2*constants_deltaG_formation['deltaG_OH^-']),  # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
'deltaG_IX': constants_deltaG_formation['deltaG_Zn(OH)3^-'] - (constants_deltaG_formation['deltaG_Zn(OH)2'] + constants_deltaG_formation['deltaG_OH^-']),  # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
'deltaG_X': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_Zn(OH)3^-'] + constants_deltaG_formation['deltaG_OH^-']),# [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
'deltaG_XI': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_Zn(OH)2'] + 2*constants_deltaG_formation['deltaG_OH^-'])# [J/mol] - Zn(OH)2 +2OH^- <--> Zn(OH)4^2-
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
pZn_values = np.array([2, 6])           # Values we iterate through for the activity of dissolved Zn - pZn = -log(c_Zn)
pZn_values_dict = {'pZn': pZn_values}   # Dictionary with values

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

    # Zn(OH)2 --> Zn(OH)3^-
    pHIX = 14 + constants_deltaG['deltaG_IX']/(constants['R']*constants['T']*np.log(10)) + 0*constants_p['pZnOH2'] - constants_p['pZnOH3']

    # Zn(OH)3^- --> Zn(OH)4^2-
    pHX = 14 + constants_deltaG['deltaG_X']/(constants['R']*constants['T']*np.log(10)) + constants_p['pZnOH3'] - constants_p['pZnOH4']

    # Zn(OH)2 --> Zn(OH)4^2-
    pHXI = 14 + constants_deltaG['deltaG_XI']/(2*constants['R']*constants['T']*np.log(10)) - (1/2)*constants_p['pZnOH4']

    # Finding intersections for the different lines
    index_I_III = np.where(EI <= EIII)[0][-1] if np.any(EI <= EIII) else None       # pH between Zn^2+ and Zn(OH)2 equilibrium
    index_III_IV = np.where(EIII <= EIV)[0][-1] if np.any(EIII <= EIV) else None    # pH between Zn(OH)2 and Zn(OH)3^- equilibrium
    index_IV_V = np.where(EIV <= EV)[0][-1] if np.any(EIV <= EV) else None          # pH between Zn(OH)3^-1 and Zn(OH)4^2- equilibrium
    index_III_V = np.where(EIII <= EV)[0][-1] if np.any(EIII <=EV) else None        # pH where Zn(OH)2 and Zn(OH)4^-2 equilibrium

    ## When pZn reaches a certain value, the domain for Zn(OH)3^- vanishes
    if pZn_values[i] >= pZn_threshold:  # If the pZn is above the threshold then we have a domain for Zn(OH)3^-
        
        # Plotting the lines of the Pourbaix diagram
        plt.plot(pH[:index_I_III], EI[:index_I_III], 'k', label='Zn$^{2+}$ - Zn')                               # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        plt.vlines(pHVIII, EI[index_I_III], 1.5, 'k', label='Zn$^{2+}$ - Zn(OH)_{2}')                           # Chemical        - Equilibrium between Zn^2+ and Zn(OH)2
        plt.plot(pH[index_I_III: index_III_IV], EIII[index_I_III: index_III_IV], 'c', label='Zn(OH)$_{2}$ - Zn')# Electrochemical - Equilibrium between Zn(OH)2 and Zn
        plt.vlines(pHIX, EIV[index_III_IV], 1.5, 'r', label='Zn(OH)$_{2}$ - Zn(OH)$_{3}^-$')                    # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)3^-
        plt.plot(pH[index_III_IV:index_IV_V], EIV[index_III_IV:index_IV_V], 'r', label='Zn(OH)$_{3}^{-}$ - Zn') # Electrochemical - Equilibrium between Zn(OH)3^- and Zn
        plt.vlines(pHX, EV[index_IV_V], 1.5, 'm', label='Zn(OH)$_{3}^-$ - Zn(OH)$_{4}^{2-}$')                   # Chemical        - Equilibrium between Zn(OH)3^- and Zn(OH)4^2-
        plt.plot(pH[index_IV_V:], EV[index_IV_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn')                         # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn                     
        
        # Adding pZn text --> Should be adjusted later
        plt.text(pH[index_I_III+10], 1, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        plt.text(pH[index_III_IV+10], 1, pZn_string, color='r', rotation = 90)  # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)3^- -- must be fixed
        plt.text(pH[index_IV_V+10], 1, pZn_string, color='m', rotation = 90)    # Adding pZn value to equilibrium between Zn(OH)3^-1 and Zn(OH)4^2 -- must be fixed

        # Adding text to different domains --> Should be adjusted later
        plt.text(2, -1.25, 'Zn(s)', fontsize=12, color='black')                                 # Zn domain
        plt.text(2, 0.5, 'Zn$^{2+}$(aq)', fontsize=12, color='black')                           # Zn^2+ domain
        plt.text(2, -0.32, 'HER', fontsize=12, color='black', rotation=-10)                     # HER line
        plt.text(2, 0.9, 'OER', fontsize=12, color='black', rotation=-10)                       # OER line
        plt.text(14, -0.55, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain
        plt.text(12, -0.4, 'Zn(OH)$^{-}_{3}$(aq)', fontsize=12, color='black', rotation=90)     # Zn(OH)3^- domain
        plt.text(8, 0.0, 'Zn(OH)$_2$(s)', fontsize=12, color='black', rotation=90)              # Zn(OH)2 domain

    else:# # If the pZn is below the threshold then we don't have a domain for Zn(OH)3^- 
        # Plotting the lines of the Pourbaix diagram
        plt.plot(pH[:index_I_III], EI[:index_I_III], 'k', label='Zn$^{2+}$ - Zn')                               # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        plt.vlines(pHVIII, EI[index_I_III], 1.5, 'k', label='Zn$^{2+}$ - Zn(OH)_{2}')                           # Chemical        - Equilibrium between Zn^2+ and Zn(OH)2
        plt.plot(pH[index_I_III: index_III_V], EIII[index_I_III: index_III_V], 'c', label='Zn(OH)$_{2}$ - Zn')# Electrochemical - Equilibrium between Zn(OH)2 and Zn
        plt.plot(pH[index_III_V:], EV[index_III_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn')                       # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn
        plt.vlines(pHXI, EV[index_III_V], 1.5, 'm', label='Zn(OH)$_{2}$ - Zn(OH)$_{4}^{2-}$')                   # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)4^2-
        
        # Adding pZn text --> Should be adjusted later
        plt.text(pH[index_I_III+10], 1, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        plt.text(pH[index_III_V+10], 1, pZn_string, color='m', rotation = 90)   # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)4^2 -- must be fixed

        # Adding text to different domains --> Should be adjusted later
        plt.text(2, -1.25, 'Zn(s)', fontsize=12, color='black')                                 # Zn domain
        plt.text(2, 0.5, 'Zn$^{2+}$(aq)', fontsize=12, color='black')                           # Zn^2+ domain
        plt.text(2, -0.32, 'HER', fontsize=12, color='black', rotation=-10)                     # HER line
        plt.text(2, 0.9, 'OER', fontsize=12, color='black', rotation=-10)                       # OER line
        plt.text(14, -0.55, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain
        plt.text(8, 0.0, 'Zn(OH)$_2$(s)', fontsize=12, color='black', rotation=90)              # Zn(OH)2 domain

plt.xlabel('pH')
plt.ylabel('Potential [V]')
plt.xlim(xmin=0, xmax=max(pH))  # Set the x-axis range
plt.ylim(ymin=-1.5, ymax=1.5)
plt.show()
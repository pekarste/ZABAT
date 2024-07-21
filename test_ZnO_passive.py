import matplotlib.pyplot as plt
import numpy as np
from Functions.Functions import vant_Hoff, deltaG_weak, deltaG_T2, E0_2

# Constants
constants = {
    'R': 8.31451,       # [J/(K*mol)]   - The universal gas constant
    'T': 25+273.15,     # [K]           - Standard temperature
    'n': 2,             # [-]           - Number of electrons transferred
    'F': 96485          # [C/mol]       - Faraday's constant
}

T = 25+273.15           # [K] - Temperature

## THERMODYNAMIC DATA ##

# Gibbs free energy of formation of species at STP [REVISED POURBAIX DIAGRAM]
constants_deltaG_formation = {
    'deltaG_Zn': 0,                     # [J/mol] - Zn(s)
    'deltaG_Zn^2+': -147.203*10**3,     # [J/mol] - Zn^2+(aq)
    'deltaG_Zn(OH)^+': -333.20*10**3,   # [J/mol] - Zn(OH)^+(aq)    -- Soluble hydroxide
    'deltaG_Zn(OH)2': -525.02*10**3,    # [J/mol] - Zn(OH)2(aq)     -- Solid hydroxide
    'deltaG_Zn(OH)2-eps':-555.82*10**3, # [J/mol] - eps-Zn(OH)2(s)  -- Solid oxide
    'deltaG_Zn(OH)3^-': -696.52*10**3,  # [J/mol] - Zn(OH)3^-(aq)   
    'deltaG_Zn(OH)4^2-': -860.59*10**3, # [J/mol] - Zn(OH)4^2-(aq)
    'deltaG_ZnO': -320.479*10**3,       # [J/mol] - ZnO
    'deltaG_H^+': 0,                    # [J/mol] - H^+(aq)
    'deltaG_H2O': -237.1*10**3,         # [J/mol] - H2O
    'deltaG_OH^-':-157.2*10**3,         # [J/mol] - OH^-(aq)
    'deltaG_H2': 0,                     # [J/mol] - H2(g) -- SI
    'deltaG_O2':0,                      # [J/mol] - O2(g) -- SI
}

# Entropy of formation of species at STP [REVISED POURBAIX DIAGRAM]
constants_S_formation = {
    'S_Zn': 41.63,          # [J/K*mol] - Zn(s)
    'S_Zn^2+': -109.8,      # [J/K*mol] - Zn^2+(aq)
    'S_Zn(OH)^+': -24,      # [J/K*mol] - Zn(OH)^+(aq)
    'S_Zn(OH)2': 42,        # [J/K*mol] - Zn(OH)2(s)
    'S_Zn(OH)2-eps': 77,    # [J/K*mol] - Solid form
    'S_Zn(OH)3^-': 40,      # [J/K*mol] - Zn(OH)3^-(aq)
    'S_Zn(OH)4^2-': 15,     # [J/K*mol] - Zn(OH)4^2-(aq)
    'S_ZnO': 43.65,         # [J/K*mol] - ZnO
    'S_H^+': 0,             # [J/K*mol] - H^+(aq)
    'S_H2O': 70,            # [J/K*mol] - H2O
    'S_OH^-':-11,           # [J/K*mol] - OH^-(aq)
    'S_H2': 131,            # [J/K*mol] - H2(g) -- SI
    'S_O2':205,             # [J/K*mol] - O2(g) -- SI
}

# Molar heat capacity of species at STP [REVISED POURBAIX DIAGRAM]
constants_Cp = {
    'Cp_Zn': np.array([21.334, 11.648*10**(-3), 0.054*10**(6)]),  # [J/K*mol] - Zn(s)
    'Cp_ZnO': np.array([45.338, 7.289*10**(-3), -0.573*10**(6)]), # [J/K*mol] - Zn^2+(aq)
    'Cp_Zn(OH)2-eps': np.array([74.27, 0, 0]),                    # [J/K*mol] - Zn(OH)^+(aq)
    'Cp_Zn^2+': np.array([-25.8, 0, 0]),                          # [J/K*mol] - Zn(OH)2(s)
    'Cp_Zn(OH)^+': np.array([10, 0, 0]),                          # [J/K*mol] - Solid form
    'Cp_Zn(OH)2': np.array([70, 0, 0]),                           # [J/K*mol] - Zn(OH)3^-(aq)
    'Cp_Zn(OH)3^-': np.array([94, 0, 0]),                         # [J/K*mol] - Zn(OH)4^2-(aq)
    'Cp_Zn(OH)4^2-': np.array([-284, 0, 0]),                      # [J/K*mol] - ZnO
    'Cp_H^+': np.array([0, 0, 0]),                                # [J/K*mol] - H^+(aq)
    'Cp_H2O': np.array([75, 0, 0]),                               # [J/K*mol] - H2O
    'Cp_OH^-': np.array([-149, 0, 0]),                            # [J/K*mol] - OH^-(aq)
    'Cp_H2': np.array([29, 0, 0]),                                # [J/K*mol] - H2(g) -- SI
    'Cp_O2': np.array([29, 0, 0]),                                # [J/K*mol] - O2(g) -- SI
}

######################################### THERMODYNAMIC CALCULATIONS #########################################

## Gibbs free energy for the reactions at STP using Hess' law
constants_deltaG = {
    'deltaG_I': constants_deltaG_formation['deltaG_Zn'] - constants_deltaG_formation['deltaG_Zn^2+'],                                                                                                       # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaG_II': constants_deltaG_formation['deltaG_Zn'] + constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)^+']) + constants_deltaG_formation['deltaG_H^+'],           # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaG_III-eps': constants_deltaG_formation['deltaG_Zn'] + 2*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)2-eps']+ 2*constants_deltaG_formation['deltaG_H^+']),# [J/mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + 2H2O  -- Solid hydroxide
    'deltaG_III': constants_deltaG_formation['deltaG_Zn'] + 2*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)2']+ 2*constants_deltaG_formation['deltaG_H^+']),        # [J/mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + 2H2O -- Soluble hydroxide
    'deltaG_III-ox': constants_deltaG_formation['deltaG_Zn'] + constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_ZnO']) + 2*constants_deltaG_formation['deltaG_H^+'],          # [J/mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn(s) + H2O       -- Solid oxide
    
    'deltaG_IV': constants_deltaG_formation['deltaG_Zn'] + 3*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)3^-'] + 3*constants_deltaG_formation['deltaG_H^+']),      # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
    'deltaG_V': constants_deltaG_formation['deltaG_Zn'] + 4*constants_deltaG_formation['deltaG_H2O'] - (constants_deltaG_formation['deltaG_Zn(OH)4^2-'] + 4*constants_deltaG_formation['deltaG_H^+']),      # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    'deltaG_VIII-eps': constants_deltaG_formation['deltaG_Zn(OH)2-eps'] - (constants_deltaG_formation['deltaG_Zn^2+'] + 2*constants_deltaG_formation['deltaG_OH^-']),                                       # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2  -- Solid hydroxide
    'deltaG_VIII': constants_deltaG_formation['deltaG_Zn(OH)2'] - (constants_deltaG_formation['deltaG_Zn^2+'] + 2*constants_deltaG_formation['deltaG_OH^-']),                                               # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2  -- Soluble hydroxide
    'deltaG_VIII-ox': (constants_deltaG_formation['deltaG_ZnO'] + constants_deltaG_formation['deltaG_H2O']) - (constants_deltaG_formation['deltaG_Zn^2+'] + 2*constants_deltaG_formation['deltaG_OH^-']),   # [J/mol] - Zn^2+ + 2OH^- <--> ZnO + H2O-- Oxide oxide

    'deltaG_IX-eps': constants_deltaG_formation['deltaG_Zn(OH)3^-'] - (constants_deltaG_formation['deltaG_Zn(OH)2-eps'] + constants_deltaG_formation['deltaG_OH^-']),                                       # [J/mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-   -- Solid hydroxide
    'deltaG_IX': constants_deltaG_formation['deltaG_Zn(OH)3^-'] - (constants_deltaG_formation['deltaG_Zn(OH)2'] + constants_deltaG_formation['deltaG_OH^-']),                                               # [J/mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-  -- Soluble hydroxide
    'deltaG_IX-ox': constants_deltaG_formation['deltaG_Zn(OH)3^-'] - (constants_deltaG_formation['deltaG_ZnO'] + constants_deltaG_formation['deltaG_OH^-'] + constants_deltaG_formation['deltaG_H2O']),     # [J/mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^- -- Soluble oxide

    'deltaG_X': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_Zn(OH)3^-'] + constants_deltaG_formation['deltaG_OH^-']),                                             # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaG_XI-eps': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_Zn(OH)2-eps'] + 2*constants_deltaG_formation['deltaG_OH^-']),                                    # [J/mol] - Zn(OH)2(s) + 2OH^- <--> Zn(OH)4^2-     -- Solid hydroxide
    'deltaG_XI': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_Zn(OH)2'] + 2*constants_deltaG_formation['deltaG_OH^-']),                                            # [J/mol] - Zn(OH)2(aq) + 2OH^- <--> Zn(OH)4^2-    -- Soluble hydroxide
    'deltaG_XI-ox': constants_deltaG_formation['deltaG_Zn(OH)4^2-'] - (constants_deltaG_formation['deltaG_ZnO'] + constants_deltaG_formation['deltaG_H2O']+ 2*constants_deltaG_formation['deltaG_OH^-']),   # [J/mol] - ZnO(s) + H2O + 2OH^- <--> Zn(OH)4^2-   -- Solid oxide
    
    'deltaG_HER': 2*constants_deltaG_formation['deltaG_H2'] - 4*constants_deltaG_formation['deltaG_H^+'],                                                                                                   # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaG_OER': 2*constants_deltaG_formation['deltaG_H2O'] - (4*constants_deltaG_formation['deltaG_H^+'] + constants_deltaG_formation['deltaG_O2']),                                                      # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaG_W': constants_deltaG_formation['deltaG_H^+'] + constants_deltaG_formation['deltaG_OH^-'] - constants_deltaG_formation['deltaG_H2O'],                                                            # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## Entropy of reaction at STP using Hess' law
constants_deltaS = {
    'deltaS_I': constants_S_formation['S_Zn'] - constants_S_formation['S_Zn^2+'],                                                                                   # [J/K*mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaS_II': constants_S_formation['S_Zn'] + constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)^+']) + constants_S_formation['S_H^+'],           # [J/K*mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaS_III-eps': constants_S_formation['S_Zn'] + 2*constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)2-eps']+ 2*constants_S_formation['S_H^+']),# [J/K*mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + 2H2O   -- Solid hydroxide
    'deltaS_III': constants_S_formation['S_Zn'] + 2*constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)2']+ 2*constants_S_formation['S_H^+']),        # [J/K*mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + 2H2O  -- Soluble hydroxide
    'deltaS_III-ox': constants_S_formation['S_Zn'] + constants_S_formation['S_H2O'] - (constants_S_formation['S_ZnO']+ 2*constants_S_formation['S_H^+']),           # [J/K*mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn(s) + H2O        -- Solid oxide
    
    'deltaS_IV': constants_S_formation['S_Zn'] + 3*constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)3^-'] + 3*constants_S_formation['S_H^+']),      # [J/K*mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
    'deltaS_V': constants_S_formation['S_Zn'] + 4*constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn(OH)4^2-'] + 4*constants_S_formation['S_H^+']),      # [J/K*mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    'deltaS_VIII-eps': constants_S_formation['S_Zn(OH)2-eps'] - (constants_S_formation['S_Zn^2+'] + 2*constants_S_formation['S_OH^-']),                             # [J/K*mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(s)    -- Solid hydroxide
    'deltaS_VIII': constants_S_formation['S_Zn(OH)2'] - (constants_S_formation['S_Zn^2+'] + 2*constants_S_formation['S_OH^-']),                                     # [J/K*mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(aq)   -- Soluble hydroxide
    'deltaS_VIII-ox': constants_S_formation['S_ZnO'] + constants_S_formation['S_H2O'] - (constants_S_formation['S_Zn^2+'] + 2*constants_S_formation['S_OH^-']),     # [J/K*mol] - Zn^2+ + 2OH^- <--> ZnO(s) + H2O  -- Solid oxide
    
    'deltaS_IX-eps': constants_S_formation['S_Zn(OH)3^-'] - (constants_S_formation['S_Zn(OH)2-eps'] + constants_S_formation['S_OH^-']),                             # [J/K*mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-      -- Solid hydroxide
    'deltaS_IX': constants_S_formation['S_Zn(OH)3^-'] - (constants_S_formation['S_Zn(OH)2'] + constants_S_formation['S_OH^-']),                                     # [J/K*mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-     -- Soluble hydroxide
    'deltaS_IX-ox': constants_S_formation['S_Zn(OH)3^-'] - (constants_S_formation['S_ZnO'] + constants_S_formation['S_H2O'] + constants_S_formation['S_OH^-']),     # [J/K*mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^-    -- Solid oxide
    
    'deltaS_X': constants_S_formation['S_Zn(OH)4^2-'] - (constants_S_formation['S_Zn(OH)3^-'] + constants_S_formation['S_OH^-']),                                   # [J/K*mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaS_XI-eps': constants_S_formation['S_Zn(OH)4^2-'] - (constants_S_formation['S_Zn(OH)2-eps'] + 2*constants_S_formation['S_OH^-']),                          # [J/K*mol] - Zn(OH)2(s) +2OH^- <--> Zn(OH)4^2-     -- Solid hydroxide
    'deltaS_XI': constants_S_formation['S_Zn(OH)4^2-'] - (constants_S_formation['S_Zn(OH)2'] + 2*constants_S_formation['S_OH^-']),                                  # [J/K*mol] - Zn(OH)2(aq) +2OH^- <--> Zn(OH)4^2-    -- Soluble hydroxide
    'deltaS_XI-ox': constants_S_formation['S_Zn(OH)4^2-'] - (constants_S_formation['S_ZnO'] + constants_S_formation['S_H2O'] + 2*constants_S_formation['S_OH^-']),  # [J/K*mol] - ZnO(s) + H2O +2OH^- <--> Zn(OH)4^2-   -- Solid oxide
    
    'deltaS_HER': 2*constants_S_formation['S_H2'] - 4*constants_S_formation['S_H^+'],                                                                               # [J/K*mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaS_OER': 2*constants_S_formation['S_H2O'] - (4*constants_S_formation['S_H^+'] + constants_S_formation['S_O2']),                                            # [J/K*mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaS_W': constants_S_formation['S_H^+'] + constants_S_formation['S_OH^-'] - constants_S_formation['S_H2O'],                                                  # [J/K*mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## Enthalpy for reaction as STP using Gibbs free energy and entropy for the reactions
constants_deltaH = {
    'deltaH_I': constants_deltaG['deltaG_I'] + constants['T']*constants_deltaS['deltaS_I'],                     # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaH_II': constants_deltaG['deltaG_II'] + constants['T']*constants_deltaS['deltaS_II'],                  # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaH_III-eps': constants_deltaG['deltaG_III-eps'] + constants['T']*constants_deltaS['deltaS_III-eps'],   # [J/mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaH_III': constants_deltaG['deltaG_III'] + constants['T']*constants_deltaS['deltaS_III'],               # [J/mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaH_III-ox': constants_deltaG['deltaG_III-ox'] + constants['T']*constants_deltaS['deltaS_III-ox'],      # [J/mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
    
    'deltaH_IV': constants_deltaG['deltaG_IV'] + constants['T']*constants_deltaS['deltaS_IV'],                  # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
    'deltaH_V': constants_deltaG['deltaG_V'] + constants['T']*constants_deltaS['deltaS_V'],                     # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    'deltaH_VIII-eps': constants_deltaG['deltaG_VIII-eps'] + constants['T']*constants_deltaS['deltaS_VIII-eps'],# [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(s)
    'deltaH_VIII': constants_deltaG['deltaG_VIII'] + constants['T']*constants_deltaS['deltaS_VIII'],            # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(aq)
    'deltaH_VIII-ox': constants_deltaG['deltaG_VIII-ox'] + constants['T']*constants_deltaS['deltaS_VIII-ox'],   # [J/mol] - Zn^2+ + 2OH^- <--> ZnO(s) + H2O
    
    'deltaH_IX-eps': constants_deltaG['deltaG_IX-eps'] + constants['T']*constants_deltaS['deltaS_IX-eps'],      # [J/mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-
    'deltaH_IX': constants_deltaG['deltaG_IX'] + constants['T']*constants_deltaS['deltaS_IX'],                  # [J/mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-
    'deltaH_IX-ox': constants_deltaG['deltaG_IX-ox'] + constants['T']*constants_deltaS['deltaS_IX-ox'],         # [J/mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^-
    
    'deltaH_X': constants_deltaG['deltaG_X'] + constants['T']*constants_deltaS['deltaS_X'],                     # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaH_XI-eps': constants_deltaG['deltaG_XI-eps'] + constants['T']*constants_deltaS['deltaS_XI-eps'],      # [J/mol] - Zn(OH)2(s) + 2OH^- <--> Zn(OH)4^2-
    'deltaH_XI': constants_deltaG['deltaG_XI'] + constants['T']*constants_deltaS['deltaS_XI'],                  # [J/mol] - Zn(OH)2(aq) + 2OH^- <--> Zn(OH)4^2-
    'deltaH_XI-ox': constants_deltaG['deltaG_XI-ox'] + constants['T']*constants_deltaS['deltaS_XI-ox'],         # [J/mol] - ZnO(s) + H2O + 2OH^- <--> Zn(OH)4^2-
    
    'deltaH_HER': constants_deltaG['deltaG_HER'] + constants['T']*constants_deltaS['deltaS_HER'],               # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaH_OER': constants_deltaG['deltaG_OER'] + constants['T']*constants_deltaS['deltaS_OER'],               # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaH_W': constants_deltaG['deltaG_W'] + constants['T']*constants_deltaS['deltaS_W'],                     # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## Reaction molar heat capacity for the reactions
constants_deltaCp = {
    'deltaCp_I': constants_Cp['Cp_Zn'] - constants_Cp['Cp_Zn^2+'],                                                                      # [J/K*mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaCp_II': constants_Cp['Cp_Zn'] + constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)^+'] + constants_Cp['Cp_H^+']),              # [J/K*mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaCp_III-eps': constants_Cp['Cp_Zn'] + 2*constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)2-eps']+ 2*constants_Cp['Cp_H^+']),   # [J/K*mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaCp_III': constants_Cp['Cp_Zn'] + 2*constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)2']+ 2*constants_Cp['Cp_H^+']),           # [J/K*mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaCp_III-ox': constants_Cp['Cp_Zn'] + constants_Cp['Cp_H2O'] - (constants_Cp['Cp_ZnO'] + 2*constants_Cp['Cp_H^+']),             # [J/K*mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn + H2O

    'deltaCp_IV': constants_Cp['Cp_Zn'] + 3*constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)3^-'] + 3*constants_Cp['Cp_H^+']),         # [J/K*mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
    'deltaCp_V': constants_Cp['Cp_Zn'] + 4*constants_Cp['Cp_H2O'] - (constants_Cp['Cp_Zn(OH)4^2-'] + 4*constants_Cp['Cp_H^+']),         # [J/K*mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    
    'deltaCp_VIII-eps': constants_Cp['Cp_Zn(OH)2-eps'] - (constants_Cp['Cp_Zn^2+'] + 2*constants_Cp['Cp_OH^-']),                        # [J/K*mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(s)
    'deltaCp_VIII': constants_Cp['Cp_Zn(OH)2'] - (constants_Cp['Cp_Zn^2+'] + 2*constants_Cp['Cp_OH^-']),                                # [J/K*mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(aq)
    'deltaCp_VIII-ox': constants_Cp['Cp_ZnO'] + constants_Cp['Cp_H2O']- (constants_Cp['Cp_Zn^2+'] + 2*constants_Cp['Cp_OH^-']),         # [J/K*mol] - Zn^2+ + 2OH^- <--> ZnO(s) + H2O
    
    'deltaCp_IX-eps': constants_Cp['Cp_Zn(OH)3^-'] - (constants_Cp['Cp_Zn(OH)2-eps'] + constants_Cp['Cp_OH^-']),                        # [J/K*mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-
    'deltaCp_IX': constants_Cp['Cp_Zn(OH)3^-'] - (constants_Cp['Cp_Zn(OH)2'] + constants_Cp['Cp_OH^-']),                                # [J/K*mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-
    'deltaCp_IX-ox': constants_Cp['Cp_Zn(OH)3^-'] - (constants_Cp['Cp_Zn(OH)2'] + constants_Cp['Cp_H2O'] + constants_Cp['Cp_OH^-']),    # [J/K*mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^-
    
    'deltaCp_X': constants_Cp['Cp_Zn(OH)4^2-'] - (constants_Cp['Cp_Zn(OH)3^-'] + constants_Cp['Cp_OH^-']),                              # [J/K*mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaCp_XI-eps': constants_Cp['Cp_Zn(OH)4^2-'] - (constants_Cp['Cp_Zn(OH)2-eps'] + 2*constants_Cp['Cp_OH^-']),                     # [J/K*mol] - Zn(OH)2(s) + 2OH^- <--> Zn(OH)4^2-  
    'deltaCp_XI': constants_Cp['Cp_Zn(OH)4^2-'] - (constants_Cp['Cp_Zn(OH)2'] + 2*constants_Cp['Cp_OH^-']),                             # [J/K*mol] - Zn(OH)2(aq) + 2OH^- <--> Zn(OH)4^2- 
    'deltaCp_XI-ox': constants_Cp['Cp_Zn(OH)4^2-'] - (constants_Cp['Cp_ZnO'] + constants_Cp['Cp_H2O'] + 2*constants_Cp['Cp_OH^-']),     # [J/K*mol] - ZnO(s) + H2O + 2OH^- <--> Zn(OH)4^2- 
    
    'deltaCp_HER': 2*constants_Cp['Cp_H2'] - 4*constants_Cp['Cp_H^+'],                                                                  # [J/K*mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaCp_OER': 2*constants_Cp['Cp_H2O'] - (4*constants_Cp['Cp_H^+'] + constants_Cp['Cp_O2']),                                       # [J/K*mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O  
    'deltaCp_W': constants_Cp['Cp_H^+'] + constants_Cp['Cp_OH^-'] - constants_Cp['Cp_H2O'],                                             # [J/K*mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## E0 at standard state
constants_E0 = {
    'EI_0': -(constants_deltaG['deltaG_I'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn^2+(aq) + 2e^- --> Zn(s)
    'EII_0': -(constants_deltaG['deltaG_II'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)^+ --> Zn(s) - Not used
    
    'EIII_0-eps': -(constants_deltaG['deltaG_III-eps'])/(constants['n']*constants['F']),# [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  solid form
    'EIII_0': -(constants_deltaG['deltaG_III'])/(constants['n']*constants['F']),        # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  soluble form
    'EIII_0-ox': -(constants_deltaG['deltaG_III-ox'])/(constants['n']*constants['F']),  # [V vs SHE] - ZnO --> Zn(s) - Not used
    
    'EIV_0': -(constants_deltaG['deltaG_IV'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)3^- --> Zn(s)
    'EV_0': -(constants_deltaG['deltaG_V'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
    'EHER_0': -(constants_deltaG['deltaG_HER'])/(2*constants['n']*constants['F']),      # [V vs SHE] - HER
    'EOER_0': -(constants_deltaG['deltaG_OER'])/(2*constants['n']*constants['F']),      # [V vs SHE] - OER
}

################################ THERMODYNAMIC DATA AT DIFFERENT TEMPERATURES ################################

## Gibbs free energy using no assumptions, except that the heat capacities are valid in this range
deltaG_new_T = {
    'deltaG_I': deltaG_T2(constants_deltaH['deltaH_I'], constants_deltaS['deltaS_I'], constants_deltaCp['deltaCp_I'], constants['T'], T),                               # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaG_II': deltaG_T2(constants_deltaH['deltaH_II'], constants_deltaS['deltaS_II'], constants_deltaCp['deltaCp_II'], constants['T'], T),                           # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaG_III-eps': deltaG_T2(constants_deltaH['deltaH_III-eps'], constants_deltaS['deltaS_III-eps'], constants_deltaCp['deltaCp_III-eps'], constants['T'], T),       # [J/mol] - Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaG_III': deltaG_T2(constants_deltaH['deltaH_III'], constants_deltaS['deltaS_III'], constants_deltaCp['deltaCp_III'], constants['T'], T),                       # [J/mol] - Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaG_III-ox': deltaG_T2(constants_deltaH['deltaH_III-ox'], constants_deltaS['deltaS_III-ox'], constants_deltaCp['deltaCp_III-ox'], constants['T'], T),           # [J/mol] - ZnO + 2e^- + 2H^+ <--> Zn + H2O
    
    'deltaG_IV': deltaG_T2(constants_deltaH['deltaH_IV'], constants_deltaS['deltaS_IV'], constants_deltaCp['deltaCp_IV'], constants['T'], T),                           # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O 
    'deltaG_V': deltaG_T2(constants_deltaH['deltaH_V'], constants_deltaS['deltaS_V'], constants_deltaCp['deltaCp_V'], constants['T'], T),                               # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    'deltaG_VIII-eps': deltaG_T2(constants_deltaH['deltaH_VIII-eps'], constants_deltaS['deltaS_VIII-eps'], constants_deltaCp['deltaCp_VIII-eps'], constants['T'], T),   # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
    'deltaG_VIII': deltaG_T2(constants_deltaH['deltaH_VIII'], constants_deltaS['deltaS_VIII'], constants_deltaCp['deltaCp_VIII'], constants['T'], T),                   # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
    'deltaG_VIII-ox': deltaG_T2(constants_deltaH['deltaH_VIII-ox'], constants_deltaS['deltaS_VIII-ox'], constants_deltaCp['deltaCp_VIII-ox'], constants['T'], T),       # [J/mol] - Zn^2+ + 2OH^- <--> ZnO + H2O
    
    'deltaG_IX-eps': deltaG_T2(constants_deltaH['deltaH_IX-eps'], constants_deltaS['deltaS_IX-eps'], constants_deltaCp['deltaCp_IX-eps'], constants['T'], T),           # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
    'deltaG_IX': deltaG_T2(constants_deltaH['deltaH_IX'], constants_deltaS['deltaS_IX'], constants_deltaCp['deltaCp_IX'], constants['T'], T),                           # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
    'deltaG_IX-ox': deltaG_T2(constants_deltaH['deltaH_IX-ox'], constants_deltaS['deltaS_IX-ox'], constants_deltaCp['deltaCp_IX-ox'], constants['T'], T),                  # [J/mol] - ZnO + H2O + OH^- <--> Zn(OH)3^-
    
    'deltaG_X': deltaG_T2(constants_deltaH['deltaH_X'], constants_deltaS['deltaS_X'], constants_deltaCp['deltaCp_X'], constants['T'], T),                               # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaG_XI-eps': deltaG_T2(constants_deltaH['deltaH_XI-eps'], constants_deltaS['deltaS_XI-eps'], constants_deltaCp['deltaCp_XI-eps'], constants['T'], T),           # [J/mol] - Zn(OH)2 + 2OH^- <--> Zn(OH)4^2-
    'deltaG_XI': deltaG_T2(constants_deltaH['deltaH_XI'], constants_deltaS['deltaS_XI'], constants_deltaCp['deltaCp_XI'], constants['T'], T),                           # [J/mol] - Zn(OH)2 + 2OH^- <--> Zn(OH)4^2-
    'deltaG_XI-ox': deltaG_T2(constants_deltaH['deltaH_XI-ox'], constants_deltaS['deltaS_XI-ox'], constants_deltaCp['deltaCp_XI-ox'], constants['T'], T),               # [J/mol] - ZnO + H2O + 2OH^- <--> Zn(OH)4^2-
    
    'deltaG_HER': deltaG_T2(constants_deltaH['deltaH_HER'], constants_deltaS['deltaS_HER'], constants_deltaCp['deltaCp_HER'], constants['T'], T),                       # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaG_OER': deltaG_T2(constants_deltaH['deltaH_OER'], constants_deltaS['deltaS_OER'], constants_deltaCp['deltaCp_OER'], constants['T'], T),                       # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaG_W': deltaG_T2(constants_deltaH['deltaH_W'], constants_deltaS['deltaS_W'], constants_deltaCp['deltaCp_W'], constants['T'], T),                               # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## Gibbs free energy using the van't Hoff approximation
deltaG_Vant_Hoff = {
    'deltaG_I': vant_Hoff(constants_deltaG['deltaG_I'], constants_deltaH['deltaH_I'], constants['T'], T),                       # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaG_II': vant_Hoff(constants_deltaG['deltaG_II'], constants_deltaH['deltaH_II'], constants['T'], T),                    # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaG_III-eps': vant_Hoff(constants_deltaG['deltaG_III-eps'], constants_deltaH['deltaH_III-eps'], constants['T'], T),     # [J/mol] - Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaG_III': vant_Hoff(constants_deltaG['deltaG_III'], constants_deltaH['deltaH_III'], constants['T'], T),                 # [J/mol] - Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaG_III-ox': vant_Hoff(constants_deltaG['deltaG_III-ox'], constants_deltaH['deltaH_III-ox'], constants['T'], T),        # [J/mol] - ZnO + 2e^- + 2H^+ <--> Zn + H2O
    
    'deltaG_IV': vant_Hoff(constants_deltaG['deltaG_IV'], constants_deltaH['deltaH_IV'], constants['T'], T),                    # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
    'deltaG_V': vant_Hoff(constants_deltaG['deltaG_V'], constants_deltaH['deltaH_V'], constants['T'], T),                       # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    'deltaG_VIII-eps': vant_Hoff(constants_deltaG['deltaG_VIII-eps'], constants_deltaH['deltaH_VIII-eps'], constants['T'], T),  # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
    'deltaG_VIII': vant_Hoff(constants_deltaG['deltaG_VIII'], constants_deltaH['deltaH_VIII'], constants['T'], T),              # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
    'deltaG_VIII-ox': vant_Hoff(constants_deltaG['deltaG_VIII-ox'], constants_deltaH['deltaH_VIII-ox'], constants['T'], T),     # [J/mol] - Zn^2+ + 2OH^- <--> ZnO + H2O
    
    'deltaG_IX-eps': vant_Hoff(constants_deltaG['deltaG_IX-eps'], constants_deltaH['deltaH_IX-eps'], constants['T'], T),        # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
    'deltaG_IX': vant_Hoff(constants_deltaG['deltaG_IX'], constants_deltaH['deltaH_IX'], constants['T'], T),                    # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
    'deltaG_IX-ox': vant_Hoff(constants_deltaG['deltaG_IX-ox'], constants_deltaH['deltaH_IX-ox'], constants['T'], T),           # [J/mol] - ZnO + H2O + OH^- <--> Zn(OH)3^-
    
    'deltaG_X': vant_Hoff(constants_deltaG['deltaG_X'], constants_deltaH['deltaH_X'], constants['T'], T),                       # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaG_XI-eps': vant_Hoff(constants_deltaG['deltaG_XI-eps'], constants_deltaH['deltaH_XI-eps'], constants['T'], T),        # [J/mol] - Zn(OH)2 + 2OH^- <--> Zn(OH)4^2-
    'deltaG_XI': vant_Hoff(constants_deltaG['deltaG_XI'], constants_deltaH['deltaH_XI'], constants['T'], T),                    # [J/mol] - Zn(OH)2 + 2OH^- <--> Zn(OH)4^2-
    'deltaG_XI-ox': vant_Hoff(constants_deltaG['deltaG_XI-ox'], constants_deltaH['deltaH_XI-ox'], constants['T'], T),           # [J/mol] - ZnO + H2O + 2OH^- <--> Zn(OH)4^2-
    
    'deltaG_HER': vant_Hoff(constants_deltaG['deltaG_HER'], constants_deltaH['deltaH_HER'], constants['T'], T),                 # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaG_OER': vant_Hoff(constants_deltaG['deltaG_OER'], constants_deltaH['deltaH_OER'], constants['T'], T),                 # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaG_W': vant_Hoff(constants_deltaG['deltaG_W'], constants_deltaH['deltaH_W'], constants['T'], T),                       # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## DeltaG using enthalpy and entropy as weak functions of temperature
deltaG_approx = {
    'deltaG_I': deltaG_weak(constants_deltaH['deltaH_I'], constants_deltaS['deltaS_I'], T),                     # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaG_II': deltaG_weak(constants_deltaH['deltaH_II'], constants_deltaS['deltaS_II'], T),                  # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaG_III-eps': deltaG_weak(constants_deltaH['deltaH_III-eps'], constants_deltaS['deltaS_III-eps'], T),   # [J/mol] - Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaG_III': deltaG_weak(constants_deltaH['deltaH_III'], constants_deltaS['deltaS_III'], T),               # [J/mol] - Zn(OH)2 + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaG_III-ox': deltaG_weak(constants_deltaH['deltaH_III-ox'], constants_deltaS['deltaS_III-ox'], T),      # [J/mol] - ZnO + 2e^- + 2H^+ <--> Zn + H2O
    
    'deltaG_IV': deltaG_weak(constants_deltaH['deltaH_IV'], constants_deltaS['deltaS_IV'], T),                  # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O 
    'deltaG_V': deltaG_weak(constants_deltaH['deltaH_V'], constants_deltaS['deltaS_V'], T),                     # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    'deltaG_VIII-eps': deltaG_weak(constants_deltaH['deltaH_VIII-eps'], constants_deltaS['deltaS_VIII-eps'], T),# [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
    'deltaG_VIII': deltaG_weak(constants_deltaH['deltaH_VIII'], constants_deltaS['deltaS_VIII'], T),            # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
    'deltaG_VIII-ox': deltaG_weak(constants_deltaH['deltaH_VIII-ox'], constants_deltaS['deltaS_VIII-ox'], T),   # [J/mol] - Zn^2+ + 2OH^- <--> ZnO + H2O
    
    'deltaG_IX-eps': deltaG_weak(constants_deltaH['deltaH_IX-eps'], constants_deltaS['deltaS_IX-eps'], T),      # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
    'deltaG_IX': deltaG_weak(constants_deltaH['deltaH_IX'], constants_deltaS['deltaS_IX'], T),                  # [J/mol] - Zn(OH)2 + OH^- <--> Zn(OH)3^-
    'deltaG_IX-ox': deltaG_weak(constants_deltaH['deltaH_IX-ox'], constants_deltaS['deltaS_IX-ox'], T),         # [J/mol] - ZnO + H2O + OH^- <--> Zn(OH)3^-
    
    'deltaG_X': deltaG_weak(constants_deltaH['deltaH_X'], constants_deltaS['deltaS_X'], T),                     # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaG_XI-eps': deltaG_weak(constants_deltaH['deltaH_XI-eps'], constants_deltaS['deltaS_XI-eps'], T),      # [J/mol] - Zn(OH)2 + 2OH^- <--> Zn(OH)4^2-
    'deltaG_XI': deltaG_weak(constants_deltaH['deltaH_XI'], constants_deltaS['deltaS_XI'], T),                  # [J/mol] - Zn(OH)2 + 2OH^- <--> Zn(OH)4^2-
    'deltaG_XI-ox': deltaG_weak(constants_deltaH['deltaH_XI-ox'], constants_deltaS['deltaS_XI-ox'], T),         # [J/mol] - ZnO + H2O + 2OH^- <--> Zn(OH)4^2-
    
    'deltaG_HER': deltaG_weak(constants_deltaH['deltaH_HER'], constants_deltaS['deltaS_HER'], T),               # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaG_OER': deltaG_weak(constants_deltaH['deltaH_OER'], constants_deltaS['deltaS_OER'], T),               # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaG_W': deltaG_weak(constants_deltaH['deltaH_W'], constants_deltaS['deltaS_W'], T),                     # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## E0 using van't approximation
E0_vant_Hoff = {
    'EI_0': -(deltaG_Vant_Hoff['deltaG_I'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn^2+(aq) --> Zn(s)
    'EII_0': -(deltaG_Vant_Hoff['deltaG_II'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)^+ --> Zn(s) - Not used
    
    'EIII_0-eps': -(deltaG_Vant_Hoff['deltaG_III-eps'])/(constants['n']*constants['F']),# [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0': -(deltaG_Vant_Hoff['deltaG_III'])/(constants['n']*constants['F']),        # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0-ox': -(deltaG_Vant_Hoff['deltaG_III-ox'])/(constants['n']*constants['F']),  # [V vs SHE] - ZnO --> Zn(s) - Not used
    
    'EIV_0': -(deltaG_Vant_Hoff['deltaG_IV'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)3^- --> Zn(s)
    'EV_0': -(deltaG_Vant_Hoff['deltaG_V'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
    
    'EHER_0':-(deltaG_Vant_Hoff['deltaG_HER'])/(2*constants['n']*constants['F']),       # [V vs SHE] - HER
    'EOER_0': -(deltaG_Vant_Hoff['deltaG_OER'])/(2*constants['n']*constants['F']),      # [V vs SHE] - OER
}

## E0 using assuming that enthalpy and entropy are weak functions of  temperature
E0_approx = {
    'EI_0': -(deltaG_approx['deltaG_I'])/(constants['n']*constants['F']),               # [V vs SHE] - Zn^2+(aq) --> Zn(s)
    'EII_0': -(deltaG_approx['deltaG_II'])/(constants['n']*constants['F']),             # [V vs SHE] - Zn(OH)^+ --> Zn(s) - Not used
    
    'EIII_0-eps': -(deltaG_approx['deltaG_III-eps'])/(constants['n']*constants['F']),   # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0': -(deltaG_approx['deltaG_III'])/(constants['n']*constants['F']),           # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0-ox': -(deltaG_approx['deltaG_III-ox'])/(constants['n']*constants['F']),      # [V vs SHE] - ZnO --> Zn(s) - Not used
    
    'EIV_0': -(deltaG_approx['deltaG_IV'])/(constants['n']*constants['F']),             # [V vs SHE] - Zn(OH)3^- --> Zn(s)
    'EV_0': -(deltaG_approx['deltaG_V'])/(constants['n']*constants['F']),               # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
    
    'EHER_0':-(deltaG_approx['deltaG_HER'])/(2*constants['n']*constants['F']),          # [V vs SHE] - HER
    'EOER_0': -(deltaG_approx['deltaG_OER'])/(2*constants['n']*constants['F']),         # [V vs SHE] - OER
}

## E0 with no assumptions, except that the heat capacities are valid for this range
E0_new = {
    'EI_0': -(deltaG_new_T['deltaG_I'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn^2+(aq) --> Zn(s)
    'EII_0': -(deltaG_new_T['deltaG_II'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)^+ --> Zn(s) - Not used
    
    'EIII_0-eps': -(deltaG_new_T['deltaG_III-eps'])/(constants['n']*constants['F']),# [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0': -(deltaG_new_T['deltaG_III'])/(constants['n']*constants['F']),        # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0-ox': -(deltaG_new_T['deltaG_III-ox'])/(constants['n']*constants['F']),  # [V vs SHE] - ZnO --> Zn(s) - Not used
    
    'EIV_0': -(deltaG_new_T['deltaG_IV'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)3^- --> Zn(s)
    'EV_0': -(deltaG_new_T['deltaG_V'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
    
    'EHER_0':-(deltaG_new_T['deltaG_HER'])/(2*constants['n']*constants['F']),       # [V vs SHE] - HER
    'EOER_0': -(deltaG_new_T['deltaG_OER'])/(2*constants['n']*constants['F'])       # [V vs SHE] - OER
}

## Using that the derivative of the standard reduction potential with respect to temperature is the change in entropy 
E0_S = {
    'EI_0': E0_2(constants_E0['EI_0'], constants_deltaS['deltaS_I'], constants['T'], T, constants['n']),                    # [V vs SHE] - Zn^2+(aq) --> Zn(s)
    'EII_0': E0_2(constants_E0['EII_0'], constants_deltaS['deltaS_II'], constants['T'], T, constants['n']),                 # [V vs SHE] - Zn(OH)^+ --> Zn(s) - Not used
    
    'EIII_0-eps': E0_2(constants_E0['EIII_0-eps'], constants_deltaS['deltaS_III-eps'], constants['T'], T, constants['n']),  # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0': E0_2(constants_E0['EIII_0'], constants_deltaS['deltaS_III'], constants['T'], T, constants['n']),              # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
    'EIII_0-ox': E0_2(constants_E0['EIII_0-ox'], constants_deltaS['deltaS_III-ox'], constants['T'], T, constants['n']),     # [V vs SHE] - ZnO --> Zn(s) - Not used
    
    'EIV_0': E0_2(constants_E0['EIV_0'], constants_deltaS['deltaS_IV'], constants['T'], T, constants['n']),                 # [V vs SHE] - Zn(OH)3^- --> Zn(s)
    'EV_0': E0_2(constants_E0['EV_0'], constants_deltaS['deltaS_V'], constants['T'], T, constants['n']),                    # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
    
    'EHER_0': E0_2(constants_E0['EHER_0'], constants_deltaS['deltaS_HER'], constants['T'], T, 2*constants['n']),            # [V vs SHE] - HER
    'EOER_0': E0_2(constants_E0['EOER_0'], constants_deltaS['deltaS_OER'], constants['T'], T, 2*constants['n'])             # [V vs SHE] - OER
}

# The pZn where the equilibrium Zn(OH)2(s) <--> Zn(OH)3^-1 and Zn(OH)3^-1 <--> Zn(OH)4^2- becomes the same and the domain for Zn(OH)3^- vanishes
pZn_threshold_eps = (deltaG_new_T['deltaG_IX-eps']-deltaG_new_T['deltaG_X'])/(constants['R']*T*np.log(10))
# The pZn where the equilibrium Zn(OH)2(aq) <--> Zn(OH)3^-1 and Zn(OH)3^-1 <--> Zn(OH)4^2- becomes the same and the domain for Zn(OH)3^- vanishes
pZn_threshold = (deltaG_new_T['deltaG_IX']-deltaG_new_T['deltaG_X'])/(constants['R']*T*np.log(10))
# The pZn where the equilibrium ZnO <--> Zn(OH)3^-1 and Zn(OH)3^-1 <--> Zn(OH)4^2- becomes the same and the domain for Zn(OH)3^- vanishes
pZn_threshold_ox = (deltaG_new_T['deltaG_IX-ox'] - deltaG_new_T['deltaG_X'])/(constants['R']*T*np.log(10))
# Threshold for stability of ZnO(s)/Zn(OH)2(aq)
lna = -(constants_deltaG_formation['deltaG_Zn(OH)2'] - (constants_deltaG_formation['deltaG_H2O'] + constants_deltaG_formation['deltaG_ZnO']))/(constants['R']*constants['T'])
a = np.exp(lna)
pa = -np.log10(a)
print(pa)

############################################## PLOTTING ##############################################

pH = np.arange(0, 16, 0.01)     # pH range

# Defining the lines for the HER and OER
EHER = E0_new['EHER_0'] - 2*constants['R']*T/(constants['n']*constants['F'])*np.log(10)*pH   # HER
EOER = E0_new['EOER_0'] - 2*constants['R']*T/(constants['n']*constants['F'])*np.log(10)*pH   # OER

# Neutral pH
Kw = np.exp(-deltaG_new_T['deltaG_W']/(constants['R']*T))
pKw = -np.log10(Kw)
pH_neutral = (1/2)*pKw

# Concentration of Zn ions
pZn_values = np.array([1, 4])           # Values we iterate through for the activity of dissolved Zn - pZn = -log(c_Zn) --- Maximum value is 6 and minimum value is 0

pZn_values_dict = {'pZn': pZn_values}   # Dictionary with values
linestyle = ['-', '--']

## Lists with the pH of the different equilibriums. Using a list because we use 2 values for pZn
x_Zn2_ZnOH2 = []

# Using soluble Zn(OH)2
x_ZnOH2_ZnOH4_2 = []
x_ZnOH2_ZnOH3_1 = []
x_ZnOH3_1_ZnOH4_2 = []
# Using solid Zn(OH)2
x_Zn2_ZnOH2_eps = []
x_ZnOH2_eps_ZnOH4_2 = []
x_ZnOH2_eps_ZnOH3_1 = []
# Using oxide form ZnO
x_Zn2_ZnO = []
x_ZnO_ZnOH4_2 = []
x_ZnO_ZnOH3_1 = []

# Initialising the figure before the for-loop
plt.figure()

## Printing the Pourbaix diagram
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
    EI = np.ones(len(pH))*E0_new['EI_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZn'])
    
    # Zn(OH)^+ --> Zn(s)
    EII = E0_new['EII_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH']+pH)

    # Zn(OH)2 --> Zn(s) ----- Solid hydroxide 
    EIII_eps = E0_new['EIII_0-eps'] - 2*constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(pH)

    # Zn(OH)2 --> Zn(s) ----- Soluble hydroxide
    EIII = E0_new['EIII_0'] - 2*constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(0*constants_p['pZnOH2']/2+pH)

    # ZnO --> Zn(s) ---- Solid oxide
    EIII_ox = E0_new['EIII_0-ox'] - 2*constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(pH)
    
    # Zn(OH)3^- --> Zn(s)
    EIV = E0_new['EIV_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH3']+3*pH)

    # Zn(OH)4^2- --> Zn(s)
    EV = E0_new['EV_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH4']+4*pH)


    ## pH calculations

    # Zn^2+ --> Zn(OH)2 ---- Solid hydroxide
    pHVIII_eps = pKw + deltaG_new_T['deltaG_VIII-eps']/(2*constants['R']*T*np.log(10)) + constants_p['pZn']/2
    x_Zn2_ZnOH2_eps.append(pHVIII_eps)
    # Zn^2+ --> Zn(OH)2 ---- Soluble hydroxide
    pHVIII = pKw + deltaG_new_T['deltaG_VIII']/(2*constants['R']*T*np.log(10)) + constants_p['pZn']/2 - constants_p['pZnOH2']/2
    x_Zn2_ZnOH2.append(pHVIII)
    # Zn^2+ --> ZnO ---- Oxide oxide
    pHVIII_ox = pKw + deltaG_new_T['deltaG_VIII-ox']/(2*constants['R']*T*np.log(10)) + constants_p['pZn']/2
    x_Zn2_ZnO.append(pHVIII_ox)

    # Zn(OH)2 --> Zn(OH)3^- ---- Using Solid Zn(OH)2
    pHIX_eps = pKw + deltaG_new_T['deltaG_IX-eps']/(constants['R']*T*np.log(10)) - constants_p['pZnOH3']
    x_ZnOH2_eps_ZnOH3_1.append(pHIX_eps)
    # Zn(OH)2 --> Zn(OH)3^- ---- Using Soluble Zn(OH)2
    pHIX = pKw + deltaG_new_T['deltaG_IX']/(constants['R']*T*np.log(10)) + constants_p['pZnOH2'] - constants_p['pZnOH3']
    x_ZnOH2_ZnOH3_1.append(pHIX)
    # ZnO --> Zn(OH)3^- ---- Using oxide ZnO
    pHIX_ox = pKw + deltaG_new_T['deltaG_IX-ox']/(constants['R']*T*np.log(10)) - constants_p['pZnOH3']
    x_ZnO_ZnOH3_1.append(pHIX_ox)

    # Zn(OH)3^- --> Zn(OH)4^2-
    pHX = pKw + deltaG_new_T['deltaG_X']/(constants['R']*T*np.log(10)) + constants_p['pZnOH3'] - constants_p['pZnOH4']
    x_ZnOH3_1_ZnOH4_2.append(pHX)

    # Zn(OH)2 --> Zn(OH)4^2- ---- Using the solid form
    pHXI_eps = pKw + deltaG_new_T['deltaG_XI-eps']/(2*constants['R']*T*np.log(10)) - (1/2)*constants_p['pZnOH4']
    x_ZnOH2_eps_ZnOH4_2.append(pHXI_eps)
    # Zn(OH)2 --> Zn(OH)4^2- ---- Using the Soluble form
    pHXI = pKw + deltaG_new_T['deltaG_XI-eps']/(2*constants['R']*T*np.log(10)) + (1/2)*constants_p['pZnOH2']- (1/2)*constants_p['pZnOH4']
    x_ZnOH2_ZnOH4_2.append(pHXI)
    # ZnO--> Zn(OH)4^2- ---- Using the oxide form
    pHXI_ox = pKw + deltaG_new_T['deltaG_XI-ox']/(2*constants['R']*T*np.log(10)) - (1/2)*constants_p['pZnOH4']
    x_ZnO_ZnOH4_2.append(pHXI_ox)

    # Finding intersections for the different lines
    index_I_III_eps = np.where(EI <= EIII_eps)[0][-1] if np.any(EI <= EIII_eps) else None   # Intersection between  Zn^2+ --> Zn and Zn(OH)2(s) --> Zn 
    index_I_III = np.where(EI <= EIII)[0][-1] if np.any(EI <= EIII) else None               # Intersection between  Zn^2+ --> Zn and Zn(OH)2(aq) --> Zn
    index_I_III_ox = np.where(EI <= EIII_ox)[0][-1] if np.any(EI <= EIII_ox) else None      # Intersection between  Zn^2+ --> Zn and ZnO(s) --> Zn
    
    index_III_eps_IV = np.where(EIII_eps <= EIV)[0][-1] if np.any(EIII_eps <= EIV) else None# Intersection between  Zn(OH)2(s) --> Zn  and Zn(OH)3^- --> Zn(s)
    index_III_IV = np.where(EIII <= EIV)[0][-1] if np.any(EIII <= EIV) else None            # Intersection between  Zn(OH)2(aq) --> Zn  and Zn(OH)3^- --> Zn(s)
    index_III_ox_IV = np.where(EIII_ox <= EIV)[0][-1] if np.any(EIII_ox <= EIV) else None   # Intersection between  ZnO(s) --> Zn  and Zn(OH)3^- --> Zn(s)

    index_IV_V = np.where(EIV <= EV)[0][-1] if np.any(EIV <= EV) else None                  # Intersection between Zn(OH)3^- --> Zn(s) and Zn(OH)4^-2 --> Zn(s)
    
    index_III_eps_V = np.where(EIII_eps <= EV)[0][-1] if np.any(EIII_eps <=EV) else None    # Intersection between Zn(OH)2(s) --> Zn(s) and Zn(OH)4^-2 --> Zn(s)
    index_III_V = np.where(EIII <= EV)[0][-1] if np.any(EIII <=EV) else None                # Intersection between Zn(OH)2(aq) --> Zn(s) and Zn(OH)4^-2 --> Zn(s)
    index_III_ox_V = np.where(EIII_ox <= EV)[0][-1] if np.any(EIII_ox <=EV) else None       # Intersection between ZnO(s) --> Zn(s) and Zn(OH)4^-2 --> Zn(s)

    # If we have any pZn values over the threshold for Zn(OH)3^- we have to add that domain
    if any(pZn_values > pZn_threshold_ox): # Checking if any of the elements pZn values are bigger than the threshold
        # Plotting the lines of the Pourbaix diagram
        plt.plot(pH[:index_I_III_ox], EI[:index_I_III_ox], 'k', label='Zn$^{2+}$ - Zn', linestyle=linestyle[i]) # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        plt.vlines(pHVIII_ox, EI[index_I_III_ox], 1.5, 'k', label='Zn$^{2+}$ - ZnO', linestyle=linestyle[i])    # Chemical        - Equilibrium between Zn^2+ and ZnO(s)

        # Adding pZn text --> Should be adjusted later
        plt.text(pH[index_I_III_ox+10], 1, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        
        # Plotting for the case where we do have Zn(OH)3^-
        if pZn_values[i] >= pZn_threshold_ox:
            plt.vlines(pHX, EV[index_IV_V], 1.5, 'm', label='Zn(OH)$_{3}^-$ - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[i])                # Chemical        - Equilibrium between Zn(OH)3^- and Zn(OH)4^2-
            plt.vlines(pHIX_ox, EIV[index_III_ox_IV], 1.5, 'r', label='ZnO - Zn(OH)$_{3}^-$', linestyle=linestyle[i])                    # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)3^-
            plt.plot(pH[index_III_ox_IV:index_IV_V], EIV[index_III_ox_IV:index_IV_V], 'r', label='Zn(OH)$_{3}^{-}$ - Zn', linestyle=linestyle[i]) # Electrochemical - Equilibrium between Zn(OH)3^- and Zn
            plt.text(pH[index_III_ox_IV+10], 1, pZn_string, color='r', rotation = 90)  # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)3^-
            plt.text(pH[index_IV_V+10], 1, pZn_string, color='m', rotation = 90)    # Adding pZn value to equilibrium between Zn(OH)3^-1 and Zn(OH)4^2
            plt.plot(pH[index_IV_V:], EV[index_IV_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[i])                         # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn   
            plt.plot(pH[index_I_III_ox: index_III_ox_IV], EIII_ox[index_I_III_ox: index_III_ox_IV], 'c', label='ZnO - Zn')                  # Electrochemical - Equilibrium between Zn(OH)2 and Zn
        
        # Plotting for the case where we do not have Zn(OH)3^-
        elif pZn_values[i]<pZn_threshold_ox:
            plt.plot(pH[index_III_ox_V:], EV[index_III_ox_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[i])     # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn
            plt.vlines(pHXI_ox, EV[index_III_ox_V], 1.5, 'm', label='ZnO - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[i])          # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)4^2-
            plt.text(pH[index_III_ox_V+10], 1, pZn_string, color='m', rotation = 90)    # Adding pZn value to equilibrium between Zn(OH)3^-1 and Zn(OH)4^2 -- must be fixed
            plt.plot(pH[index_I_III_ox: index_III_ox_V], EIII_ox[index_I_III_ox: index_III_ox_V], 'c', label='ZnO - Zn')        # Electrochemical - Equilibrium between Zn(OH)2 and Zn

        # If we are at the end, we can add the text so it isn't added two times
        if i == len(pZn_values)-1:
            # Adding text to different domains --> Should be adjusted later
            plt.text(5, -1.25, 'Zn(s)', fontsize=12, color='black')                                 # Zn domain
            plt.text(np.min(x_Zn2_ZnO)-3, 0.45, 'Zn$^{2+}$(aq)', fontsize=12, color='black')        # Zn^2+ domain
            plt.text(np.max([np.max(x_ZnOH3_1_ZnOH4_2), np.max(x_ZnO_ZnOH4_2)])+0.25, -0.75, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain
            plt.text(np.min(x_ZnO_ZnOH3_1)+0.25, -0.4, 'Zn(OH)$^{-}_{3}$(aq)', fontsize=12, color='black', rotation=90)  # Zn(OH)3^- domain
            plt.text(np.max(x_Zn2_ZnO)+1, -0.25, 'ZnO(s)', fontsize=12, color='black', rotation=90) # ZnO domain

    # If none of the pZn values are above the threshold for Zn(OH)3^- we don't have to add it
    else:
        # Plotting the lines of the Pourbaix diagram
        plt.plot(pH[:index_I_III_ox], EI[:index_I_III_ox], 'k', label='Zn$^{2+}$ - Zn', linestyle=linestyle[i])                               # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        plt.vlines(pHVIII_ox, EI[index_I_III_ox], 1.5, 'k', label='Zn$^{2+}$ - ZnO', linestyle=linestyle[i])                           # Chemical        - Equilibrium between Zn^2+ and Zn(OH)2
        plt.plot(pH[index_I_III_ox: index_III_ox_V], EIII_ox[index_I_III_ox: index_III_ox_V], 'c', label='ZnO - Zn')# Electrochemical - Equilibrium between Zn(OH)2 and Zn
        plt.plot(pH[index_III_ox_V:], EV[index_III_ox_V:], 'm', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[i])                       # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn
        plt.vlines(pHXI_ox, EV[index_III_ox_V], 1.5, 'm', label='ZnO - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[i])                   # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)4^2-
        
        # Adding pZn text --> Should be adjusted later
        plt.text(pH[index_I_III_ox+10], 1, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        plt.text(pH[index_III_ox_V+10], 1, pZn_string, color='m', rotation = 90)   # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)4^2 -- must be fixed

        if i == len(pZn_values)-1:
            # Adding text to different domains --> Should be adjusted later
            plt.text(5, -1.25, 'Zn(s)', fontsize=12, color='black')                                 # Zn domain
            plt.text(np.min(x_Zn2_ZnO)-3, 0.45, 'Zn$^{2+}$(aq)', fontsize=12, color='black')      # Zn^2+ domain
            plt.text(np.max(x_ZnO_ZnOH4_2)+0.25, -0.75, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain
            plt.text(np.max(x_Zn2_ZnO)+1, -0.25, 'ZnO(s)', fontsize=12, color='black', rotation=0)              # Zn(OH)2 domain 

# Adding things that are present in the diagram for all cases, and just have to be added ones
plt.plot(pH,EHER, '--')     # Line for the HER
plt.plot(pH,EOER, '--')     # Line for the OER
plt.text(2, -0.32, 'HER', fontsize=12, color='black', rotation=-15)  # HER line
plt.text(2, 0.85, 'OER', fontsize=12, color='black', rotation=-15)   # OER line
plt.text(1, 1.30, 'T = ' + str(T-273.15) + '$^{o}C$',fontsize=12, color='black', bbox=dict(facecolor='white', alpha=1))
plt.vlines(pH_neutral, ymin=-1.5, ymax=1.5, colors='k', label='Neutral pH', alpha=0.40)             # Adding a line for neutral
plt.xlabel('pH  /  []')
plt.ylabel('Potential - E  /  V')
plt.xlim(xmin=0, xmax=max(pH))  # Set the x-axis range
plt.ylim(ymin=-1.5, ymax=1.5)
plt.show()
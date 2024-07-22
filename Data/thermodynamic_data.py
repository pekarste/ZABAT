import numpy as np

# Constants
constants = {
    'R': 8.31451,       # [J/(K*mol)]   - The universal gas constant
    'T': 25+273.15,     # [K]           - Standard temperature
    'n': 2,             # [-]           - Number of electrons transferred
    'F': 96485          # [C/mol]       - Faraday's constant
}

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
delta_r_G = {
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
    
    'deltaG_XII': constants_deltaG_formation['deltaG_Zn(OH)2'] - (constants_deltaG_formation['deltaG_H2O'] + constants_deltaG_formation['deltaG_ZnO']),                                                     # [J/mol] - ZnO(s) + H2O(l) <--> Zn(OH)2(aq)

    'deltaG_HER': 2*constants_deltaG_formation['deltaG_H2'] - 4*constants_deltaG_formation['deltaG_H^+'],                                                                                                   # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaG_OER': 2*constants_deltaG_formation['deltaG_H2O'] - (4*constants_deltaG_formation['deltaG_H^+'] + constants_deltaG_formation['deltaG_O2']),                                                      # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaG_W': constants_deltaG_formation['deltaG_H^+'] + constants_deltaG_formation['deltaG_OH^-'] - constants_deltaG_formation['deltaG_H2O'],                                                            # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## Entropy of reaction at STP using Hess' law
delta_r_S = {
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
    
    'deltaS_XII': constants_S_formation['S_Zn(OH)2'] - (constants_S_formation['S_H2O'] + constants_S_formation['S_ZnO']),                                           # [J/K*mol] - ZnO(s) + H2O(l) <--> Zn(OH)2(aq)

    'deltaS_HER': 2*constants_S_formation['S_H2'] - 4*constants_S_formation['S_H^+'],                                                                               # [J/K*mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaS_OER': 2*constants_S_formation['S_H2O'] - (4*constants_S_formation['S_H^+'] + constants_S_formation['S_O2']),                                            # [J/K*mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaS_W': constants_S_formation['S_H^+'] + constants_S_formation['S_OH^-'] - constants_S_formation['S_H2O'],                                                  # [J/K*mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## Enthalpy for reaction as STP using Gibbs free energy and entropy for the reactions
delta_r_H = {
    'deltaH_I': delta_r_G['deltaG_I'] + constants['T']*delta_r_S['deltaS_I'],                     # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
    'deltaH_II': delta_r_G['deltaG_II'] + constants['T']*delta_r_S['deltaS_II'],                  # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
    
    'deltaH_III-eps': delta_r_G['deltaG_III-eps'] + constants['T']*delta_r_S['deltaS_III-eps'],   # [J/mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaH_III': delta_r_G['deltaG_III'] + constants['T']*delta_r_S['deltaS_III'],               # [J/mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + H2O
    'deltaH_III-ox': delta_r_G['deltaG_III-ox'] + constants['T']*delta_r_S['deltaS_III-ox'],      # [J/mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
    
    'deltaH_IV': delta_r_G['deltaG_IV'] + constants['T']*delta_r_S['deltaS_IV'],                  # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
    'deltaH_V': delta_r_G['deltaG_V'] + constants['T']*delta_r_S['deltaS_V'],                     # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
    
    'deltaH_VIII-eps': delta_r_G['deltaG_VIII-eps'] + constants['T']*delta_r_S['deltaS_VIII-eps'],# [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(s)
    'deltaH_VIII': delta_r_G['deltaG_VIII'] + constants['T']*delta_r_S['deltaS_VIII'],            # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(aq)
    'deltaH_VIII-ox': delta_r_G['deltaG_VIII-ox'] + constants['T']*delta_r_S['deltaS_VIII-ox'],   # [J/mol] - Zn^2+ + 2OH^- <--> ZnO(s) + H2O
    
    'deltaH_IX-eps': delta_r_G['deltaG_IX-eps'] + constants['T']*delta_r_S['deltaS_IX-eps'],      # [J/mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-
    'deltaH_IX': delta_r_G['deltaG_IX'] + constants['T']*delta_r_S['deltaS_IX'],                  # [J/mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-
    'deltaH_IX-ox': delta_r_G['deltaG_IX-ox'] + constants['T']*delta_r_S['deltaS_IX-ox'],         # [J/mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^-
    
    'deltaH_X': delta_r_G['deltaG_X'] + constants['T']*delta_r_S['deltaS_X'],                     # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
    
    'deltaH_XI-eps': delta_r_G['deltaG_XI-eps'] + constants['T']*delta_r_S['deltaS_XI-eps'],      # [J/mol] - Zn(OH)2(s) + 2OH^- <--> Zn(OH)4^2-
    'deltaH_XI': delta_r_G['deltaG_XI'] + constants['T']*delta_r_S['deltaS_XI'],                  # [J/mol] - Zn(OH)2(aq) + 2OH^- <--> Zn(OH)4^2-
    'deltaH_XI-ox': delta_r_G['deltaG_XI-ox'] + constants['T']*delta_r_S['deltaS_XI-ox'],         # [J/mol] - ZnO(s) + H2O + 2OH^- <--> Zn(OH)4^2-
    
    'deltaH_XII': delta_r_G['deltaG_XII'] + constants['T']*delta_r_S['deltaS_XII'],               # [J/mol] - ZnO(s) + H2O(l) <--> Zn(OH)2(aq)

    'deltaH_HER': delta_r_G['deltaG_HER'] + constants['T']*delta_r_S['deltaS_HER'],               # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaH_OER': delta_r_G['deltaG_OER'] + constants['T']*delta_r_S['deltaS_OER'],               # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
    'deltaH_W': delta_r_G['deltaG_W'] + constants['T']*delta_r_S['deltaS_W'],                     # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## Reaction molar heat capacity for the reactions
delta_r_Cp = {
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
    
    'deltaCp_XII': constants_Cp['Cp_Zn(OH)2'] - (constants_Cp['Cp_H2O'] + constants_Cp['Cp_ZnO']),                                      # [J/K*mol] - ZnO(s) + H2O(l) <--> Zn(OH)2(aq)

    'deltaCp_HER': 2*constants_Cp['Cp_H2'] - 4*constants_Cp['Cp_H^+'],                                                                  # [J/K*mol] - 4H^+ + 4e^- <--> 2H2(g)
    'deltaCp_OER': 2*constants_Cp['Cp_H2O'] - (4*constants_Cp['Cp_H^+'] + constants_Cp['Cp_O2']),                                       # [J/K*mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O  
    'deltaCp_W': constants_Cp['Cp_H^+'] + constants_Cp['Cp_OH^-'] - constants_Cp['Cp_H2O'],                                             # [J/K*mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
}

## E0 at standard state
constants_E0 = {
    'EI_0': -(delta_r_G['deltaG_I'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn^2+(aq) + 2e^- --> Zn(s)
    'EII_0': -(delta_r_G['deltaG_II'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)^+ --> Zn(s) - Not used
    
    'EIII_0-eps': -(delta_r_G['deltaG_III-eps'])/(constants['n']*constants['F']),# [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  solid form
    'EIII_0': -(delta_r_G['deltaG_III'])/(constants['n']*constants['F']),        # [V vs SHE] - Zn(OH)2(s) --> Zn(s) -  soluble form
    'EIII_0-ox': -(delta_r_G['deltaG_III-ox'])/(constants['n']*constants['F']),  # [V vs SHE] - ZnO --> Zn(s) - Not used
    
    'EIV_0': -(delta_r_G['deltaG_IV'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)3^- --> Zn(s)
    'EV_0': -(delta_r_G['deltaG_V'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
    'EHER_0': -(delta_r_G['deltaG_HER'])/(2*constants['n']*constants['F']),      # [V vs SHE] - HER
    'EOER_0': -(delta_r_G['deltaG_OER'])/(2*constants['n']*constants['F']),      # [V vs SHE] - OER
}
'''
POURBAIX DIAGRAM OF Zn

This script is generating a Pourbaix diagram of Zn for the ZABAT project.
It is using thermodynamic data listed by Beverskog et al [1] and SI Chemical Data [2].
Furthermore, it is also temperature dependent and meant to generate a Pourbaix diagram 
for the temperature range 25-100 degrees.

Made by: Pål Emil England Karstensen and Sidsel Meli Hanetho
'''

import matplotlib.pyplot as plt
import numpy as np
from Functions.Functions import vant_Hoff, deltaG_weak, deltaG_T2, E0_2, add_polygon, E_OER_HER
from Data.thermodynamic_data import *   # Imports all the thermodynamic data as well as the constants
import streamlit as st


def PourbaixDiagram(pZn, temperature):

    T = temperature          # [K] - Temperature
    # Concentration of Zn ions
    pZn_value = pZn     # Activity of dissolved Zn - pZn = -log(c_Zn) --- Maximum value is 6 and minimum value is 0

    ################################ THERMODYNAMIC DATA AT DIFFERENT TEMPERATURES ################################

    ## Gibbs free energy using no assumptions, except that the heat capacities are valid in this range
    deltaG_new_T = {
        'deltaG_I': deltaG_T2(delta_r_H['deltaH_I'], delta_r_S['deltaS_I'], delta_r_Cp['deltaCp_I'], constants['T'], T),                               # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
        'deltaG_II': deltaG_T2(delta_r_H['deltaH_II'], delta_r_S['deltaS_II'], delta_r_Cp['deltaCp_II'], constants['T'], T),                           # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
        
        'deltaG_III-eps': deltaG_T2(delta_r_H['deltaH_III-eps'], delta_r_S['deltaS_III-eps'], delta_r_Cp['deltaCp_III-eps'], constants['T'], T),       # [J/mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
        'deltaG_III': deltaG_T2(delta_r_H['deltaH_III'], delta_r_S['deltaS_III'], delta_r_Cp['deltaCp_III'], constants['T'], T),                       # [J/mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + H2O
        'deltaG_III-ox': deltaG_T2(delta_r_H['deltaH_III-ox'], delta_r_S['deltaS_III-ox'], delta_r_Cp['deltaCp_III-ox'], constants['T'], T),           # [J/mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn + H2O
        
        'deltaG_IV': deltaG_T2(delta_r_H['deltaH_IV'], delta_r_S['deltaS_IV'], delta_r_Cp['deltaCp_IV'], constants['T'], T),                           # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O 
        'deltaG_V': deltaG_T2(delta_r_H['deltaH_V'], delta_r_S['deltaS_V'], delta_r_Cp['deltaCp_V'], constants['T'], T),                               # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
        
        'deltaG_VIII-eps': deltaG_T2(delta_r_H['deltaH_VIII-eps'], delta_r_S['deltaS_VIII-eps'], delta_r_Cp['deltaCp_VIII-eps'], constants['T'], T),   # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
        'deltaG_VIII': deltaG_T2(delta_r_H['deltaH_VIII'], delta_r_S['deltaS_VIII'], delta_r_Cp['deltaCp_VIII'], constants['T'], T),                   # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
        'deltaG_VIII-ox': deltaG_T2(delta_r_H['deltaH_VIII-ox'], delta_r_S['deltaS_VIII-ox'], delta_r_Cp['deltaCp_VIII-ox'], constants['T'], T),       # [J/mol] - Zn^2+ + 2OH^- <--> ZnO + H2O
        
        'deltaG_IX-eps': deltaG_T2(delta_r_H['deltaH_IX-eps'], delta_r_S['deltaS_IX-eps'], delta_r_Cp['deltaCp_IX-eps'], constants['T'], T),           # [J/mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-
        'deltaG_IX': deltaG_T2(delta_r_H['deltaH_IX'], delta_r_S['deltaS_IX'], delta_r_Cp['deltaCp_IX'], constants['T'], T),                           # [J/mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-
        'deltaG_IX-ox': deltaG_T2(delta_r_H['deltaH_IX-ox'], delta_r_S['deltaS_IX-ox'], delta_r_Cp['deltaCp_IX-ox'], constants['T'], T),               # [J/mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^-
        
        'deltaG_X': deltaG_T2(delta_r_H['deltaH_X'], delta_r_S['deltaS_X'], delta_r_Cp['deltaCp_X'], constants['T'], T),                               # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
        
        'deltaG_XI-eps': deltaG_T2(delta_r_H['deltaH_XI-eps'], delta_r_S['deltaS_XI-eps'], delta_r_Cp['deltaCp_XI-eps'], constants['T'], T),           # [J/mol] - Zn(OH)2(s) + 2OH^- <--> Zn(OH)4^2-
        'deltaG_XI': deltaG_T2(delta_r_H['deltaH_XI'], delta_r_S['deltaS_XI'], delta_r_Cp['deltaCp_XI'], constants['T'], T),                           # [J/mol] - Zn(OH)2(aq) + 2OH^- <--> Zn(OH)4^2-
        'deltaG_XI-ox': deltaG_T2(delta_r_H['deltaH_XI-ox'], delta_r_S['deltaS_XI-ox'], delta_r_Cp['deltaCp_XI-ox'], constants['T'], T),               # [J/mol] - ZnO(s) + H2O + 2OH^- <--> Zn(OH)4^2-
        
        'deltaG_XII': deltaG_T2(delta_r_H['deltaH_XII'],delta_r_S['deltaS_XII'], delta_r_Cp['deltaCp_XII'], constants['T'], T),                        # [J/mol] - ZnO(s) + H2O(l) <--> Zn(OH)2(aq)

        'deltaG_HER': deltaG_T2(delta_r_H['deltaH_HER'], delta_r_S['deltaS_HER'], delta_r_Cp['deltaCp_HER'], constants['T'], T),                       # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
        'deltaG_OER': deltaG_T2(delta_r_H['deltaH_OER'], delta_r_S['deltaS_OER'], delta_r_Cp['deltaCp_OER'], constants['T'], T),                       # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
        'deltaG_W': deltaG_T2(delta_r_H['deltaH_W'], delta_r_S['deltaS_W'], delta_r_Cp['deltaCp_W'], constants['T'], T),                               # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
    }

    ## Gibbs free energy using the van't Hoff approximation
    deltaG_Vant_Hoff = {
        'deltaG_I': vant_Hoff(delta_r_G['deltaG_I'], delta_r_H['deltaH_I'], constants['T'], T),                       # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
        'deltaG_II': vant_Hoff(delta_r_G['deltaG_II'], delta_r_H['deltaH_II'], constants['T'], T),                    # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
        
        'deltaG_III-eps': vant_Hoff(delta_r_G['deltaG_III-eps'], delta_r_H['deltaH_III-eps'], constants['T'], T),     # [J/mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
        'deltaG_III': vant_Hoff(delta_r_G['deltaG_III'], delta_r_H['deltaH_III'], constants['T'], T),                 # [J/mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + H2O
        'deltaG_III-ox': vant_Hoff(delta_r_G['deltaG_III-ox'], delta_r_H['deltaH_III-ox'], constants['T'], T),        # [J/mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn + H2O
        
        'deltaG_IV': vant_Hoff(delta_r_G['deltaG_IV'], delta_r_H['deltaH_IV'], constants['T'], T),                    # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O
        'deltaG_V': vant_Hoff(delta_r_G['deltaG_V'], delta_r_H['deltaH_V'], constants['T'], T),                       # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
        
        'deltaG_VIII-eps': vant_Hoff(delta_r_G['deltaG_VIII-eps'], delta_r_H['deltaH_VIII-eps'], constants['T'], T),  # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
        'deltaG_VIII': vant_Hoff(delta_r_G['deltaG_VIII'], delta_r_H['deltaH_VIII'], constants['T'], T),              # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2
        'deltaG_VIII-ox': vant_Hoff(delta_r_G['deltaG_VIII-ox'], delta_r_H['deltaH_VIII-ox'], constants['T'], T),     # [J/mol] - Zn^2+ + 2OH^- <--> ZnO + H2O
        
        'deltaG_IX-eps': vant_Hoff(delta_r_G['deltaG_IX-eps'], delta_r_H['deltaH_IX-eps'], constants['T'], T),        # [J/mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-
        'deltaG_IX': vant_Hoff(delta_r_G['deltaG_IX'], delta_r_H['deltaH_IX'], constants['T'], T),                    # [J/mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-
        'deltaG_IX-ox': vant_Hoff(delta_r_G['deltaG_IX-ox'], delta_r_H['deltaH_IX-ox'], constants['T'], T),           # [J/mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^-
        
        'deltaG_X': vant_Hoff(delta_r_G['deltaG_X'], delta_r_H['deltaH_X'], constants['T'], T),                       # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
        
        'deltaG_XI-eps': vant_Hoff(delta_r_G['deltaG_XI-eps'], delta_r_H['deltaH_XI-eps'], constants['T'], T),        # [J/mol] - Zn(OH)2(s) + 2OH^- <--> Zn(OH)4^2-
        'deltaG_XI': vant_Hoff(delta_r_G['deltaG_XI'], delta_r_H['deltaH_XI'], constants['T'], T),                    # [J/mol] - Zn(OH)2(aq) + 2OH^- <--> Zn(OH)4^2-
        'deltaG_XI-ox': vant_Hoff(delta_r_G['deltaG_XI-ox'], delta_r_H['deltaH_XI-ox'], constants['T'], T),           # [J/mol] - ZnO(s) + H2O + 2OH^- <--> Zn(OH)4^2-
        
        'deltaG_XII': vant_Hoff(delta_r_G['deltaG_XII'], delta_r_H['deltaH_XII'], constants['T'], T),                 # [J/mol] - ZnO(s) + H2O(l) <--> Zn(OH)2(aq)

        'deltaG_HER': vant_Hoff(delta_r_G['deltaG_HER'], delta_r_H['deltaH_HER'], constants['T'], T),                 # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
        'deltaG_OER': vant_Hoff(delta_r_G['deltaG_OER'], delta_r_H['deltaH_OER'], constants['T'], T),                 # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
        'deltaG_W': vant_Hoff(delta_r_G['deltaG_W'], delta_r_H['deltaH_W'], constants['T'], T),                       # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
    }

    ## DeltaG using enthalpy and entropy as weak functions of temperature
    deltaG_approx = {
        'deltaG_I': deltaG_weak(delta_r_H['deltaH_I'], delta_r_S['deltaS_I'], T),                     # [J/mol] - Zn^2+ + 2e^- <--> Zn(s)
        'deltaG_II': deltaG_weak(delta_r_H['deltaH_II'], delta_r_S['deltaS_II'], T),                  # [J/mol] - Zn(OH)^+ + 2e^- + H^+ <--> Zn(s) + H2O
        
        'deltaG_III-eps': deltaG_weak(delta_r_H['deltaH_III-eps'], delta_r_S['deltaS_III-eps'], T),   # [J/mol] - Zn(OH)2(s) + 2e^- + 2H^+ <--> Zn(s) + H2O
        'deltaG_III': deltaG_weak(delta_r_H['deltaH_III'], delta_r_S['deltaS_III'], T),               # [J/mol] - Zn(OH)2(aq) + 2e^- + 2H^+ <--> Zn(s) + H2O
        'deltaG_III-ox': deltaG_weak(delta_r_H['deltaH_III-ox'], delta_r_S['deltaS_III-ox'], T),      # [J/mol] - ZnO(s) + 2e^- + 2H^+ <--> Zn + H2O
        
        'deltaG_IV': deltaG_weak(delta_r_H['deltaH_IV'], delta_r_S['deltaS_IV'], T),                  # [J/mol] - Zn(OH)3^- + 2e^- + 3H^+ <--> Zn(s) + 3H2O 
        'deltaG_V': deltaG_weak(delta_r_H['deltaH_V'], delta_r_S['deltaS_V'], T),                     # [J/mol] - Zn(OH)4^2- + 2e^- + 4H^+ <--> Zn(s) + 4H2O
        
        'deltaG_VIII-eps': deltaG_weak(delta_r_H['deltaH_VIII-eps'], delta_r_S['deltaS_VIII-eps'], T),# [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(s)
        'deltaG_VIII': deltaG_weak(delta_r_H['deltaH_VIII'], delta_r_S['deltaS_VIII'], T),            # [J/mol] - Zn^2+ + 2OH^- <--> Zn(OH)2(aq)
        'deltaG_VIII-ox': deltaG_weak(delta_r_H['deltaH_VIII-ox'], delta_r_S['deltaS_VIII-ox'], T),   # [J/mol] - Zn^2+ + 2OH^- <--> ZnO(s) + H2O
        
        'deltaG_IX-eps': deltaG_weak(delta_r_H['deltaH_IX-eps'], delta_r_S['deltaS_IX-eps'], T),      # [J/mol] - Zn(OH)2(s) + OH^- <--> Zn(OH)3^-
        'deltaG_IX': deltaG_weak(delta_r_H['deltaH_IX'], delta_r_S['deltaS_IX'], T),                  # [J/mol] - Zn(OH)2(aq) + OH^- <--> Zn(OH)3^-
        'deltaG_IX-ox': deltaG_weak(delta_r_H['deltaH_IX-ox'], delta_r_S['deltaS_IX-ox'], T),         # [J/mol] - ZnO(s) + H2O + OH^- <--> Zn(OH)3^-
        
        'deltaG_X': deltaG_weak(delta_r_H['deltaH_X'], delta_r_S['deltaS_X'], T),                     # [J/mol] - Zn(OH)3^- + OH^- <--> Zn(OH)4^-2
        
        'deltaG_XI-eps': deltaG_weak(delta_r_H['deltaH_XI-eps'], delta_r_S['deltaS_XI-eps'], T),      # [J/mol] - Zn(OH)2(s) + 2OH^- <--> Zn(OH)4^2-
        'deltaG_XI': deltaG_weak(delta_r_H['deltaH_XI'], delta_r_S['deltaS_XI'], T),                  # [J/mol] - Zn(OH)2(aq) + 2OH^- <--> Zn(OH)4^2-
        'deltaG_XI-ox': deltaG_weak(delta_r_H['deltaH_XI-ox'], delta_r_S['deltaS_XI-ox'], T),         # [J/mol] - ZnO(s) + H2O + 2OH^- <--> Zn(OH)4^2-
        
        'deltaG_XII': deltaG_weak(delta_r_H['deltaH_XII'], delta_r_S['deltaS_XII'], T),               # [J/mol] - ZnO(s) + H2O(l) <--> Zn(OH)2(aq)

        'deltaG_HER': deltaG_weak(delta_r_H['deltaH_HER'], delta_r_S['deltaS_HER'], T),               # [J/mol] - 4H^+ + 4e^- <--> 2H2(g)
        'deltaG_OER': deltaG_weak(delta_r_H['deltaH_OER'], delta_r_S['deltaS_OER'], T),               # [J/mol] - O2(g) + 4e^- + 4H^+ <--> 2H2O
        'deltaG_W': deltaG_weak(delta_r_H['deltaH_W'], delta_r_S['deltaS_W'], T),                     # [J/mol] - H2O(l) <--> H^+(aq) + OH^-(aq)
    }

    ## E0 using van't approximation
    E0_vant_Hoff = {
        'EI_0': -(deltaG_Vant_Hoff['deltaG_I'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn^2+(aq) --> Zn(s)
        'EII_0': -(deltaG_Vant_Hoff['deltaG_II'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)^+ --> Zn(s) - Not used
        
        'EIII_0-eps': -(deltaG_Vant_Hoff['deltaG_III-eps'])/(constants['n']*constants['F']),# [V vs SHE] - Zn(OH)2(s) --> Zn(s) 
        'EIII_0': -(deltaG_Vant_Hoff['deltaG_III'])/(constants['n']*constants['F']),        # [V vs SHE] - Zn(OH)2(aq) --> Zn(s)
        'EIII_0-ox': -(deltaG_Vant_Hoff['deltaG_III-ox'])/(constants['n']*constants['F']),  # [V vs SHE] - ZnO(s) --> Zn(s)
        
        'EIV_0': -(deltaG_Vant_Hoff['deltaG_IV'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)3^- --> Zn(s)
        'EV_0': -(deltaG_Vant_Hoff['deltaG_V'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
        
        'EHER_0':-(deltaG_Vant_Hoff['deltaG_HER'])/(2*constants['n']*constants['F']),       # [V vs SHE] - HER
        'EOER_0': -(deltaG_Vant_Hoff['deltaG_OER'])/(2*constants['n']*constants['F']),      # [V vs SHE] - OER
    }

    ## E0 using assuming that enthalpy and entropy are weak functions of  temperature
    E0_approx = {
        'EI_0': -(deltaG_approx['deltaG_I'])/(constants['n']*constants['F']),               # [V vs SHE] - Zn^2+(aq) --> Zn(s)
        'EII_0': -(deltaG_approx['deltaG_II'])/(constants['n']*constants['F']),             # [V vs SHE] - Zn(OH)^+ --> Zn(s)
        
        'EIII_0-eps': -(deltaG_approx['deltaG_III-eps'])/(constants['n']*constants['F']),   # [V vs SHE] - Zn(OH)2(s) --> Zn(s) 
        'EIII_0': -(deltaG_approx['deltaG_III'])/(constants['n']*constants['F']),           # [V vs SHE] - Zn(OH)2(aq) --> Zn(s) 
        'EIII_0-ox': -(deltaG_approx['deltaG_III-ox'])/(constants['n']*constants['F']),     # [V vs SHE] - ZnO(s) --> Zn(s)
        
        'EIV_0': -(deltaG_approx['deltaG_IV'])/(constants['n']*constants['F']),             # [V vs SHE] - Zn(OH)3^- --> Zn(s)
        'EV_0': -(deltaG_approx['deltaG_V'])/(constants['n']*constants['F']),               # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
        
        'EHER_0':-(deltaG_approx['deltaG_HER'])/(2*constants['n']*constants['F']),          # [V vs SHE] - HER
        'EOER_0': -(deltaG_approx['deltaG_OER'])/(2*constants['n']*constants['F']),         # [V vs SHE] - OER
    }

    ## E0 with no assumptions, except that the heat capacities are valid for this range
    E0_new = {
        'EI_0': -(deltaG_new_T['deltaG_I'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn^2+(aq) --> Zn(s)
        'EII_0': -(deltaG_new_T['deltaG_II'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)^+ --> Zn(s)
        
        'EIII_0-eps': -(deltaG_new_T['deltaG_III-eps'])/(constants['n']*constants['F']),# [V vs SHE] - Zn(OH)2(s) --> Zn(s)
        'EIII_0': -(deltaG_new_T['deltaG_III'])/(constants['n']*constants['F']),        # [V vs SHE] - Zn(OH)2(aq) --> Zn(s)
        'EIII_0-ox': -(deltaG_new_T['deltaG_III-ox'])/(constants['n']*constants['F']),  # [V vs SHE] - ZnO --> Zn(s) - Not used
        
        'EIV_0': -(deltaG_new_T['deltaG_IV'])/(constants['n']*constants['F']),          # [V vs SHE] - Zn(OH)3^- --> Zn(s)
        'EV_0': -(deltaG_new_T['deltaG_V'])/(constants['n']*constants['F']),            # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
        
        'EHER_0':-(deltaG_new_T['deltaG_HER'])/(2*constants['n']*constants['F']),       # [V vs SHE] - HER
        'EOER_0': -(deltaG_new_T['deltaG_OER'])/(2*constants['n']*constants['F'])       # [V vs SHE] - OER
    }

    ## Using that the derivative of the standard reduction potential with respect to temperature is the change in entropy 
    E0_S = {
        'EI_0': E0_2(constants_E0['EI_0'], delta_r_S['deltaS_I'], constants['T'], T, constants['n']),                    # [V vs SHE] - Zn^2+(aq) --> Zn(s)
        'EII_0': E0_2(constants_E0['EII_0'], delta_r_S['deltaS_II'], constants['T'], T, constants['n']),                 # [V vs SHE] - Zn(OH)^+ --> Zn(s)
        
        'EIII_0-eps': E0_2(constants_E0['EIII_0-eps'], delta_r_S['deltaS_III-eps'], constants['T'], T, constants['n']),  # [V vs SHE] - Zn(OH)2(s) --> Zn(s)
        'EIII_0': E0_2(constants_E0['EIII_0'], delta_r_S['deltaS_III'], constants['T'], T, constants['n']),              # [V vs SHE] - Zn(OH)2(aq) --> Zn(s)
        'EIII_0-ox': E0_2(constants_E0['EIII_0-ox'], delta_r_S['deltaS_III-ox'], constants['T'], T, constants['n']),     # [V vs SHE] - ZnO(s) --> Zn(s) - Not used
        
        'EIV_0': E0_2(constants_E0['EIV_0'], delta_r_S['deltaS_IV'], constants['T'], T, constants['n']),                 # [V vs SHE] - Zn(OH)3^- --> Zn(s)
        'EV_0': E0_2(constants_E0['EV_0'], delta_r_S['deltaS_V'], constants['T'], T, constants['n']),                    # [V vs SHE] - Zn(OH)4^2- --> Zn(s)
        
        'EHER_0': E0_2(constants_E0['EHER_0'], delta_r_S['deltaS_HER'], constants['T'], T, 2*constants['n']),            # [V vs SHE] - HER
        'EOER_0': E0_2(constants_E0['EOER_0'], delta_r_S['deltaS_OER'], constants['T'], T, 2*constants['n'])             # [V vs SHE] - OER
    }

    # The pZn where the equilibrium Zn(OH)2(s) <--> Zn(OH)3^-1 and Zn(OH)3^-1 <--> Zn(OH)4^2- becomes the same and the domain for Zn(OH)3^- vanishes
    pZn_threshold_eps = (deltaG_new_T['deltaG_IX-eps']-deltaG_new_T['deltaG_X'])/(constants['R']*T*np.log(10))
    # The pZn where the equilibrium Zn(OH)2(aq) <--> Zn(OH)3^-1 and Zn(OH)3^-1 <--> Zn(OH)4^2- becomes the same and the domain for Zn(OH)3^- vanishes
    pZn_threshold = (deltaG_new_T['deltaG_IX']-deltaG_new_T['deltaG_X'])/(constants['R']*T*np.log(10))
    # The pZn where the equilibrium ZnO <--> Zn(OH)3^-1 and Zn(OH)3^-1 <--> Zn(OH)4^2- becomes the same and the domain for Zn(OH)3^- vanishes
    pZn_threshold_ox = (deltaG_new_T['deltaG_IX-ox'] - deltaG_new_T['deltaG_X'])/(constants['R']*T*np.log(10))
    # Threshold for stability of ZnO(s)/Zn(OH)2(aq) (pZn higher than this and we have no passivation with ZnO)
    pZn_threshold_pass = (deltaG_new_T['deltaG_XII'])/(constants['R']*T*np.log(10))

    #print(f'Threshold where ZnOH3^- dissapears (when ZnOH2 is dominating): {pZn_threshold}')
    #print(f'Threshold where ZnOH3^- dissapears (when ZnO is dominating): {pZn_threshold_ox}')
    #print(f'Threshold where ZnOH2 becomes the most dominant: {pZn_threshold_pass}')

    # Neutral pH
    Kw = np.exp(-deltaG_new_T['deltaG_W']/(constants['R']*T))
    pKw = -np.log10(Kw)
    pH_neutral = (1/2)*pKw

    ############################################## PLOTTING ##############################################

    pH = np.arange(0, 16, 0.01)     # pH range

    # Defining the lines for the HER and OER
    slope = - 2*constants['R']*T/(constants['n']*constants['F'])*np.log(10)
    EHER = E_OER_HER(E_0 = E0_new['EHER_0'], slope = - 2*constants['R']*T/(constants['n']*constants['F'])*np.log(10), pH=pH)   # HER
    EOER = E_OER_HER(E_0 = E0_new['EOER_0'], slope = - 2*constants['R']*T/(constants['n']*constants['F'])*np.log(10), pH=pH)   # OER

    pZn_values_dict = {'pZn': pZn_value}   # Dictionary with values
    linestyle = ['-', '--']

    ## Printing the Pourbaix diagram

    pZn_string = 'pZn = {}'.format(pZn_value) # The string printed for the different pZn values in the plot. 

    # The amount of Zn in solution
    constants_p = {         # Setting all values to be the same because they are the dominating species within their domain
    'pZn': pZn_value,       # pZn = -log(Zn^2+)
    'pZnOH': pZn_value,     # pZnOH = -log(Zn(OH)^+)
    'pZnOH2': pZn_value,    # pZnOH2 = -log(Zn(OH)2)
    'pZnOH3': pZn_value,    # pZnOH3 = -log(Zn(OH)3^-)
    'pZnOH4': pZn_value,    # pZnOH4 = -log(Zn(OH)4^-2)
    }

    ## E calculations

    # Zn^2+ --> Zn(s)
    EI = np.ones(len(pH))*E0_new['EI_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZn'])

    # Zn(OH)^+ --> Zn(s)
    EII = E0_new['EII_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH']+pH)

    # Zn(OH)2(s) --> Zn(s) ----- Solid hydroxide 
    EIII_eps = E0_new['EIII_0-eps'] - 2*constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(pH)

    # Zn(OH)2(aq) --> Zn(s) ----- Soluble hydroxide
    EIII = E0_new['EIII_0'] - 2*constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH2']/2+pH)

    # ZnO(s) --> Zn(s) ---- Solid oxide
    EIII_ox = E0_new['EIII_0-ox'] - 2*constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(pH)

    # Zn(OH)3^- --> Zn(s)
    EIV = E0_new['EIV_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH3']+3*pH)

    # Zn(OH)4^2- --> Zn(s)
    EV = E0_new['EV_0'] - constants['R']*T*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH4']+4*pH)


    ## pH calculations

    # Zn^2+ --> Zn(OH)2(s) ---- Solid hydroxide
    pHVIII_eps = pKw + deltaG_new_T['deltaG_VIII-eps']/(2*constants['R']*T*np.log(10)) + constants_p['pZn']/2
    # Zn^2+ --> Zn(OH)2(aq) ---- Soluble hydroxide
    pHVIII = pKw + deltaG_new_T['deltaG_VIII']/(2*constants['R']*T*np.log(10)) + constants_p['pZn']/2 - constants_p['pZnOH2']/2
    # Zn^2+ --> ZnO(s) ---- Oxide oxide
    pHVIII_ox = pKw + deltaG_new_T['deltaG_VIII-ox']/(2*constants['R']*T*np.log(10)) + constants_p['pZn']/2

    # Zn(OH)2(s) --> Zn(OH)3^- ---- Using Solid Zn(OH)2
    pHIX_eps = pKw + deltaG_new_T['deltaG_IX-eps']/(constants['R']*T*np.log(10)) - constants_p['pZnOH3']
    # Zn(OH)2(aq) --> Zn(OH)3^- ---- Using Soluble Zn(OH)2
    pHIX = pKw + deltaG_new_T['deltaG_IX']/(constants['R']*T*np.log(10)) + constants_p['pZnOH2'] - constants_p['pZnOH3']
    # ZnO(s) --> Zn(OH)3^- ---- Using oxide ZnO
    pHIX_ox = pKw + deltaG_new_T['deltaG_IX-ox']/(constants['R']*T*np.log(10)) - constants_p['pZnOH3']

    # Zn(OH)3^- --> Zn(OH)4^2-
    pHX = pKw + deltaG_new_T['deltaG_X']/(constants['R']*T*np.log(10)) + constants_p['pZnOH3'] - constants_p['pZnOH4']

    # Zn(OH)2(s) --> Zn(OH)4^2- ---- Using the solid form
    pHXI_eps = pKw + deltaG_new_T['deltaG_XI-eps']/(2*constants['R']*T*np.log(10)) - (1/2)*constants_p['pZnOH4']
    # Zn(OH)2(aq) --> Zn(OH)4^2- ---- Using the Soluble form
    pHXI = pKw + deltaG_new_T['deltaG_XI-eps']/(2*constants['R']*T*np.log(10)) + (1/2)*constants_p['pZnOH2']- (1/2)*constants_p['pZnOH4']
    # ZnO(s) --> Zn(OH)4^2- ---- Using the oxide form
    pHXI_ox = pKw + deltaG_new_T['deltaG_XI-ox']/(2*constants['R']*T*np.log(10)) - (1/2)*constants_p['pZnOH4']


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

    # Initialising the figure
    fig, ax = plt.subplots()

    ## If the concentration of pZn is higher than the threshold for having Zn(OH)3^- (meaning pZn is lower than pZn_threshold)
    if pZn_value < pZn_threshold_ox:
        # Plotting the lines of the Pourbaix diagram
        ax.plot(pH[:index_I_III_ox], EI[:index_I_III_ox], 'k', label='Zn$^{2+}$ - Zn', linestyle=linestyle[0])          # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        ax.vlines(pHVIII_ox, EI[index_I_III_ox], 1.5, 'k', label='Zn$^{2+}$ - ZnO', linestyle=linestyle[0])             # Chemical        - Equilibrium between Zn^2+ and ZnO
        ax.plot(pH[index_I_III_ox: index_III_ox_V], EIII_ox[index_I_III_ox: index_III_ox_V], 'k', label='ZnO - Zn')     # Electrochemical - Equilibrium between ZnO and Zn
        ax.plot(pH[index_III_ox_V:], EV[index_III_ox_V:], 'k', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[0])  # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn(s)
        ax.vlines(pHXI_ox, EV[index_III_ox_V], 1.5, 'k', label='ZnO - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[0])       # Chemical        - Equilibrium between ZnO(s) and Zn(OH)4^2-
        # Colouring the Zn(s) domain
        text_coord_Zn = [np.mean([np.min(pH), pH_neutral]),
                            np.mean([EI[0], -1.5])]
        if EV[-1] <= -1.5:
            add_polygon(ax=ax, 
                        vertices=[(np.min(pH), -1.5), (np.max(pH), EV[-1]), (pHXI_ox, EV[index_III_ox_V]), (pHVIII_ox, EI[index_I_III_ox]), (np.min(pH), EI[0])],
                        colour='black', alpha=0.25, label='Zn(s)', text_coord=text_coord_Zn, text='Zn(s)', text_rotation=0)
        elif EV[-1] > -1.5:
            add_polygon(ax=ax, 
                        vertices=[(np.min(pH), -1.5), (np.max(pH), -1.5), (np.max(pH), EV[-1]), (pHXI_ox, EV[index_III_ox_V]), (pHVIII_ox, EI[index_I_III_ox]), (np.min(pH), EI[0])],
                        colour='black', alpha=0.25, label='Zn(s)', text_coord=text_coord_Zn, text='Zn(s)', text_rotation=0)

        # Colouring the ZnO(s) domain
        text_coord_ZnO = [np.mean([pHVIII_ox, pHXI_ox])-0.7, 
                        np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHVIII_ox, pHXI_ox])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHVIII_ox, pHXI_ox]))])-0.044]
        add_polygon(ax, 
                    [(pHVIII_ox, EIII_ox[index_I_III_ox]), (pHXI_ox, EIII_ox[index_III_ox_V]), (pHXI_ox, 1.5), (pHVIII_ox, 1.5)],
                    colour='black', alpha=0.25, label='ZnO(s)', text_coord=text_coord_ZnO, text='ZnO(s)', text_rotation=0)
        
        # Adding pZn text 
        #ax.text(pH[index_I_III_ox+15], 0.9, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        #ax.text(pH[index_III_ox_V+15], 0.9, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)4^2 -- must be fixed

        # Adding labels to different domains
        # Zn^2+ domain 
        ax.text(np.min([np.mean([np.min(pH), pH_neutral]), np.mean([np.min(pH), pHVIII_ox])])-1,
                np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([np.min(pH), pHVIII_ox])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([np.min(pH), pHVIII_ox]))])-0.05,
                'Zn$^{2+}$(aq)', fontsize=12, color='black')
        # Zn(OH)4^2- domain
        ax.text(np.mean([pHXI_ox, np.max(pH)])-0.25,
                np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHXI_ox, np.max(pH)])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHXI_ox, np.max(pH)]))])-0.4,
                'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain

    ## If the concentration of Zn species is low enough to produce Zn(OH)3^-, but still high enough for passivation
    elif pZn_value >= pZn_threshold_ox and pZn_value < pZn_threshold_pass: # Checking if pZn is bigger than the threshold
        # Plotting the lines of the Pourbaix diagram
        ax.plot(pH[:index_I_III_ox], EI[:index_I_III_ox], 'k', label='Zn$^{2+}$ - Zn', linestyle=linestyle[0])                      # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        ax.vlines(pHVIII_ox, EI[index_I_III_ox], 1.5, 'k', label='Zn$^{2+}$ - ZnO', linestyle=linestyle[0])                         # Chemical        - Equilibrium between Zn^2+ and ZnO(s)
        ax.vlines(pHX, EV[index_IV_V], 1.5, 'k', label='Zn(OH)$_{3}^-$ - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[0])                # Chemical        - Equilibrium between Zn(OH)3^- and Zn(OH)4^2-
        ax.vlines(pHIX_ox, EIV[index_III_ox_IV], 1.5, 'k', label='ZnO - Zn(OH)$_{3}^-$', linestyle=linestyle[0])                    # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)3^-
        ax.plot(pH[index_III_ox_IV:index_IV_V], EIV[index_III_ox_IV:index_IV_V], 'k', label='Zn(OH)$_{3}^{-}$ - Zn', linestyle=linestyle[0]) # Electrochemical - Equilibrium between Zn(OH)3^- and Zn
        ax.plot(pH[index_IV_V:], EV[index_IV_V:], 'k', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[0])                      # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn   
        ax.plot(pH[index_I_III_ox: index_III_ox_IV], EIII_ox[index_I_III_ox: index_III_ox_IV], 'k', label='ZnO - Zn')               # Electrochemical - Equilibrium between Zn(OH)2 and Zn
        
        # Colouring the Zn(s) domain
        text_coord_Zn = [np.mean([np.min(pH), pH_neutral]),
                            np.mean([EI[0], -1.5])]
        if EV[-1] <= -1.5:
            add_polygon(ax=ax, 
                        vertices=[(np.min(pH), -1.5), (np.max(pH), EV[-1]), (pHX, EV[index_IV_V]), (pHIX_ox, EIV[index_III_ox_IV]),  (pHVIII_ox, EI[index_I_III_ox]), (np.min(pH), EI[0])],
                        colour='black', alpha=0.25, label='Zn(s)', text_coord=text_coord_Zn, text='Zn(s)', text_rotation=0)
        elif EV[-1] > -1.5:
            add_polygon(ax=ax, 
                        vertices=[(np.min(pH), -1.5), (np.max(pH), -1.5), (np.max(pH), EV[-1]), (pHX, EV[index_IV_V]), (pHIX_ox, EIV[index_III_ox_IV]),  (pHVIII_ox, EI[index_I_III_ox]), (np.min(pH), EI[0])],
                        colour='black', alpha=0.25, label='Zn(s)', text_coord=text_coord_Zn, text='Zn(s)', text_rotation=0)
        # Colouring the ZnO(s) domain
        text_coord_ZnO = [np.mean([pHVIII_ox, pHIX_ox])-0.7, 
                        np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHVIII_ox, pHIX_ox])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHVIII_ox, pHIX_ox]))])-0.044]
        add_polygon(ax, 
                    [(pHVIII_ox, EIII_ox[index_I_III_ox]), (pHIX_ox, EIII_ox[index_III_ox_IV]), (pHIX_ox, 1.5), (pHVIII_ox, 1.5)],
                    colour='black', alpha=0.25, label='ZnO(s)', text_coord=text_coord_ZnO, text='ZnO(s)', text_rotation=0)

        # Adding pZn text 
        #ax.text(pH[index_I_III_ox+15], 0.9, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        #ax.text(pH[index_IV_V+15], 0.9, pZn_string, color='k', rotation = 90)       # Adding pZn value to equilibrium between Zn(OH)3^-1 and Zn(OH)4^2
        # Checking if there is room for the letters and sets the label accordingly 
        # if pHX - pHIX_ox < 0.5:
        #     ax.text(pH[index_III_ox_IV-50], 0.9, pZn_string, color='k', rotation = 90)     # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)3^-
        # elif pHX - pHIX_ox >= 0.5:
        #     ax.text(pH[index_III_ox_IV+15], 0.9, pZn_string, color='k', rotation = 90)     # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)3^-
        
        # Adding text to different domains

        # Zn^2+ domain
        ax.text(np.min([np.mean([np.min(pH), pH_neutral]), np.mean([np.min(pH), pHVIII_ox])])-1,
                np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([np.min(pH), pHVIII_ox])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([np.min(pH), pHVIII_ox]))])-0.05,
                'Zn$^{2+}$(aq)', fontsize=12, color='black')
        # Zn(OH)4^2- domain
        ax.text(np.mean([pHX, np.max(pH)])-0.25,
                np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHX, np.max(pH)])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHX, np.max(pH)]))])-0.4,
                'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)   # Zn(OH)4^2- domain
        
        # Zn(OH)3^- domain 
        # (Checks if it is room for the label inside the domain)
        if pHX-pHIX_ox<1:   
            # Text doesn't fit, use annotation with an arrow
            text_pos = (15, E_OER_HER(E_0=E0_new['EOER_0'], slope=slope, pH=15)+0.2)
            arraw_head_pos = (np.mean([pHIX_ox, pHX]), np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHIX_ox, pHX])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHIX_ox, pHX]))]))
            ax.annotate('Zn(OH)$^{-}_{3}$(aq)', xy=arraw_head_pos, xytext= text_pos,
                        fontsize=12, color='black', rotation=90,
                        arrowprops=dict(facecolor='black', arrowstyle='->'))
        elif pHX-pHIX_ox>=1:
            # If it is room
            ax.text(np.mean([pHIX_ox, pHX])-0.25,
                    np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHIX_ox, pHX])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHIX_ox, pHX]))])-0.4,
                    'Zn(OH)$^{-}_{3}$(aq)', fontsize=12, color='black', rotation=90) 

    ## For the case when the concentration of Zn species is so low that we don't have passivating nature (Produce soluble Zn(OH)2(aq) instead of ZnO(s))
    elif pZn_value >= pZn_threshold_pass:
        # Plotting the lines of the Pourbaix diagram
        ax.plot(pH[:index_I_III], EI[:index_I_III], 'k', label='Zn$^{2+}$ - Zn', linestyle=linestyle[0])                                # Electrochemical - Equilibrium between Zn^2+ and Zn(s) 
        ax.vlines(pHVIII, EI[index_I_III], 1.5, 'k', label='Zn$^{2+}$ - Zn(OH)2(aq)', linestyle=linestyle[0])                           # Chemical        - Equilibrium between Zn^2+ and ZnO(s)
        ax.vlines(pHX, EV[index_IV_V], 1.5, 'k', label='Zn(OH)$_{3}^-$ - Zn(OH)$_{4}^{2-}$', linestyle=linestyle[0])                    # Chemical        - Equilibrium between Zn(OH)3^- and Zn(OH)4^2-
        ax.vlines(pHIX, EIV[index_III_IV], 1.5, 'k', label='Zn(OH)2(aq) - Zn(OH)$_{3}^-$', linestyle=linestyle[0])                      # Chemical        - Equilibrium between Zn(OH)2 and Zn(OH)3^-
        ax.plot(pH[index_III_IV:index_IV_V], EIV[index_III_IV:index_IV_V], 'k', label='Zn(OH)$_{3}^{-}$ - Zn', linestyle=linestyle[0])  # Electrochemical - Equilibrium between Zn(OH)3^- and Zn
        ax.plot(pH[index_IV_V:], EV[index_IV_V:], 'k', label='Zn(OH)$_{4}^{2-}$ - Zn', linestyle=linestyle[0])                          # Electrochemical - Equilibrium between Zn(OH)4^2- and Zn   
        ax.plot(pH[index_I_III: index_III_IV], EIII[index_I_III: index_III_IV], 'k', label='Zn(OH)2(aq) - Zn')                          # Electrochemical - Equilibrium between Zn(OH)2 and Zn
        # Colouring the Zn(s) domain
        text_coord_Zn = [np.mean([np.min(pH), pH_neutral]),
                            np.mean([EI[0], -1.5])]
        if EV[-1] <= -1.5:
            add_polygon(ax=ax, 
                        vertices=[(np.min(pH), -1.5), (np.max(pH), EV[-1]), (pHX, EV[index_IV_V]), (pHIX, EIV[index_III_IV]),  (pHVIII, EI[index_I_III]), (np.min(pH), EI[0])],
                        colour='black', alpha=0.25, label='Zn(s)', text_coord=text_coord_Zn, text='Zn(s)', text_rotation=0)
        elif EV[-1] > -1.5:
            add_polygon(ax=ax, 
                        vertices=[(np.min(pH), -1.5), (np.max(pH), -1.5), (np.max(pH), EV[-1]), (pHX, EV[index_IV_V]), (pHIX, EIV[index_III_IV]),  (pHVIII, EI[index_I_III]), (np.min(pH), EI[0])],
                        colour='black', alpha=0.25, label='Zn(s)', text_coord=text_coord_Zn, text='Zn(s)', text_rotation=0)

        # # Colouring the ZnO(s) domain
        text_coord_ZnOH2 = [np.mean([pHVIII, pHIX]), 
                        np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHVIII, pHIX])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHVIII, pHIX]))])-0.4]
        add_polygon(ax, 
                    [(pHVIII, EIII[index_I_III]), (pHIX, EIII[index_III_IV]), (pHIX, 1.5), (pHVIII, 1.5)],
                    colour='black', alpha=0, label='Zn(OH)$_{2}$(aq)', text_coord=text_coord_ZnOH2, text='Zn(OH)$_{2}$(aq)', text_rotation=90)

        # Adding pZn text
        #ax.text(pH[index_I_III+15], 0.9, pZn_string, color='k', rotation = 90)   # Adding pZn value to equilibrium between Zn^2+ and Zn(OH)2 -- must be fixed
        #ax.text(pH[index_III_IV+15], 0.9, pZn_string, color='k', rotation = 90)  # Adding pZn value to equilibrium between Zn(OH)2 and Zn(OH)3^-
        #ax.text(pH[index_IV_V+15], 0.9, pZn_string, color='k', rotation = 90)    # Adding pZn value to equilibrium between Zn(OH)3^-1 and Zn(OH)4^2
        
        # Adding text to different domains
        
        # Zn^2+ domain
        ax.text(np.min([np.mean([np.min(pH), pH_neutral]), np.mean([np.min(pH), pHVIII])])-1,
                np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([np.min(pH), pHVIII_ox])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([np.min(pH), pHVIII_ox]))])-0.05,
                'Zn$^{2+}$(aq)', fontsize=12, color='black') 
        # Zn(OH)4^2- domain
        ax.text(np.mean([pHX, np.max(pH)])-0.25,
                np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHX, np.max(pH)])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHX, np.max(pH)]))])-0.4,
                'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)
        # Zn(OH)3^- domain
        ax.text(np.mean([pHIX, pHX])-0.15,
                np.mean([E_OER_HER(E0_new['EHER_0'], slope ,pH=np.mean([pHIX, pHX])), E_OER_HER(E0_new['EOER_0'], slope , pH=np.mean([pHIX, pHX]))])-0.4,
                'Zn(OH)$^{-}_{3}$(aq)', fontsize=12, color='black', rotation=90)  # Zn(OH)3^- domain

    # Adding things that are present in the diagram for all cases
    ax.plot(pH,EHER, '--', c='k')     # Line for the HER
    ax.plot(pH,EOER, '--', c='k')     # Line for the OER
    ax.text(2, -0.35, 'HER', fontsize=12, color='black', rotation=-15)  # HER line
    ax.text(2, 0.85, 'OER', fontsize=12, color='black', rotation=-15)   # OER line
    ax.text(1, 1.30, 'T = ' + str(T-273.15) + '$^{o}C$',fontsize=12, color='black', bbox=dict(facecolor='white', alpha=1))
    ax.vlines(pH_neutral, ymin=-1.5, ymax=1.5, colors='k', label='Neutral pH', alpha=0.40, linestyle=linestyle[1])             # Adding a line for neutral
    ax.set_xlabel('pH  /  []')
    ax.set_ylabel('Potential - E  /  V')
    ax.set_xlim(xmin=0, xmax=max(pH))   # Set the x-axis range
    ax.set_ylim(ymin=-1.5, ymax=1.5)    # Set the y-axis range
    st.pyplot(fig=fig)

    st.write(f'At T = {temperature - 273.15} °C we have that:')
    
    st.write('**Passivating effects**')  
    st.write(r'At this temperature, ZnO will be passivating as long as pZn $\leqslant %.2f$' %pZn_threshold_pass)
    st.write(r'After that, the dominating species will be Zn(OH)$_{2}$(aq), and there will not be any passivating effects')

    st.write('**Suppressed domain**')
    st.write(r'At this tempertaure, the formation of Zn(OH)$_{3}^{-}$(aq) is supressed until pZn $< %.2f$ . ' %pZn_threshold_ox)
    st.write(r'After that, the domain for Zn(OH)$_{3}^{-}$(aq) will keep expanding until it reaches an equilibrium with Zn(OH)$_{2}$(aq) at pZn $\approx %.2f$' %pZn_threshold_pass)
    
    st.write('**Neutral pH of water**')
    st.write(r'At this temperature, the neutral $pH$ of water is pH$_{Neutral} = %.2f$' %pH_neutral)


def main():
    # Setting the title for the app
    st.title('Pourbaix Diagram of Zn')

    # Writing an introduction
    Introduction='''Pourbaix Diagrams are a way of displaying the most stable compound/species as a function of **pH** and 
            potential, **E**, at a specific temperature and pressure. It is purely based on thermodynamics and equilibria,
            so it does not tell anything about the rate of which processes occur. Still, it is a popular tool within the
            field of corrosion. Pourbaix Diagrams are usually illustrating a compound (like Zn in this case) in contact
            with water. Since it is in contact with water, the lines for the Hydrogen Evolution Reaction (HER) and 
            Oxygen Evolution Reaction (OER) are also normally included. The shaded areas in the diagram shows the domains
            where we have **immune** Zn(s) and **passivating** ZnO(s) regions, meaning no corrosion.
            '''
    # Creating a subheader
    subheader1 = 'Equilibria'

    # Writing a description for the subheader
    multi1 = '''**Vertical lines** represents a purely chemical equilibrium, meaning that there is no transfer of electrons (no change in oxidation states).
            They are therefore independent of the potential. An example of this is the equilibrium between Zn$^{2+}$(aq) and ZnO(s)  
            \\
            $$ Zn^{2+}(aq) + H_{2}O = ZnO(s) + 2H^{+}(aq) $$  
            \\
            **Horisontal lines** represents an electrochemical equilibrium which is independent of the **pH**, and only a function of the potential.  
            An example of this would be the equilibrium between Zn(s) and Zn$^{2+}$  
            \\
            $$ Zn^{2+} + 2e^{-} = Zn(s) $$  
            \\
            **Sloped lines** represents an electrochemical equilibrium which is also dependent of the **pH**.  
            An example of this would be the equilibrium between the Zn(s) and ZnO(s)  
            \\
            $$ ZnO(s) + 2e^{-} + 2H^{+}(aq) = Zn(s) + H_{2}O(l) $$  
            '''
    
    # Creating a subheader
    subheader2 = 'Concentration and temperature'

    # Writing a description for the subheader
    multi2 = ''' **Concentration**  
            In addition to the **pH** and **potential** defining the equilibria, the amount of dissolved species also shifts equilibria with dissolved Zn-species.
            This is taken into account by the quantity **pZn**. \\
            An example of such an equilibrium is the equilibrium between Zn(s) and Zn$^{2+}$(aq)
            \\
            $$ Zn^{2+}(aq) + 2e^{-} = Zn(s) $$
            \\
            Using the Nernst-equation to express the equilibrium reduction potential for the reactiong gives
            \\
            $$ E^{Rev}(T) = E^{0}(T) - \\frac{RT \\ln{10}}{nF}pZn^{2+} $$
            \\
            In this domiain, most of the dissolved species in the solution will come from Zn$^{2+}$ ions, pZn $\\approx$ pZn$^{2+}$.
            \\
            In this diagram we approximate the activity to be equal to the concentration 
            (divided by a reference concentration 1M). The **pZn** can be visualised like this  
            \\
            $$ pZn = -log(a_{Zn}) \\approx -log\\left(\\frac{c_{Zn}}{c^{0}} \\right) $$  
            \\
            This will for instance imply that if we have a total concentration of Zn-species dissolved to be c$_{Zn}$ = 10$^{-6}$M, then **pZn** = 6.
            A value of 6 is usually regarded as the corrosion limit, and a value of 8 is often used as a measure for ultra pure water [1]
            \\
            \\
            **Temperature**  
            In order to produce this Pourbaix Diagram, thermodynamic data from Beverskog *et al.* [1] was used along with supplementary data
            from SI Chemical Data [2]. By assuming that the heat capacities provided by Beverskog *et al.* [1] are valid within the temperature
            range between 25-100 °C, the temperature dependence could also be implemented. The temperature dependence were implemented by assuming...  
            \\
            $$ \Delta_{r}G^{0}(T_{2}) = \Delta_{r}H^{0}(T_{2}) - T_{2}\Delta_{r}S^{0}(T_{2}) $$   
            \\
            Where $\Delta_{r}H^{0}(T_{2})$ and $\Delta_{r}S^{0}(T_{2})$ are found by  
            \\
            $$ \Delta_{r}H^{0}(T_{2}) = \Delta_{r}H^{0}(T_{1}) + \int_{T_{1}}^{T_{2}} \Delta_{r}C_{p} \,dT  \quad\land\quad \Delta_{r}S^{0}(T_{2}) = \Delta_{r}S^{0}(T_{1}) + \int_{T_{1}}^{T_{2}} \\frac{\Delta_{r}C_{p}}{T} \,dT$$
            \\
            \\
            Here, $\Delta_{r}G$, $\Delta_{r}H$, and $\Delta_{r}S$ represents Gibbs free energy, enthalpy, and entropy, for the reaction.
            T$_{2}$ is symbolising another temperature than the reference temperature T_{1}$ = 25°C. $\Delta_{r}C_{p}$ is the
            heat capacity for the reaction. One can connect the standard reduction reaction for an electrochemical reaction to the Gibbs free energy 
            for the reaction at a specific temperature by  
            \\
            $$ E^{0}(T) = -\\frac{\Delta_{r}G^{0}(T)}{nF} $$
            \\
            \\
            The neutral **pH** of water is also a function of temperature and is represented by the vertical dotted grey line.
            '''
    # Using streamlit to print the text
    st.write(Introduction)
    ## Inserting the sliders and temperature input

    # text input for the pZn with two decimal resolution/step size
    st.subheader('Input values')
    pZn = st.slider('Select a value for **pZn**', min_value=0.0, max_value=8.0, value=3.0, step=0.01)

    # Text input for temperature
    temperature_input = st.text_input('Enter a temperature (25 - 100°C): ', value=25)
    try:
        temperature = float(temperature_input) + 273.15

        # Check if temperature is outside of range
        if temperature <25+273.15 or temperature >100+273.15:
            st.error('Please enter a temperature within the specifiec range')
            return 
    
    except ValueError:
        st.error('Please enter a valid temperature. Must be a positive number between 25-100 °C. Use "." as comma.')
    
        return
    
    ## Plotting the Pourbvaix Diagram
    PourbaixDiagram(pZn=pZn, temperature=temperature)

    st.subheader(subheader1)
    st.markdown(multi1)
    st.subheader(subheader2)
    st.markdown(multi2)

if __name__ == "__main__":
    main()

    references = ''' **References**  
    \\
    [1] B. Beverskog and I. Puigdomenech, “Revised pourbaix diagrams for zinc at 25–300 °c,” Corrosion Science, vol. 39, no. 1, pp. 107–114, 1997. 
    \\
    [2] A. Blackman, “Aylward and findlay’s si chemical data,” 2014. 
    '''
    st.markdown(references)



## References 
#
# [1] B. Beverskog and I. Puigdomenech, “Revised pourbaix diagrams for zinc at 25–300 °c,” Corrosion Science, vol. 39, no. 1, pp. 107–114, 1997.
# [2] A. Blackman, “Aylward and findlay’s si chemical data,” 2014.
#
#
########################################################

### 
#Running the streamlit app:
# - Open Command line (seach for cmd in the start meny)
# - Navigate to the directory with the venv (enter: cd OneDrive - SINTEF\Documents\Prosjekter\ZABAT)
# - Activate the environment (it has streamlit) (enter venv-zabat\Scripts\activate)
# - Navigate to the directory with the script (enter: cd GitHub\ZABAT)
# - Run the script with streamlit (enter: streamlit run PourbaixDiagram_Zn_Streamlit.py or python -m streamlit run PourbaixDiagram_Zn_Streamlit)
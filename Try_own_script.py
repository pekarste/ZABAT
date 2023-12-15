# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 10:22:43 2023

@author: sidselh

Based on chat_gpt.py
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

class ChemicalEquilibrium:
    def __init__(self):
        self.equilibrium_constants = {'ZnOH4': 10**18, 'ZnOH3': 10**13.7, 'HCO3': 10**9.56, 'H2CO3': 10**(6.33),\
                                      'pCO2': 10**(-1.55), 'CO3': 10**(-2.0), 'ZnCO3': 10**(-10), 'HF': 10**3.3, 'HF2': 10**3.3,\
                                          'ZnF': 10**0.8, 'ZnO': 10**(-15.96), 'ZnOH': 10**5.0, 'ZnOH2_sat': 10**(-14.82),\
                                              'ZnOH2_aq': 10**8.3, 'H2O': 10**4.18, 'KOH': 10**(-0.2)}  # Initialize with a non-zero value
 
    def system_of_equations(self, concentrations):
        zn_concentration, oh4_concentration, oh3_concentration = concentrations
        eq1 = zn_concentration * oh4_concentration**4 * - self.equilibrium_constants['ZnOH4']
        eq2 = zn_concentration * oh3_concentration**3 - self.equilibrium_constants['ZnOH3']
        eq3 = zn_concentration + 4 * oh4_concentration + 3*oh3_concentration - total_concentration # Why are there 4 and 3 multiples here? IOs not this describing the total concentration of Zn?
        return [eq1, eq2, eq3]
    
#######################################################################################################################
    # def soe(self, concentration, pH):
    #     C_Zn_ion, C_ZnOH3, C_ZnOH4, C_ZnOH, C_Zn_tot = concentration      # Array with the different forms of Zn
        
    #     K_ZnOH4 = self.equilibrium_constants['ZnOH4']           # Equilibrium constant for formation or dissolution of Zn(OH)4
    #     K_ZnOH3 = self.equilibrium_constants['ZnOH3']           # Equilibrium constant for formation or dissolution of Zn(OH)3
    #     Ksp_ZnOH2 = self.equilibrium_constants['ZnOH2_sat']     # Equilibrium constant for formation or dissolution of Zn(OH)2
    #     K_ZnOH = self.equilibrium_constants['ZnOH']             # Equilibrium constant for formation or dissolution of ZnOH+
    #     Ksp_ZnO = self.equilibrium_constants['ZnO']             # Equilibirum constant for formation or dissolution of ZnO
    #     C_OH = 10**(-14 + pH)                                   # Concentration of Hydroxide 

    #     eq1 = C_Zn_ion*C_OH**4 - K_ZnOH4*C_ZnOH4        # Dissolution/formation of Zn(OH)4
    #     eq2 = C_Zn_ion*C_OH**3 - K_ZnOH3*C_ZnOH3        # Dissolution/formation of Zn(OH)3
    #     eq3 = C_Zn_ion*C_OH**2 - Ksp_ZnOH2              # Dissolution/formation of Zn(OH2)
    #     eq4 = C_Zn_ion*C_OH - K_ZnOH*C_ZnOH             # Dissolution/formation of ZnOH
    #     eq5 = C_Zn_ion*C_OH**2 - Ksp_ZnO                # Dissolution/formation of ZnO
    #     eq6 = C_Zn_ion + C_ZnOH4 + C_ZnOH3 + C_ZnOH + C_ZnO + C_ZnOH2 - C_Zn_tot   # The total concentration of Zn species

    #     return [eq1, eq2, eq3, eq4, eq5, eq6]
########################################################################################################################

    # Calculating for hydroxide species
    def calculate_equilibrium_constants(self, ph_values):
        equilibrium_constants = {'ZnOH4': [], 'ZnOH3': [], 'ZnO': []}
        zn_concentrations = []
        oh4_concentrations = []
        oh3_concentrations = []

        for ph in ph_values:
            global total_concentration
            total_concentration = 10 ** (-13.96 + ph)  # Total concentration of H⁺ and OH⁻ --- Perhaps say -14 + pH or something

            # Initial guess for concentrations
            initial_guess = [total_concentration / 5, total_concentration / 5, total_concentration / 5]
            # initial_concentrations = {'KOH': 7.0, 'KF': 1.4, 'K2CO3': 1.4}  # Molar concentration, M = mol/l
            # i_g = [initial_concentrations['KOH'], initial_concentrations['KF'], initial_concentrations['K2CO3']]
            
            # Update the instance variables for equilibrium constants
            options = {'maxfev': 10000, 'xtol': 1e6}
            concentrations = fsolve(self.system_of_equations, initial_guess, **options)
            self.equilibrium_constants['ZnOH4'] = concentrations[0] * concentrations[1]**4
            self.equilibrium_constants['ZnOH3'] = concentrations[0] * concentrations[2]**3

            # Solve the system of equations again to get concentrations
            concentrations = fsolve(self.system_of_equations, initial_guess, **options)

            equilibrium_constants['ZnOH4'].append(self.equilibrium_constants['ZnOH4'])
            equilibrium_constants['ZnOH3'].append(self.equilibrium_constants['ZnOH3'])
            zn_concentrations.append(concentrations[0])
            oh4_concentrations.append(concentrations[1])
            oh3_concentrations.append(concentrations[2])

        return equilibrium_constants, zn_concentrations, oh4_concentrations, oh3_concentrations

# Initialize ChemicalEquilibrium
equilibrium_system = ChemicalEquilibrium()

# pH range from 7.0 to 15.0 at 0.1 increments
ph_range = np.arange(7.0, 15.1, 0.1)

# Calculate equilibrium constants and concentrations over the pH range
equilibrium_constants, zn_concentrations, oh4_concentrations, oh3_concentrations = \
    equilibrium_system.calculate_equilibrium_constants(ph_range)
    
# Additional calculations for concentrations
c_H = 10**(-ph_range)
c_OH = 10**(-13.96) / c_H
c_Zn = equilibrium_constants['ZnO']

# Calculating concentrations
# c_HCO3 = c_H * equilibrium_constants['CO3'] * 10**9.56
# c_H2CO3 = c_H * c_HCO3 * 10**6.33
# p_CO2 = c_H2CO3 * 10**(-1.55)
# c_CO3 = equilibrium_constants['K2CO3'] * 10**(-2.0)     # Need to check Ksp or Keq
# c_ZnCO3 = c_Zn * c_CO3 * 10**(-10)
# c_HF = c_H * c_F * 10**3.3
# c_HF2 = c_HF * c_F * 10**0.86
# c_ZnF+ = c_Zn * c_F- * 10**0.8
# c_ZnO = c_Zn * c_OH**2 * 10**(-15.96)   # Change to ZnO_sat
# c_ZnOH+ = c_Zn * c_OH * 10**5
# c_ZnOH2_sat = c_Zn * c_OH**2 * 10**(-14.82)
# c_ZnOH2_aq = c_Zn * c_OH**2 * 10**8.3
# c_ZnOH3 = c_Zn * c_OH**3 * 10**13.7
# c_ZnOH4 = c_Zn * c_OH**4 * 10**18
# c_H2O = c_H * c_OH * 10**4.18
# c_KOH = c_K * c_OH * 10**(-0.2)   Need to check this!
# KF = K+ + F-

# Plot the results
plt.figure(figsize=(10, 6))

# Plot concentrations of Zn²⁺ and OH⁻
plt.plot(ph_range, zn_concentrations, label='[Zn²⁺]')
plt.plot(ph_range, oh4_concentrations, label='[Zn(OH)₄²⁻]')
plt.plot(ph_range, oh3_concentrations, label='[Zn(OH)₃⁻]')

# Set plot properties
plt.xlabel('pH')
plt.ylabel('Concentration (M)')
plt.title('Concentration vs pH')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()

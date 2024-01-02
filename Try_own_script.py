# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 10:22:43 2023

@author: sidselh

Based on chat_gpt.py
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

class Zn_solution:
    def __init__(self, initial_concentrations):

        '''
        INPUT:
        initial_concentration: contains the initial concentration of Zn^2+, KOH, K2CO3, and KF in an array in mol/L

        Currently, KOH is not used since it would dominate the pH, so the pH is instead set to vary.
        It is also assumed full dissociation of K2CO3 and KF.
        Concentration of Zn^+ is assumed to be a natural occuring concentration and usually very low
        '''
        
        # Dictionary with equilibrium constants
        self.equilibrium_constants = {  'ZnOH4': 10**18, 'ZnOH3': 10**13.7,'ZnOH2_sat': 10**(-14.82),'ZnOH2_aq': 10**8.3, 'ZnOH': 10**5.0, 'ZnO': 10**(-15.96), \
                                        'ZnCO3': 10**(-10), 'H2CO3': 10**(6.33), 'HCO3': 10**(9.56),\
                                        'pCO2': 10**(-1.55), 'HF': 10**3.3, 'HF2': 10**0.86,\
                                        'ZnF': 10**0.8}  # Initialize with a non-zero value
        #'CO3': 10**(-2.0), , 'H2O': 10**4.18, 'KOH': 10**(-0.2)
        
        # Initial concentration of species
        self.c_Zn_2 = initial_concentrations[0]     # Setting the initial concentration of Zn^2+ for the system
        self.c_KOH = initial_concentrations[1]      # Setting the initial concentration of KOH for the system -- Not used
        self.c_K2CO3 = initial_concentrations[2]    # Setting the initial concentration of K2CO3 for the system
        self.c_KF = initial_concentrations[3]       # Setting the initial concentration of KF

        self.c_Zn_tot = self.c_Zn_2                 # The total concentration of Zn species is always equal to the initial concentration of Zn^2+
        self.c_COx_tot = self.c_K2CO3               # The total concentration of COx species is always equal to the initial concentration of CO3^2-
        self.c_K_tot = self.c_KOH + 2*self.c_K2CO3 + self.c_KF
        self.c_F_tot = self.c_KF                    # The total concentration of F species is always equal to the initial concentration of KF

        # Defining the pH
        self.pH_range = np.arange(0, 15.1, 0.1)     # The pH range for the study

        # Defining matrix to keep the concentrations of calculated species
        self.num_species = 21                       # Total number of species in the system
        self.concentration_matrix = np.zeros((self.num_species, len(self.pH_range)))# Matrix for storing the concentrations for different pH values

    def distribute_Zn_solution_species(self, x, pH):
        '''
        INPUT:
        x: An array containing first estimates of the Zn^2+, CO3^2-, K+, and F- (based on the initial value) -- Needed since fsolve is used later
        pH: An array of the pH range

        This part contains the thermodynamic equilibrium equations describing the distribution of the different
        species based on the estimated guesses from x. This is used since fsolve is used later to solve the system

        RETURNS:
        concentration_array: An array of the different concentrations from the system of equilibriums. 
        '''
        
        # Empty array to store solutions
        concentration_array = np.zeros(self.num_species)

        # Concentration of protons and hydroxide from pH
        c_H = 10**(-pH)
        c_OH = 10**(-13.96)/c_H

        # Values from x needed to solve the system
        c_Zn_2 = x[0]
        c_CO3_2 = x[1]
        c_K_1 = x[2]
        c_F_1 = x[3]

        ## System of equations

        # Equilibrium constants
        K_CO2 = self.equilibrium_constants['pCO2']              # Equilibrium constant for formation of CO2
        K_H2CO3 = self.equilibrium_constants['H2CO3']           # Equilibrium constant for formation of H2CO3
        K_HCO3_1 = self.equilibrium_constants['HCO3']           # Equilibrium constant for formation of HCO3^-

        #K_K2CO3 = self.equilibrium_constants['K2CO3']          #- Missing -- K2CO3 gives intial conditions
        K_ZnCO3 = self.equilibrium_constants['ZnCO3']           # Equilibrium constant for formation of ZnCO3

        K_HF2 = self.equilibrium_constants['HF2']               # Equilibrium constant for formation of HF2^-
        K_HF = self.equilibrium_constants['HF']                 # Equilibrium constant for formation of HF
        K_ZnF_1 = self.equilibrium_constants['ZnF']             # Equilibrium constant for formation of ZnF

        K_ZnOH4 = self.equilibrium_constants['ZnOH4']           # Equilibrium constant for formation of Zn(OH)4
        K_ZnOH3 = self.equilibrium_constants['ZnOH3']           # Equilibrium constant for formation of Zn(OH)3
        Ksp_ZnOH2 = self.equilibrium_constants['ZnOH2_sat']     # Equilibrium constant for formation of Zn(OH)2
        K_ZnOH2 = self.equilibrium_constants['ZnOH2_aq']        # Equilibrium constant for formation of Zn(OH)2
        K_ZnOH = self.equilibrium_constants['ZnOH']             # Equilibrium constant for formation of ZnOH+
        Ksp_ZnO = self.equilibrium_constants['ZnO']             # Equilibirum constant for formation of ZnO

        #K_KOH = self.equilibrium_constants['KOH']               # Equilibrium constant for dissociation of KOH
        #K_KF = self.equilibrium_constants['KF']                 # Equilibrium constant for dissociation of KF

        # Carbon acid in water
        c_HCO3_1 = c_CO3_2*c_H*K_HCO3_1            # Formation of HCO3^-1 from CO3^-2 
        c_H2CO3 = c_HCO3_1*c_H*K_H2CO3             # Formation of H2CO3 from HCO3^-1 
        c_CO2 = c_H2CO3*K_CO2                      # Formation of H2CO3

        # Carbonates from salts
        #c_K2CO3 = c_CO3_2*c_K_1**2/K_K2CO3         # Dissolution of K2CO3 -- Complete dissolution?
        c_ZnCO3 = c_CO3_2*c_Zn_2*K_ZnCO3            # Formation of ZnCO3

        # Fluorine species
        c_HF = c_H*c_F_1*K_HF                       # Formation of HF
        c_HF2 = c_HF*c_F_1*K_HF2                    # Formation of HF2
        c_ZnF_1 = c_Zn_2*c_F_1*K_ZnF_1              # Formation of ZnF^-1

        # Zn complexes and salts
        c_ZnOH4 = K_ZnOH4*(c_Zn_2*c_OH**4)          # Formation of Zn(OH)4^-2
        c_ZnOH3 = K_ZnOH3*(c_Zn_2*c_OH**3)          # Formation of Zn(OH)3^-1
        c_ZnOH2 = K_ZnOH2*(c_Zn_2*c_OH**2)          # Formation of Zn(OH)2
        c_ZnOH = K_ZnOH*(c_Zn_2*c_OH)               # Formation of ZnOH^+1
        c_ZnO = Ksp_ZnO*(c_Zn_2*c_OH**2)            # Formation of ZnO

        # Potassium salts
        #c_KOH = c_K_1*c_OH/K_KOH                  # Formation of KOH -- Complete dissolution? - Get from x
        #c_KF = c_F_1*c_K_1/K_KF                   # Formation of KF -- Complete dissolution? - Get from x


        # Initilise concentration array
        concentration_array[0] = c_Zn_2
        concentration_array[1] = c_ZnOH4
        concentration_array[2] = c_ZnOH3
        concentration_array[3] = c_ZnOH2
        concentration_array[4] = c_ZnOH
        concentration_array[5] = c_ZnO
        concentration_array[6] = c_ZnCO3
        concentration_array[7] = c_ZnF_1

        concentration_array[8] = c_CO2
        concentration_array[9] = c_H2CO3
        concentration_array[10] = c_HCO3_1
        concentration_array[11] = c_CO3_2

        concentration_array[12] = c_K_1
        #concentration_array[13] = c_KF
        #concentration_array[14] = c_K2CO3
        #concentration_array[15] = c_KOH
        concentration_array[16] = c_F_1
        concentration_array[17] = c_HF
        concentration_array[18] = c_HF2

        concentration_array[19] = c_H
        concentration_array[20] = c_OH


        return concentration_array
 
    def conservation_Zn_solution(self, x, pH):
        '''
        INPUT
        x: An array containing first estimates of the Zn^2+, CO3^2-, K+, and F- (based on the initial value) -- Needed since fsolve is used later
        pH: An array of the pH range

        This part calculates the conservatuion of different species and uses fsolve solve the system

        RETURN:
        equation_array: An array of the total concentration of Zn-, CO3-, F-, and K-species to match the initial concentration
        '''

        # Getting the concentration array from the system of equilibriums
        concentration_array = self.distribute_Zn_solution_species(x, pH)

        # Extracting the different concentrations from the concentration_array
        c_Zn_2 = concentration_array[0]
        c_ZnOH4 = concentration_array[1]
        c_ZnOH3 = concentration_array[2]
        c_ZnOH2 = concentration_array[3]
        c_ZnOH = concentration_array[4]
        c_ZnO = concentration_array[5]
        c_ZnCO3 = concentration_array[6]
        c_ZnF_1 = concentration_array[7]

        c_CO2 = concentration_array[8]
        c_H2CO3 = concentration_array[9]
        c_HCO3_1 = concentration_array[10]
        c_CO3_2 = concentration_array[11]

        c_K_1 = concentration_array[12]
        #c_KF = concentration_array[13]
        #c_K2CO3 = concentration_array[14]
        #c_KOH = concentration_array[15]
        c_F_1 = concentration_array[16]
        c_HF = concentration_array[17]
        c_HF2 = concentration_array[18]

        c_H = concentration_array[19]
        c_OH = concentration_array[20]


        # Defining total concentrations which describes conservation
        c_Zn_tot = c_Zn_2 + c_ZnOH4 + c_ZnOH3 + c_ZnOH2 + c_ZnOH + c_ZnO + c_ZnCO3 + c_ZnF_1
        c_COx_tot = c_CO2 + c_CO3_2 + c_HCO3_1 + c_H2CO3
        c_K_tot = c_K_1 #+ c_KF + 2*c_K2CO3 + c_KOH
        c_F_tot = c_F_1 + c_HF + c_HF2 + c_ZnF_1 #+ c_KF 

        ## Conservation equations
        equation_array = np.zeros(4)
        equation_array[0] = c_Zn_tot - self.c_Zn_tot   # Total concentration of Zn species
        equation_array[1] = c_COx_tot - self.c_COx_tot # Total concentration of COx species
        equation_array[2] = c_F_tot - self.c_F_tot     # Total concentration of F species
        equation_array[3] = c_K_tot - self.c_K_tot     # Total concentration of K species

        return equation_array
    
    def calculate_Zn_solution_concentrations(self):
        '''
        This part uses fsolve to solve the system of equilibrium equations based on initial guesses. 
        Using lambda x permits the use of x in the prior methods

        It appends the solutions to the different equilibria the empty concentration_matrix in the init.
        '''

        for i in range(len(self.pH_range)):
            c_Zn_2_0, c_CO3_2_0, c_K_1_0, c_F_1_0 = self.c_Zn_tot, self.c_COx_tot, self.c_K_tot, self.c_F_tot

            x0 = np.array([c_Zn_2_0, c_CO3_2_0, c_K_1_0, c_F_1_0])
            options = {'maxfev': 10000}
            x = fsolve(lambda x: self.conservation_Zn_solution(x, self.pH_range[i]), x0, **options)

            self.concentration_matrix[:,i] = self.distribute_Zn_solution_species(x, self.pH_range[i])

    def plot_Zn_species_distribution(self):
        '''
        This part plots the concentration distribution for all Zn-species
        '''

        # Zn-species
        plt.figure()
        plt.plot(self.pH_range, self.concentration_matrix[0, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[1, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[2, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[3, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[4, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[5, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[6, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[7, : ], linewidth = 3)
        plt.hlines(self.c_Zn_2, min(self.pH_range)-0.75, max(self.pH_range)+.75, 'k', '--')
        plt.xlim(min(self.pH_range)-0.75, max(self.pH_range)+0.75)
        plt.title('Zn - ion species')
        plt.xlabel('pH / [-]')
        plt.ylabel('Concentration / [M]')
        plt.legend(['Zn2+', 'Zn(OH)4^2-', 'Zn(OH)3^-', 'Zn(OH)2 (aq)', 'Zn(OH)^+', 'ZnO', 'ZnCO3', 'ZnF^+'])
        plt.show()

    def plot_COx_species_distribution(self):
        '''
        This part plots the concentration distribution for all COx species
        '''

        # Carbonates
        plt.figure()
        plt.plot(self.pH_range, self.concentration_matrix[8, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[9, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[10, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[11, : ], linewidth = 3)
        plt.hlines(self.c_K2CO3, min(self.pH_range)-0.75, max(self.pH_range)+0.75, 'k', '--')
        plt.xlim(min(self.pH_range)-0.75, max(self.pH_range)+0.75)
        plt.title('COx - species')
        plt.xlabel('pH / [-]')
        plt.ylabel('Concentration / [M]')
        plt.legend(['CO2', 'H2CO3', 'HCO3^-', 'CO3^2-'])
        plt.show()

    def plot_F_species_distribution(self):
        '''
        This part plots the concentration distribution for all F-species
        '''
        # Zn-species
        plt.figure()
        plt.plot(self.pH_range, self.concentration_matrix[16, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[17, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[18, : ], linewidth = 3)
        plt.plot(self.pH_range, self.concentration_matrix[7, : ], linewidth = 3)
        plt.hlines(self.c_KF, min(self.pH_range)-0.75, max(self.pH_range)+0.75, 'k', '--')
        plt.xlim(min(self.pH_range)-0.75, max(self.pH_range)+0.75)
        plt.title('F - ion species')
        plt.xlabel('pH / [-]')
        plt.ylabel('Concentration / [M]')
        plt.legend(['F^-', 'HF', 'HF^2-', 'ZnF^+'])
        plt.show()

# Initialize ChemicalEquilibrium
c_Zn_2_0 = 10**(-6)
c_KOH_0 = 7
c_K2CO3_0 = 1.4
c_KF_0 = 1.4

# Initialises and solving the system
initial_concentrations = np.array([c_Zn_2_0, c_KOH_0, c_K2CO3_0, c_KF_0])   # Initial concentrations
Zn_solution_system = Zn_solution(initial_concentrations)                    # Initialises the Zn-solution class
Zn_solution_system.calculate_Zn_solution_concentrations()                   # Calculates the concentration distributions
Zn_solution_system.plot_Zn_species_distribution()                           # Plots the Zn-species concentration distribution
Zn_solution_system.plot_COx_species_distribution()                          # Plots the COx-species concentration distribution
Zn_solution_system.plot_F_species_distribution()                            # Plots the F-species concentration distribution
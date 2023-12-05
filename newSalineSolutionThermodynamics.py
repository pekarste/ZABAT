# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:03:21 2023

@author: sidselh
"""

# -*- coding: utf-8 -*-
import numpy as np                                      # Importing numpy for using 
import matplotlib.pyplot as plt                         # Importing matplotlib.pyplot for plotting
from scipy.optimize import fsolve                       # Importing fsolve from scipy.optimze for numerical solving power


# Defining an object SalineSolution 
class SalineSolution:
    def __init__(self, NaCl, c_COXt, p_CO2, H_CO2):     # Initialising the class with input self, NaCl, c_COXt, p_CO2, H_CO2
        
        '''
        NaCl: Concentration of NaCl
        c_COXt: Total concentration of COx containing species
        p_CO2: Partial pressure of CO2
        H_CO2: Concentration of H_CO2?
        '''
        
        self.NaCl = NaCl                                # Initialise the NaCl concentration from the input "NaCl"
        self.c_COXt = c_COXt                            # Initialise the total concentration of COx containing species from the input "c_COXt"
        self.p_CO2 = p_CO2                              # Initialise the CO2 pressure from the input "p_CO2"
        self.H_CO2 = H_CO2                              # Initialise the HCO2 concentration from the input "H_CO2"
        self.pH_range = np.arange(0, 14.1, 0.1)         # Initialies an array for the pH = {0,14} with incerements of 0.1
        self.c_Clt = NaCl                               # Initialise the total concentration of Cl- ions by setting the total concentration equal to the input concentration of "NaCl"
        self.c_Nat = NaCl                               # initialise the total concentration of Na+ ions by setting the total concentration equal to the input concentration of "NaCl"
        self.c_CO2_sat = H_CO2 * p_CO2                  # Initialise the saturation concentration of CO2 by multiplying the input concentrations of "H_CO2" with "p_CO2"

        # Initialize concentration matrix
        self.c = np.zeros((9, len(self.pH_range)))      # Making an array of zeros with shape=(9, len(self.pH_range)) to store values. Makes a matrix with 9 rows and len(self.pH_range) columns

    # Method for calculating the concentration of different species as a function of pH.
    def distribute_saline_solution_species(self, x, pH):# Method for calculating the distribution of saline species
        
        '''
        x: Input variable ...?
        pH: This is the pH as pH = {0, 14} with incerements of 0.1
        '''
        
        c = np.zeros(9)                                 # c is an empty array with 9 entries

        c_H = 10**(-pH)                                 # Defining the concentration of H+ to be 10^(-pH)
        c_OH = 10**(-13.96) / c_H                       # Defining the concentration of OH- to be K_W/c_H

        c_Cl = x[0]                                     # First entry in x is the concentration oc Cl- ions
        c_Na = x[1]                                     # Second entry in x is the concentration of Na+ ions
        c_CO3 = x[2]                                    # Third entry in x is the concentration of CO3^2- ions

        ## Acids

        # Constants
        Ka_1 = 10**(-6.35)                              # Ka_1 for H2CO3 - Found in SI Chemical data
        Ka_2 = 10**(-10.33)                             # Ka_2 for H2CO3 (or Ka_1 for HCO3-) - Found in SI Chemical data
        K_3 = 1.7*10**(3)                               # Equilibrium constant for reaction ebtween CO2 and H2O - Don't found myself, but can be calculated from thermodynamic data
        Kcomplex = 10**(-0.2)                           # Kcomplex for the complexion of NaOH - I don't know where this is from

        # Species
        c_HCO3 = c_H * c_CO3 /Ka_2                      # Using law of mass action on concentration of HCO3-
        c_H2CO3 = c_H * c_HCO3 /Ka_1                    # Using law of mass action on concentration of H2CO3
        c_CO2 = c_H2CO3 * K_3                           # Using law of mass action on concentration oc CO2 in water

        # Na Complexes
        c_NaOH = c_Na * c_OH * Kcomplex                 # Law of mass action on concentration of NaOH

        # Initialise concentration array
        c[0] = c_Cl                                     # First entry in "c" is the concentration Cl' ions
        c[1] = c_Na                                     # Second entry in "c" is the concentration Na+ ions
        c[2] = c_CO3                                    # Third entry in "c" is the concentration of CO3^2- ions
        c[3] = c_HCO3                                   # Fourth entry in "c" is the concentration of HCO3^- ions
        c[4] = c_H2CO3                                  # Fifth entry in "c" is the concentration H2CO3
        c[5] = c_CO2                                    # Sixth entry in "c" is the concentration of CO2 in water
        c[6] = c_NaOH                                   # Seventh entry in "c" is the concentration of NaOH complex
        c[7] = c_H                                      # Eight entry in "c" is the proton concentration
        c[8] = c_OH                                     # Ninth entry in "c" is the hydroxide concentration

        return c                                        # Returns the concentration arra "c" 

    # Method for calculating the conservation of saline species as a function of pH
    def conserve_saline_solution_species(self, x, pH):
        c = self.distribute_saline_solution_species(x, pH)

        c_Cl = c[0]
        c_Na = c[1]
        c_CO3 = c[2]
        c_HCO3 = c[3]
        c_H2CO3 = c[4]
        c_CO2 = c[5]
        c_NaOH = c[6]
        c_H = c[7]
        c_OH = c[8]

        ClT = c_Cl
        COXT = c_CO2 + c_CO3 + c_HCO3 + c_H2CO3
        NaT = c_Na + c_NaOH

        results = np.zeros(3)
        results[0] = ClT - self.c_Clt
        results[1] = NaT - self.c_Nat
        results[2] = COXT - self.c_COXt

        return results

    def calculate_concentrations(self):
        for i in range(len(self.pH_range)):
            c_Cl0, c_Na0, c_CO30 = self.c_Clt, self.c_Nat, self.c_COXt

            x0 = np.array([c_Cl0, c_Na0, c_CO30])
            options = {'maxfev': 10000}
            x = fsolve(lambda x: self.conserve_saline_solution_species(x, self.pH_range[i]), x0, **options)

            self.c[:, i] = self.distribute_saline_solution_species(x, self.pH_range[i])

    def plot_species_distribution(self):
        ind = np.where((self.c[5, :] - self.c_CO2_sat) > 0)[0]
        sat_pH = self.pH_range[ind[-1]]
        y_sat = [np.min(self.c[5, :]), np.max(self.c[5, :])]

      
        # NaCl Species Distribution
        plt.figure()
        plt.plot(self.pH_range, self.c[0, :], linewidth=3)
        plt.plot(self.pH_range, self.c[1, :], linewidth=3)
        plt.plot(self.pH_range, self.c[6, :], linewidth=3)
        plt.title('NaCl Species Distribution')
        plt.xlabel('pH / -')
        plt.ylabel('Concentration / mol L$^{-1}$')
        plt.legend(['Cl', 'Na', 'NaOH'])
        plt.show()

        # CO_x Species Distribution
        plt.figure()
        plt.plot(self.pH_range, self.c[2, :], linewidth=3)
        plt.plot(self.pH_range, self.c[3, :], linewidth=3)
        plt.plot(self.pH_range, self.c[4, :], linewidth=3)
        plt.plot(self.pH_range, self.c[5, :], linewidth=3)
        plt.plot([sat_pH, sat_pH], y_sat, linewidth=3, color='k', linestyle='--')
        plt.title('CO_x Species Distribution')
        plt.xlabel('pH / -')
        plt.ylabel('Concentration / mol L$^{-1}$')
        plt.legend(['CO$_3$', 'HCO$_3$', 'H$_2$CO$_3$', 'CO$_2$'])
        plt.show()

# Example usage:
NaCl_value = 1.02
c_COXt_value = 0.002
p_CO2_value = 0.06
H_CO2_value = 0.031

saline_system = SalineSolution(NaCl_value, c_COXt_value, p_CO2_value, H_CO2_value)
saline_system.calculate_concentrations()
saline_system.plot_species_distribution()


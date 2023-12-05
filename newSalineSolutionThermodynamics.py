# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:03:21 2023

@author: sidselh
"""

# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


class SalineSolution:
    def __init__(self, NaCl, c_COXt, p_CO2, H_CO2):
        self.NaCl = NaCl
        self.c_COXt = c_COXt
        self.p_CO2 = p_CO2
        self.H_CO2 = H_CO2
        self.pH_range = np.arange(0, 14.1, 0.1)
        self.c_Clt = NaCl
        self.c_Nat = NaCl
        self.c_CO2_sat = H_CO2 * p_CO2

        # Initialize concentration matrix
        self.c = np.zeros((9, len(self.pH_range)))

    def distribute_saline_solution_species(self, x, pH):
        c = np.zeros(9)

        c_H = 10**(-pH)
        c_OH = 10**(-13.96) / c_H

        c_Cl = x[0]
        c_Na = x[1]
        c_CO3 = x[2]

        # Acids
        c_HCO3 = c_H * c_CO3 * 10**10
        c_H2CO3 = c_H * c_HCO3 * 10**6.16
        c_CO2 = c_H2CO3 * 1.7e3

        # Na Complexes
        c_NaOH = c_Na * c_OH * 10**-0.2

        c[0] = c_Cl
        c[1] = c_Na
        c[2] = c_CO3
        c[3] = c_HCO3
        c[4] = c_H2CO3
        c[5] = c_CO2
        c[6] = c_NaOH
        c[7] = c_H
        c[8] = c_OH

        return c

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

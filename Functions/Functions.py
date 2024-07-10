# Using van't Hoff to calculate the Gibbs free energy at a different temperature than the reference temperature (at which deltaH_1 is measured) which is T1
import numpy as np
def vant_Hoff(deltaG_1, deltaH_1, T1, T2):
    '''
    Assumptions:
    Assumes constant pressure

    Input:
    deltaG_1 : Gibbs free energy for reaction at T1 - [J/mol]
    deltaH_1 : Reaction enthalpy at T1              - [J/mol]
    T1 : The reference temperature                  - [K] - Set to be 25 degrees -- Can perhaps define it within the function and not have it as an input parameter
    T2 : The new temperature                        - [K]

    Output:
    DeltaG_2 : Gibbs free energy for reaction at T2 - [J/mol]
    '''
    return deltaG_1*(T2/T1) + deltaH_1*(1 - (T2/T1))

def deltaG_weak(deltaH_1, deltaS_1, T2):
    '''
    Assumptions:
    Assumes that both the enthalpy and entropy of the reaction are weak functions of temperature

    Input:
    deltaH_1 : The reaction enthalpy at T1              - [J/mol]
    deltaS_1 : The reaction entropy at T1               - [J/K*mol]
    T2 : The new temperature                            - [K]

    Output:
    DeltaG_2 : Gibbs free ennergy for reaction at T2    - [J/mol]
    '''

    return deltaH_1 - T2*deltaS_1

def deltaG_T2(deltaH_1, deltaS_1, delta_Cp, T1, T2):
    '''
    Assumptions:
    Assumes that the heat capacities described by REVISED POURBIAX is accurate enought to make a difference

    Input:
    deltaH_1 : The reaction enthalpy at T1              - [J/mol]
    deltaS_1 : The reaction entropy at T1               - [J/K*mol]
    T1 : The old temperature (25 degrees)               - [K]
    T2 : The new temperature                            - [K]
    
    Output:
    deltaG_2 : The Gibbs free energy for reaction at T2 - [J/mol]
    '''

    deltaH_2 = deltaH_1 + (delta_Cp[0]*(T2 - T1) + delta_Cp[1]*(T2**2 - T1**2) - delta_Cp[2]*(1/T2 - 1/T1))
    deltaS_2 = deltaS_1 + (delta_Cp[0]*np.log(T2/T1) + delta_Cp[1]*(T2 - T1) - (1/2)*delta_Cp[2]*(1/(T2**2) - 1/(T1**2)))
    return deltaH_2 - T2*deltaS_2

def E0_2(E0_1, deltaS_1, T1, T2, n):
    '''
    Assumptions:
    Assumes that the entropy is a weak function of temperature

    Input:
    E0_1 : The standard reduction potential at T1   - [V vs SHE]
    deltaS_1 : The reaction entropy at T1           - [J/K*mol]
    T1 : The old temperature                        - [K]
    T2 : The new temperature                        - [K]
    n : The number of electrons                     - [-]

    Output:
    E0_2 : The standard reduction potential at T2   - [V vs SHE]
    '''

    return E0_1 + (deltaS_1/(n*96485))*(T2 - T1)

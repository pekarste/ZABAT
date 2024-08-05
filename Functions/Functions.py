import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Using van't Hoff to calculate the Gibbs free energy for a reaction
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

# Using the approximation that the enthalpy and the entropy are weak functions of temperature on the temperature range
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

# Using no assumptions except that the heat capaicites are valid within the temperature range
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

# Calculating the standard reduction potential using that the derivative of the standrad reduction potential with respect to temperature is the entropy for the reaction
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

# Adding polygons to fill in the shapes in the Pourbaix diagram
def add_polygon(ax, vertices, colour, alpha, label, text_coord, text, text_rotation=0):
    '''
    Fills the specified domain with colour using matplotlib.patches.Polygon

    Inputs:
    ax: The Axes targed for the polygon
    vertices: The vertices for the figure (using coordinates from the plot)
    colour: The colour of the domain
    alpha: Transparency (0-1)
    label: Label for the region (is not shown, more for keeping track)
    text_coord: The text coordinates
    text: The actual text for the domain
    text_rotation: The rotation of the text (if any)

    Outputs:
    polygon: Filled polygon in the plot
    
    '''
    polygon = Polygon(vertices, closed=True, fill=True, color=colour, alpha=alpha)
    ax.add_patch(polygon)
    ax.text(*text_coord, text, fontsize=12, color='black', rotation = text_rotation)
    ax.plot([], [], color=colour, alpha=alpha, label=label)

def E_OER_HER(E_0, slope, pH):
    '''
    Calculates the lines for the OER and HER for the Pourbaix diagram using hte Nernst equation

    Inputs:
    E_0: The standard reduction potential for either the OER or HER
    slope: The slope for the line
    pH: The pH at where we want to calculate the potential

    Output:
    E: The electrochemical ptoential from the Nernst equation    
    '''
    return E_0 + slope*pH
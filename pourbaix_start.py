import matplotlib.pyplot as plt
import numpy as np

# Constants
constants = {
'R': 8.31451, # [J/(K*mol)],
'T': 298.16,  # [K],
'n': 2,       # Number of electrons transferred
'F': 96485    # [C/mol]
}

# E0 constants
constants_E0 = {
'EI_0': -0.76,      # [V] - Zn^2+(aq) --> Zn(s)
'EII_0': 0.2,       # [V] - Zn(OH)^+ --> Zn(s) 
'EIII_0': -0.4225,  # [V] - Zn(OH)2(s) --> Zn(s) -  assumed solid phase, not aq.
'EIV_0': -0.1,      # [V] - Zn(OH)3^- --> Zn(s)
'EV_0': 0.278,      # [V] - Zn(OH)4^2- --> Zn(s)
'EVI_0': 1,         # [V] - ZnO --> Zn(s)
'EHER_0': 0,        # [V] - HER
'EOER_0': 1.23      # [V] -  OER
}

# deltaG constants
constants_deltaG = {
'deltaG_VIII': -94.3*10**3, # [J/mol]
'deltaG_IX': -16.5*10**3,   # [J/mol]
'deltaG_X': -6.9*10**3      # [J/mol]
}

# Plotting

pH = np.arange(0, 17, 0.01)     # pH range
EHER = constants_E0['EHER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH   # HER
EOER = constants_E0['EOER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH   # OER
slope = -2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)
angle = np.arctan(slope)


plt.figure()
plt.plot(pH,EHER, '--')     # Line for the HER
plt.plot(pH,EOER, '--')     # Line for the OER
pZn_values = np.array([0, 2, 4])# Values we iterate through for the activity of dissolved Zn
pZn_values_dict = {'pZn': pZn_values} # Maybe make this a dictionary as Dict = {'pZn': pZn_values} Or something...

for i in range(len(pZn_values)):

    pZn_string = f'pZn =  {int(pZn_values_dict['pZn'][i])}' # The string printed for the different pZn values in the plot. 

    # Concentration constants
    constants_p = {             # Setting all values to be the same
    'pZn': pZn_values[i],       # pZn = -log(Zn^2+)
    'pZnOH': pZn_values[i],     # pZnOH = -log(Zn(OH)^+)
    'pZnOH2': pZn_values[i],    # pZnOH2 = -log(Zn(OH)2)
    'pZnOH3': pZn_values[i],    # pZnOH3 = -log(Zn(OH)3^-)
    'pZnOH4': pZn_values[i]     # pZnOH4 = -log(Zn(OH)4^-2)
    }

    # E calculations

    # Zn^2+ --> Zn(s)
    EI = np.ones(len(pH))*constants_E0['EI_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZn'])
    
    E_Zn_Zn2 = np.ones(len(pH))*EI
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

    # pH calculations

    # Zn^2+ --> Zn(OH)2
    pHVIII = 14 + constants_deltaG['deltaG_VIII']/(2*constants['R']*constants['T']*np.log(10)) + constants_p['pZn']/2 - 0*constants_p['pZnOH2']/2

    # Zn(OH)2 --> Zn(OH)3^-
    pHIX = 14 + constants_deltaG['deltaG_IX']/(constants['R']*constants['T']*np.log(10)) + constants_p['pZnOH2'] - constants_p['pZnOH3']

    # Zn(OH)3^- --> Zn(OH)4^2-
    pHX = 14 + constants_deltaG['deltaG_X']/(constants['R']*constants['T']*np.log(10)) + constants_p['pZnOH3'] - constants_p['pZnOH4']

    # Since we have to decide some lines in order to get a cleaner plot, we decide the pH intersections from the sloped lines
    index_I_III = np.where(EI <= EIII)[0][-1] if np.any(EI <= EIII) else None       # Finds the closest intersection between EI and EIII 
    index_III_IV = np.where(EIII <= EIV)[0][-1] if np.any(EIII <= EIV) else None    # Finds the closest intersection between EIII and EIV
    index_IV_V = np.where(EIV <= EV)[0][-1] if np.any(EIV <= EV) else None          # Finds the closest intersection between EIV and EV

    # For å fikse unøyaktige skjæringspunkter kan man enten sjekke om man kan indeksere gjennom pH, eller gjøre linjestørrelsen litt tykkere

    plt.plot(pH[:index_I_III], EI[:index_I_III], 'k', label='Zn$^{2+}$ - Zn')
    #plt.vlines(pHVIII, -1.5, 1.5, 'k', label='Zn$^{2+}$ - Zn')
    #plt.plot(pH, EII,'k', label='Zn(OH)$^+$ - Zn')
    plt.plot(pH[index_I_III: index_III_IV], EIII[index_I_III: index_III_IV], 'c', label='Zn(OH)$_2$ - Zn')
    plt.plot(pH[index_III_IV:index_IV_V], EIV[index_III_IV:index_IV_V], 'r', label='Zn(OH)$_3^-$ - Zn')
    plt.plot(pH[index_IV_V:], EV[index_IV_V:], 'm', label='Zn(OH)$_4^{2-}$ - Zn')
    plt.vlines(pH[index_I_III], EI[index_I_III], 1.5, 'c', label = 'pH = ')
    #plt.vlines(pHIX, -1.5, 1.5, 'k', label='Zn(OH)$_{2}$ - Zn(OH)$_3^-$')
    plt.vlines(pH[index_III_IV], EIII[index_III_IV], 1.5, 'r', label='test2')
    #plt.vlines(pHX, -1.5, 1.5, 'k', label='Zn(OH)$_{3}^-$ - Zn(OH)$_4^{2-}$')
    plt.vlines(pH[index_IV_V], EIV[index_IV_V], 1.5, 'm', label = 'test3')
    plt.text(pH[index_I_III+10], 1, pZn_string, color='k', rotation = 90)   # Adding pZn values to vertical lines for different pZn values-- must be fixed
    plt.text(pH[index_III_IV+10], 1, pZn_string, color='r', rotation = 90)  # Adding pZn values to vertical lines for different pZn values-- must be fixed
    plt.text(pH[index_IV_V+10], 1, pZn_string, color='m', rotation = 90)    # Adding pZn values to vertical lines for different pZn values-- must be fixed
#plt.legend()

plt.text(2, -1.25, 'Zn(s)', fontsize=12, color='black')
plt.text(2, 0.5, 'Zn$^{2+}$(aq)', fontsize=12, color='black')
plt.text(2, -0.32, 'HER', fontsize=12, color='black', rotation=-10)
plt.text(2, 0.9, 'OER', fontsize=12, color='black', rotation=-10)
plt.text(14, -0.55, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)
plt.text(12, -0.4, 'Zn(OH)$^{-}_{3}$(aq)', fontsize=12, color='black', rotation=90)
plt.text(8, 0.0, 'Zn(OH)$_2$(s)', fontsize=12, color='black', rotation=90)


plt.xlabel('pH')
plt.ylabel('Potential [V]')
plt.xlim(xmin=0, xmax=max(pH))  # Set the x-axis range
plt.ylim(ymin=-1.5, ymax=1.5)
plt.show()


'''
Comment:
We have used Erbe's notes, and the Pourbaix notes to fill in each other, while in reality there are differences that should be considered. I am convinced that by using one or the other we are bound to fix the issues
'''
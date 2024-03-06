import matplotlib.pyplot as plt
import numpy as np

# Constants
constants = {
'R': 8.31451, # [J/(K*mol)],
'T': 298.16, # [K],
'n': 2,
'F': 96485 # [C/mol]
}

# E0 constants
constants_E0 = {
'EI_0': -0.76, # [V]
'EII_0': 0.2, # [V]
'EIII_0': -0.4225, # [V], assumed solid phase, not aq.
'EIV_0': -0.1, # [V]
'EV_0': 0.278, # [V]
'EVI_0': 1, # [V]
'EHER_0': 0, # [V]
'EOER_0': 1.23 # [V]
}

# deltaG constants
constants_deltaG = {
'deltaG_VIII': -94.3*10**3, # [J/mol]
'deltaG_IX': -16.5*10**3, # [J/mol]
'deltaG_X': -6.9*10**3 # [J/mol]
}

# Plotting
# pH range
pH = np.arange(0, 16, 0.01)
EHER = constants_E0['EHER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH
EOER = constants_E0['EOER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH

plt.figure()
plt.plot(pH,EHER, '--')
plt.plot(pH,EOER, '--')
pZn_values = np.array([0, 2, 4])

for i in range(len(pZn_values)):

    # Concentration constants
    constants_p = {
    'pZn': pZn_values[i],
    'pZnOH': pZn_values[i],
    'pZnOH2': pZn_values[i],
    'pZnOH3': pZn_values[i],
    'pZnOH4': pZn_values[i]
    }

    # E calculations

    # Zn^2+ --> Zn(s)
    EI = constants_E0['EI_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZn'])
    E_Zn_Zn2 = np.ones(len(pH))*EI
    # Zn(OH)^+ --> Zn(s)
    EII = constants_E0['EII_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH']+pH)

    # Zn(OH)2 --> Zn(s)
    EIII = constants_E0['EIII_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(0*constants_p['pZnOH2']+2*pH)

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

    index_I_III = np.where(E_Zn_Zn2 <= EIII)[0][-1] if np.any(E_Zn_Zn2 <= EIII) else None
    index_III_IV = np.where(EIII <= EIV)[0][-1] if np.any(EIII <= EIV) else None
    index_IV_V = np.where(EIV <= EV)[0][-1] if np.any(EIV <= EV) else None

    # For å fikse unøyaktige skjæringspunkter kan man enten sjekke om man kan indeksere gjennom pH, eller gjøre linjestørrelsen litt tykkere

    plt.plot(pH[:index_I_III], E_Zn_Zn2[:index_I_III], 'k', label='Zn$^{2+}$ - Zn')
    #plt.plot(pH, E_Zn_Zn2, 'g', label='Zn$^{2+}$ - Zn')
    #plt.plot(pH, EII,'k', label='Zn(OH)$^+$ - Zn')
    plt.plot(pH[index_I_III: index_III_IV], EIII[index_I_III: index_III_IV], 'c', label='Zn(OH)$_2$ - Zn')
    plt.plot(pH[index_III_IV:index_IV_V], EIV[index_III_IV:index_IV_V], 'r', label='Zn(OH)$_3^-$ - Zn')
    plt.plot(pH[index_IV_V:], EV[index_IV_V:], 'm', label='Zn(OH)$_4^{2-}$ - Zn')
    #plt.plot(pH, EIII, 'k', label='Zn(OH)$_2$ - Zn')
    #plt.plot(pH, EIV, 'r', label='Zn(OH)$_3^-$ - Zn')
    #plt.plot(pH, EV, 'b', label='Zn(OH)$_4^{2-}$ - Zn')
    #plt.plot(pH, EVI, 'k', label='ZnO - Zn')
    #plt.vlines(pHVIII, -1.5, 1.5, 'k', label='Zn$^{2+}$ - Zn(OH)$_2$')
    plt.vlines(pH[index_I_III], E_Zn_Zn2[index_I_III], 1.5, 'k', label = 'test')
    #plt.vlines(pHIX, -1.5, 1.5, 'k', label='Zn(OH)$_{2}$ - Zn(OH)$_3^-$')
    plt.vlines(pH[index_III_IV], EIII[index_III_IV], 1.5, 'r', label='test2')
    #plt.vlines(pHX, -1.5, 1.5, 'k', label='Zn(OH)$_{3}^-$ - Zn(OH)$_4^{2-}$')
    plt.vlines(pH[index_IV_V], EIV[index_IV_V], 1.5, 'm', label = 'test3')
#plt.legend()
plt.text(2, -1.25, 'Zn(s)', fontsize=12, color='black')
plt.text(2, 0.5, 'Zn$^{2+}$(aq)', fontsize=12, color='black')
plt.text(2, -0.40, 'HER', fontsize=12, color='black', rotation=-30)
plt.text(2, 0.75, 'OER', fontsize=12, color='black', rotation=-30)
plt.text(14, -0.55, 'Zn(OH)$^{2-}_{4}$(aq)', fontsize=12, color='black', rotation=90)
plt.text(12, -0.4, 'Zn(OH)$^{-}_{3}$(aq)', fontsize=12, color='black', rotation=90)
plt.text(8, 0.0, 'Zn(OH)$_2$(s)', fontsize=12, color='black', rotation=90)
plt.xlabel('pH')
plt.ylabel('Potential [V]')
plt.xlim(xmin=0, xmax=max(pH))  # Set the x-axis range
plt.ylim(ymin=-1.5, ymax=1.5)
plt.show()
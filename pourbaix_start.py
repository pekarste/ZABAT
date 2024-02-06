import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import scipy

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

# Concentration constants
constants_p = {
'pZn': 0,
'pZnOH': 0,
'pZnOH2': 0,
'pZnOH3': 0,
'pZnOH4': 0
}

# deltaG constants
constants_deltaG = {
'deltaG_VIII': -94.3*10**3, # [J/mol]
'deltaG_IX': 16.5*10**3, # [J/mol]
'deltaG_X': -6.9*10**3 # [J/mol]
}

# pH range
pH = np.arange(0, 15, 0.1)

# E calculations
EI = constants_E0['EI_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZn'])
EII = constants_E0['EII_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH']+pH)
EIII = constants_E0['EIII_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH2']+2*pH)
EIV = constants_E0['EIV_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH3']+3*pH)
EV = constants_E0['EV_0'] - constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(constants_p['pZnOH4']+4*pH)
EVI = constants_E0['EVI_0'] - 2*constants['R']*constants['T']*np.log(10)/(constants['n']*constants['F'])*(pH)
EHER = constants_E0['EHER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH
EOER = constants_E0['EOER_0'] - 2*constants['R']*constants['T']/(constants['n']*constants['F'])*np.log(10)*pH

# pH calculations
pHVIII = 14 + constants_deltaG['deltaG_VIII']/(2*constants['R']*constants['T']*np.log(10)) + constants_p['pZn']/2 - constants_p['pZnOH2']/2
pHIX = 14 + constants_deltaG['deltaG_IX']/(constants['R']*constants['T']*np.log(10)) + constants_p['pZnOH2'] - constants_p['pZnOH3']
pHX = 14 + constants_deltaG['deltaG_X']/(constants['R']*constants['T']*np.log(10)) + constants_p['pZnOH3'] - constants_p['pZnOH4']

# print(pHVIII, pHIX, pHX)
# print(EHER)

# Plotting
plt.figure()
plt.plot(pH,EHER, '--')
plt.plot(pH,EOER, '--')
plt.hlines(constants_E0['EI_0'], pH[0], pH[-1], label='Zn$^{2+}$ - Zn')
plt.plot(pH, EII,'k', label='Zn(OH)$^+$ - Zn')
plt.plot(pH, EIII, 'k', label='Zn(OH)$_2$ - Zn')
plt.plot(pH, EIV, 'k', label='Zn(OH)$_3^-$ - Zn')
plt.plot(pH, EV, 'k', label='Zn(OH)$_4^{2-}$ - Zn')
plt.plot(pH, EVI, 'k', label='ZnO - Zn')
plt.vlines(pHVIII, -1.5, 1.5, 'k', label='Zn$^{2+}$ - Zn(OH)$_2$')
plt.vlines(pHIX, -1.5, 1.5, 'k', label='Zn(OH)$_{2}$ - Zn(OH)$_3^-$')
plt.vlines(pHX, -1.5, 1.5, 'k', label='Zn(OH)$_{3}^-$ - Zn(OH)$_4^{2-}$')
plt.legend()
plt.xlabel('pH')
plt.ylabel('Potential [V]')
plt.xlim(xmin=0, xmax=17)  # Set the x-axis range
plt.ylim(ymin=-1.5, ymax=1.5)
plt.show()

# test = 1.23+(1/F)*(-237*10**3+157*10**3)
# print(test)


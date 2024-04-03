import numpy as np
R = 8.314
T = 25+273.15
F = 96485

G_Zn_2plus = -35814 # cal
G_ZnOH_plus = -78700 # cal

G_VIII = -94.3*10**3# J*mol^-1

pKa_VIII = 14 + G_VIII/(2*8.314*298*np.log(10))
print((2*R*T)*np.log(10)/(2*F))
print((4*R*T)*np.log(10)/(2*F))
#print(pKa_VIII)
#print(10.96/2)
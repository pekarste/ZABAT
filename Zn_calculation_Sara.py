from scipy.optimize import fsolve
import math

def equations(p):
    x, y = p
    Ksp = 10**(-15.66)
    K = 10**17.66
    eq1 = (x-y)*(6+2*x-4*y)**2 - Ksp
    eq2 = (y)/((x-y)*(6+2*x-4*y)**4) - K
    return [eq1, eq2]

x, y =  fsolve(equations, (0.001, 0.001))

print(10**(17.66)*10**(-15.66)*6**2)
print(10**(-15.66)/6**2)

# importing the required modules
import math
import matplotlib.pyplot as plt
import numpy as np
import sympy

#from sympy import symbols

#from sympy import *

en = float(input("photon energy (eV) = "))
#cff = float(input("constant focal distance = "))
m = float(input("diffraction order = "))
R = float(input("distance object - grating (mm) = "))
R_prime = float(input("distance grating - image (mm) = "))
a0 = float(input("groove density (l/mm) = "))
#a1 = float(input("groove density a1(l/mm^2)"))
#a2 = float(input("groove density a2(l/mm^3)"))
#a3 = float(input("groove density a3(l/mm^4)"))

l_mm = 1240*1e-6/en

A = m*l_mm*a0*R_prime
B = R+R_prime
delta = (A**2+B**2-A**2*B/R_prime)**0.5

x1 = (-A+delta)/B
y1 = A/R_prime+x1
x2 = (-A-delta)/B
y2 = A/R_prime+x2

print('x1 =', x1, 'x2 = ', x2, 'y1 = ', y1, 'y2 = ', y2)

if (math.fabs(x1)>1):
    print(" The solution is outside the ArcSin domain ")
else:
    print("\n")
    alpha = math.asin(x1)
    #print(Convert_to_deg(alpha))
    #print(m * Lambda * D0,  first_solution)
    y1 = A/R_prime+x1
    y2 = A/R_prime+x2


    beta1 = math.asin(y1)

    #image size is calculated straight into microns:

#print(A, B, l_m, x1, x2)
alpha_deg = alpha/math.pi*180

beta_deg = beta/math.pi*180


print('alpha (rad) =',alpha,', alpha (deg) =', alpha_deg,)
print('beta (rad) =',beta,', beta (deg) =', beta_deg,)


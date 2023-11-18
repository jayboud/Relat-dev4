import numpy as np
import sympy as sp
import scipy
sp.init_printing()


r,rs,p2,a,d,theta = sp.symbols('r,r_s,p2,a,d,theta')


# changements de variable et simplifications
rs = 0
r = a*sp.sin(theta)

# variables intermediaires
p2 = r*r + a*a*sp.cos(theta)**2
d = r*r + a*a - r*rs

# tenseur metrique
g_00 = (1 - r*rs/p2)
g_03 = 2*a*r*rs*sp.sin(theta)**2/p2
g_11 = -p2/d
g_22 = -p2/d
g_33 = -(r*r + a*a + a*a*r*rs*sp.sin(theta)**2/p2)*sp.sin(theta)**2
g = [g_00,g_03,g_11,g_22,g_33]

for element in g:
    sp.simplify(element)

print(g)

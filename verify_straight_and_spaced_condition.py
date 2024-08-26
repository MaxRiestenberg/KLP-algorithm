# Verify the straight-and-spaced criteria
# Author: Max Riestenberg

import numpy as np
# import matplotlib.pyplot as plt
import time

# Load the saved files and read off the numbers
# Note that if you just ran "compute_straight_and_spaced_parameters" the file name has a _temp at the end, so you need to change it.
canglesP=np.load('cosangleP_pairs.npy')
canglesM=np.load('cosangleM_pairs.npy')
spacing=np.load('spacing_pairs.npy') 

print('The smallest cosine of a zeta-zeta angle is')
print(np.min(canglesP))
print('The smallest cosine of a iotazeta-iotazeta angle is')
print(np.min(canglesM))
print('The largest iotazeta-zeta angle between midpoints is at most')
print(np.arccos(np.min(canglesP))+np.arccos(np.min(canglesM)))
print('The smallest d_alpha spacing parameter is')
print(np.min(spacing))

# Now we check the conditions of Theorem 4.2.
dim=3; 
zeta0=np.sqrt(1/(2*dim*(dim-1))); # a constant depending on dimension
print(f'zeta0 = {zeta0}')
parameter=0.7; # set a parameter strictly between 0.5 and 1 to determine the auxiliary paramter epsilon_aux

epsmax=np.arccos(-1/(dim-1)); 
eps=np.arccos(np.min(canglesP))+np.arccos(np.min(canglesM)); 
eps=np.arccos(0.87)+np.arccos(0.86)
epsaux=parameter*eps+(1-parameter)*epsmax; 
print(f'epsaux = {epsaux}')
s=np.min(spacing); 

delta1=np.arccosh(np.sqrt(dim/((1+(dim-1)*np.cos(epsaux)))))+10.0**(-10); 
delta2=np.arccosh(np.sqrt(dim/((1+(dim-1)*np.cos(epsaux-eps)))))+10.0**(-10); 
delta3=np.arccosh(np.sqrt(dim/((1+(dim-1)*np.cos(2*epsaux-eps)))))+10.0**(-10); 
delta4=delta3+np.expm1(delta3)*np.exp(-s)+10.0**(-10); 

print(f'delta1={delta1}')
print(f'delta2={delta2}')
print(f'delta3={delta3}')
print(f'delta4={delta4}')

print(f'Item 1a: {np.cos(2*epsaux-eps) > -1/(dim-1)}')
print(f'Item 1b: {np.cos(epsaux+(delta4*zeta0)/np.sinh(s-delta4))> -1/(dim-1)}')
print(f'Item 2a: {(1-dim)*np.cos(epsaux)+dim/(np.cosh(delta1))**2 <= 1}')
print(f'Item2b: {(1-dim)*np.cos(2*epsaux-eps)+dim/(np.cosh(delta3))**2 <= 1}')
print(f'Item3: {(1-dim)*np.cos(epsaux-eps)+dim/(np.cosh(delta2))**2 <= 1}')
print(f'Item4: {np.exp(delta1-s)-np.exp(-s)<= delta2}')
print(f'Item5: {2*delta4<s}')

print('If all items are true, then the example passes the check!')
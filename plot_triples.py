# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 13:11:50 2021

@author: Max Riestenberg
"""

## Plot pairs ##

import numpy as np
import matplotlib.pyplot as plt

# np.save('cosangle_pairs.npy', cangles) # save
# np.save('spacing_pairs.npy', spacing) # save
cangles = np.load('cosangle_triples.npy') # load
spacing = np.load('spacing_triples.npy') # load
print(spacing)

## Plot the result ##
# We will make a plot with three horizontal lines.
# The horizontal lines indicate whether we have passed the straight and spaced check.

numbery=np.cos(2*np.pi/3)
numberr=np.cos(np.pi-(1/40))
numberg=(4/(1+spacing))-1
print(numberg)

plt.figure(figsize=(5,4))
plt.plot(cangles,color='gray')
plt.title("$ \cos ({\\angle}_{m_2} (m_1,m_3))$")
plt.xlabel("Words of length 9")
plt.ylim([-1.1,0])
plt.axhline(y=numberg, color='black', linestyle=':') 
plt.axhline(y=numberr, color='black', linestyle='-')
plt.axhline(y=numbery, color='black', linestyle='--')
plt.savefig('words9.eps', format='eps')
plt.show()

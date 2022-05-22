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
cangles = np.load('cosangle_pairs.npy') # load
spacing = np.load('spacing_pairs.npy') # load

## Plot the result ##
# We will make a plot with three horizontal lines.
# The horizontal lines indicate whether we have passed the straight and spaced check.

# numbery is a hypothetical bound from eps(zetamod)
numbery=np.cos(np.pi/6)
# numberr is the bound we need to apply my l2g
numberr=np.cos(1/80)
# numberg is the bound we need to apply the l2g from hyperbolic space
numberg=spacing

# Make the plot

plt.figure(figsize=(5,4))
plt.plot(cangles,color='gray')
plt.title("$ \cos ({\\angle}_{m_1} (p,m_2))$")
plt.xlabel("Words of length 8")
plt.ylim([0,1.1])
plt.axhline(y=numberg, color='black', linestyle=':')
plt.axhline(y=numberr, color='black', linestyle='-')
plt.axhline(y=numbery, color='black', linestyle='--')
plt.savefig('words8.eps', format='eps')
plt.show()
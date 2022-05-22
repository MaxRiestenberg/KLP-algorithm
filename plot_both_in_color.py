# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 12:51:43 2021

@author: Max Riestenberg
"""

## Plot pairs and triples in color for my website ##
## The code is essentially done, but I need to run it on the data of triples. ##

import numpy as np
import matplotlib.pyplot as plt

# np.save('cosangle_pairs.npy', cangles) # save
# np.save('spacing_pairs.npy', spacing) # save
cangles_p = np.load('cosangle_pairs.npy') # load
spacing_p = np.load('spacing_pairs.npy') # load

cangles_t = np.load('cosangle_triples.npy') # load
spacing_t = np.load('spacing_triples.npy') # load

## Plot the result ##
# We will make a plot with three horizontal lines.
# The horizontal lines indicate whether we have passed the straight and spaced check.

# numbery is a hypothetical bound from eps(zetamod)
numbery_p=np.cos(np.pi/6)
# numberr is the bound we need to apply my l2g
numberr_p=np.cos(1/80)
# numberg is the bound we need to apply the l2g from hyperbolic space
numberg_p=spacing_p

numbery_t=np.cos(2*np.pi/3)
numberr_t=np.cos(np.pi-(1/40))
numberg_t=(4/(1+spacing_t))-1


# Make the plot

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(13,4))
fig.suptitle('Checking the straight and spaced condition')

ax1.plot(cangles_p)
ax1.set_title("$ \cos ({\\angle}_{m_1} (p,m_2))$")
ax1.set_xlabel("Words of length 8")
ax1.set_ylim([0,1.1])
ax1.axhline(y=numberg_p, color='g', linestyle='-')
ax1.axhline(y=numberr_p, color='r', linestyle='-')
ax1.axhline(y=numbery_p, color='y', linestyle='-')

ax2.plot(cangles_t)
ax2.set_title("$ \cos ({\\angle}_{m_2} (m_1,m_3))$")
ax2.set_xlabel("Words of length 9")
ax2.set_ylim([-1.1,0])
ax2.axhline(y=numberg_t, color='g', linestyle='-') 
ax2.axhline(y=numberr_t, color='r', linestyle='-')
ax2.axhline(y=numbery_t, color='y', linestyle='-')

plt.show()
plt.savefig('words8and9.eps', format='eps')
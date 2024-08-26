#!/usr/bin/env python
# coding: utf-8

# @author: Max Riestenberg

# Check the straight and spaced condition for a surface group of genus 2 acting on SL(3,R)/SO(3).
# Step 1: Define generators "a,b,c,d" of the surface group in SO(2,1).
# Step 2: Choose a natural number "k". Enumerate words of length 2k.
# Step 3: For a word w of length 2k, decompose: w=w1*w2
#          -Choose a basepoint p=identity coset. 
#          -Measure the distance from m1=mid(p,w1^(-1)*p) to m2=mid(p,w2*p)
#          -Measure the zeta-angle and iota-zeta-angle based at m1 between p and m2.
#          -The worst angles between midpoints are the sum the worst measurements from the previous step.
# Step 4: Check the straight and spaced condition. (Do this with "verify_straight_and_spaced_condition.py")
# Conclusion: In this example, the check passes for words of length 2k=8. 

# import numpy and plt
import numpy as np
# import matplotlib.pyplot as plt
import time

## STEP 1 ##

# define the first generator a and its inverse A
# "cot(pi/8) = cosh(t/2)"
t=2*np.arccosh(1/np.tan(np.pi/8))

c=np.cosh(t)
s=np.sinh(t)

a=np.array([[c,0,s],[0,1,0],[s,0,c]])
A=np.array([[c,0,-s],[0,1,0],[-s,0,c]])

# define b and B
theta=np.pi/4
cc=np.cos(theta)
ss=np.sin(theta)
rot=np.array([[cc,ss,0],[-ss,cc,0],[0,0,1]])
rotinv=np.array([[cc,-ss,0],[ss,cc,0],[0,0,1]])
b=rot@a@rotinv
B=rot@A@rotinv

# define c and C
theta=np.pi/2
cc=np.cos(theta)
ss=np.sin(theta)
rot=np.array([[cc,ss,0],[-ss,cc,0],[0,0,1]])
rotinv=np.array([[cc,-ss,0],[ss,cc,0],[0,0,1]])
c=rot@a@rotinv
C=rot@A@rotinv

# define d and D
theta=3*(np.pi/4)
cc=np.cos(theta)
ss=np.sin(theta)
rot=np.array([[cc,ss,0],[-ss,cc,0],[0,0,1]])
rotinv=np.array([[cc,-ss,0],[ss,cc,0],[0,0,1]])
d=rot@a@rotinv
D=rot@A@rotinv

# Optional: check the relation. the check passes!
# print(a@d@C@b@A@D@c@B)

## Define important geometric functions for the computation.
def dist(mat1,mat2):
    """compute the Riemannian distance between two symmetric positive definite matrices"""
    """ Here we take the Riemannian metric induced by twice the trace form. This matches my convention in the paper"""
    a1, k1 = np.linalg.eigh(mat1)
    g1 = k1@np.diag(a1**(-1/2))@np.transpose(k1) # inverse square root of mat1. g1 translates mat1 to the identity.
    x1 = g1@mat2@np.transpose(g1) # We now apply g1 to mat2, and measure the distance of the result to the identity
    vvd, garbage = np.linalg.eigh(x1) # compute via eigh
    vvd = np.abs(vvd)
    vvd = np.log(vvd)/2 #strictly speaking, the vvd is half of the logs of singular values, because of the orbit map...
    return np.sqrt(2)*np.linalg.norm(vvd)

def dalpha(mat1,mat2):
    """compute the pseudometric d_alpha"""
    """Here we use singular values of x^-1/2 y x^-1/2"""
    a1, k1 = np.linalg.eigh(mat1)
    g1 = k1@np.diag(a1**(-1/2))@np.transpose(k1) # inverse square root of mat1
    x1 = g1@mat2@np.transpose(g1) # g1 translates mat1 to the identity. We now apply it to mat2, and measure the distance of the result to the identity
    vvd, garbage = np.linalg.eigh(x1) # compute via eigh
    return (np.log(vvd[-1])-np.log(vvd[-2]))/2

def zangle(mat1,mat2,z):
    """compute zeta-angle at basepoint between p and q"""
    garbage,k1=np.linalg.eigh(mat1)
    garbage,k2=np.linalg.eigh(mat2)
    return np.trace(k1@z@np.transpose(k1)@k2@z@np.transpose(k2))/np.trace(z@z)

# Make a dictionary
test_vars = {}
test_vars['a'] = locals()['a']
test_vars['b'] = locals()['b']
test_vars['c'] = locals()['c']
test_vars['d'] = locals()['d']
test_vars['A'] = locals()['A']
test_vars['B'] = locals()['B']
test_vars['C'] = locals()['C']
test_vars['D'] = locals()['D']

# Using kbmag in GAP, we obtain a text file of all geodesic words of length at most 9 for this presentation.
# Open and read the text file of all geodesic words of length at most 9. 
# Opening this file takes me about 20 seconds.
# The path on my tablet is C:\\Users\\riest\\OneDrive\\Desktop\\length_9_wds.txt 
# Note: you need to change this path by hand to run this code!
filename = 'C:\\Users\\riest\\OneDrive\\Desktop\\length_9_wds.txt'
print(f'Opening file {filename}')
with open(filename) as f:
    lines=f.readlines()
print('I opened it!')

# A few helpful landmarks in the textfile.
first_length6=22289;
first_length7=155577;
first_length8=1085905;
first_length9=7579441;
last_lineplus1=52903257;

# In this for loop we iterate over the words of length 8 and compute the cosines of angles between midpoints.
# Running this for words of length 8 takes about 35 minutes on my tablet.

# First we set the bounds of the for loop. We will let i be in range(imin,imax). It does imin,...,imax-1 then stops.
imin=first_length8;
imax=first_length9;
length=imax-imin;

# Next we initialize the vector "cangles" with zeros. 
# The for loop below will compute the cosines of angles of midpoints and store them in cangles.
canglesP=np.zeros(length);
canglesM=np.zeros(length);
spacing=np.zeros(length);

# The zeta angle depends on a choice of diagonal matrix in positive Weyl chamber
# define zeta the diagonal matrix (d-1,-1,...,-1)
zeta=np.diag(np.array([2,-1,-1]));
iotazeta=np.diag(np.array([1,1,-2]));
# define trz2 = Tr(zeta.zeta) # Note that the scale of the metric won't matter for the way we use this later.
trz2=np.trace(zeta@zeta);

# find start time
start = time.time()
print(f"For loop start time: {start}")
print(f"For loop length: {length}")

# Now we actually do the computation. Could be optimized by re-using SVD computations.
for i in range(imin,imin+length):
    current_line=lines[i]
    swapped_line=current_line.swapcase()
    w1=test_vars[swapped_line[3]]@test_vars[swapped_line[2]]@test_vars[swapped_line[1]]@test_vars[swapped_line[0]]
    w2=test_vars[current_line[4]]@test_vars[current_line[5]]@test_vars[current_line[6]]@test_vars[current_line[7]]
    point1=w1@np.transpose(w1) # the translate of the basepoint by w1
    a1,k1 = np.linalg.eigh(point1) # eigh quickly computes eigenvalues of a real symmetric matrix 
    a1=np.abs(a1) # numerical issues may lead to negative eigenvalues
    mid1 = k1@np.diag(a1**(1/2))@np.transpose(k1) # the square root of hht
    mid1_inv_rt = k1@np.diag(a1**(-1/4))@np.transpose(k1) # this isometry translates mid1 to the origin
    aa2, kk2 = np.linalg.eigh(w2@np.transpose(w2)) # second point, svd
    aa2=np.abs(aa2) # numerical issues may lead to negative eigenvalues
    mid2 = kk2@np.diag(aa2**(1/2))@np.transpose(kk2) # square root of w2.Transpose(w2)
    spacing[i-imin]= dalpha(mid1,mid2)
    mid2_by_m1 = mid1_inv_rt@mid2@mid1_inv_rt # translate mid2 by the isometry taking mid1 to the origin. lets us measure angle at the origin
    canglesP[i-imin]=zangle(np.linalg.inv(mid1),mid2_by_m1,zeta)
    canglesM[i-imin]=zangle(np.linalg.inv(mid1),mid2_by_m1,iotazeta)

# calculate how long the for loop took
for_time = time.time()-start
print(f'For loop duration: {round(for_time)} seconds = {round(for_time/60)} minutes')

# If this is True, then some entry of cangles is nan, so we have a numerical issue.
if np.any(np.isnan(canglesP)):
    print('One of the angles is NAN')
if np.any(np.isnan(canglesM)):
    print('One of the angles is NAN')  

# Save the result
np.save('cosangleP_pairs_temp.npy', canglesP)
np.save('cosangleM_pairs_temp.npy', canglesM)
np.save('spacing_pairs_temp.npy', spacing)

# Print the result
print(f'The smallest cosine of a zeta-zeta angle is {np.min(canglesP)}')
print(f'The smallest cosine of a iotazeta-iotazeta angle is {np.min(canglesM)}')
print(f'The biggest iotazeta-zeta angle is at most {np.arccos(np.min(canglesP))+np.arccos(np.min(canglesM))}')
print(f'The smallest d_alpha spacing parameter is {np.min(spacing)}')

# Now use the other Python file to verify the criteria
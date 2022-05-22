#!/usr/bin/env python
# coding: utf-8

# @author: Max Riestenberg

# Check the straight and spaced condition for a surface group of genus 2 acting on the hyperbolic plane.
# Step 1: Define generators "a,b,c,d" of the surface group in SO(2,1).
# Step 2: Choose a natural number "k". Enumerate words of length 2k.
# Step 3: For a word w of length 2k, decompose: w=w1*w2
#          -Choose a basepoint "p" in hyperbolic space. 
#          -Measure the distance from m1=mid(p,w1^(-1)*p) to m2=mid(p,w2*p)
#          -Measure the angle based at m1 between p and m2.
#          -The angles between midpoints are half as straight as the worst measurement from the previous step.
# Step 4: Check the straight and spaced condition.
# Conclusion: 

# import numpy and plt
import numpy as np
import matplotlib.pyplot as plt
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

# define a basepoint.
p=np.array([[0],[0],[1]])

## Define important geometric functions for the computation. 
## We use the hyperboloid model of the hyperbolic plane. 

def form(vec1,vec2):
    """define a bilinear form on R3 with signature (2,1)"""
    return vec1[0,0]*vec2[0,0]+vec1[1,0]*vec2[1,0]-vec1[2,0]*vec2[2,0]

def norm(vector):
    """return the norm squared of a vector with respect to chosen bilinear form of signature (2,1)"""
    return form(vector,vector)

def dist(v,w):
    """return the hyperbolic distance from v to w"""
    """here, the metric is normalized to that we have half the frobenius norm on \mathfrak{p}"""
    return np.arccosh(-form(v,w))

def midp(v,w):
    """return the midpoint on H2 from v to w"""
    foo=(v+w)/2
    return foo/np.sqrt(-norm(foo))

# def angle(p,q,r):
#     """return the angle based at p from q to r"""
#     vec1= q+form(p,q)*p
#     vec1=vec1/(np.sqrt(norm(vec1)))
#     vec2= r+form(p,r)*p
#     vec2=vec2/(np.sqrt(norm(vec2)))
#     return np.arccos(form(vec1,vec2))

def cosangle(p,q,r):
    """return cosine of an angle based at p"""
    qq=q/form(q,p)
    rr=r/form(r,p)
    return form(qq+p,rr+p)/np.sqrt(norm(qq+p)*norm(rr+p))

def cosangle2(p,q,r):
    """return cosine of an angle based at p. Computed in an alternate way. Hoping it helps with numerical issues."""
    fqp=form(q,p)
    frp=form(r,p)
    vec1=q+fqp*p
    vec2=r+frp*p
    return form(vec1,vec2)/np.sqrt(norm(vec1)*norm(vec2))

def sscheck(g,h):
    """Check the straight and spaced condition at p,gp,hp"""
    m1=midp(p,g@p)
    m2=midp(p,h@p)
    spacing_factor=np.sqrt((1-form(m1,m2))/2)
    coss_half_angle=cosangle(m1,p,m2)
    return spacing_factor*coss_half_angle

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


# Open and read the text file of all geodesic words of length at most 9. Thanks Teddy and Wendy!
# Opening this file takes about 20 seconds
# The path on my tablet is C:\\Users\\riest\\OneDrive\\Desktop\\length_9_wds.txt
with open('C:\\Users\\John Maxwell\\Desktop\\length_9_wds.txt') as f:
    lines=f.readlines()
    

# A few helpful landmarks in the textfile.
first_length6=22289;
first_length7=155577;
first_length8=1085905;
first_length9=7579441;
last_lineplus1=52903257;

# In this for loop we iterate over the words of length 8 and compute the cosines of angles between midpoints.
# Running this for words of length 8 took 373 seconds on my desktop
# Running this for words of length 9 took 3786 seconds on my desktop

# First we set the bounds of the for loop. We will let i be in range(imin,imax), which stops at imax-1.
imin=first_length8;
imax=first_length9;

# Next we initialize the vector "cangles" with zeros. 
# The for loop below will compute the cosines of angles of midpoints and store them in cangles.
cangles=np.zeros(imax-imin)
spacing=0;

# find start time
start = time.time()

# Now we actually do the computation
for i in range(imin,imax):
    current_line=lines[i]
    swapped_line=current_line.swapcase()
    w1=test_vars[swapped_line[3]]@test_vars[swapped_line[2]]@test_vars[swapped_line[1]]@test_vars[swapped_line[0]]
    w2=test_vars[current_line[4]]@test_vars[current_line[5]]@test_vars[current_line[6]]@test_vars[current_line[7]]
    m1=midp(p,w1@p)
    m2=midp(p,w2@p)
    spacing=max(spacing,((1-form(m1,m2))/2)**(-1/2) )
    cangles[i-imin]=cosangle2(m1,p,m2)

# calculate how long the for loop took
for_time = time.time()-start
print('The for loop took ')
print(for_time)
print('seconds')

# If this is True, then some entry of cangles is nan, so we have a numerical issue.

if np.any(np.isnan(cangles)):
    print('One of the angles is NAN')

print('The spacing parameter s_min is ')
print(spacing)

np.save('cosangle_pairs.npy', cangles) # save
np.save('spacing_pairs.npy', spacing) # save
# new_num_arr = np.load('data.npy') # load

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

plt.plot(cangles)
plt.title("$ \cos ({\\angle}_{m_1} (p,m_2))$")
plt.xlabel("words of length 8")
plt.ylim([0,1.1])
plt.axhline(y=numberg, color='g', linestyle='-')
plt.axhline(y=numberr, color='r', linestyle='-')
plt.axhline(y=numbery, color='y', linestyle='-')
plt.show()

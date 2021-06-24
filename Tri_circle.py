import numpy as np


def dir_vec(A,B):
  return B-A

def norm_vec(A,B):
  return np.matmul(omat, dir_vec(A,B))

#Area using Hero's formula
def tri_hero(a,b,c):
  s = (a+b+c)/2
  area = np.sqrt(s*(s-a)*(s-b)*(s-c))
  return area
def tri_section(A,B,k):
  V = (k*A+B)/(k+1)
  return V

#Generate line points
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB

#Triangle vertices
def tri_vert(a,b,c):
  p = (a**2 + c**2-b**2 )/(2*a)
  q = np.sqrt(c**2-p**2)
  A = np.array([p,q]) 
  B = np.array([0,0]) 
  C = np.array([a,0]) 
  return A,B,C

#Generate line points
#def line_dir_pt(m,A,k1,k2):
# len =10
# x_AB = np.zeros((2,len))
# lam_1 = np.linspace(k1,k2,len)
# for i in range(len):
# temp1 = A + lam_1[i]*m
# x_AB[:,i]= temp1.T
# return x_AB

def line_dir_pt(m,A, dim):
  len = 10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,10,len)
  for i in range(len):
    temp1 = A + lam_1[i]*m
    x_AB[:,i]= temp1.T
  return x_AB


#Generate line points
#def line_gen(A,B):
# len =10
# x_AB = np.zeros((2,len))
# lam_1 = np.linspace(0,1,len)
# for i in range(len):
# temp1 = A + lam_1[i]*(B-A)
# x_AB[:,i]= temp1.T
# return x_AB

#Foot of the Altitude
def alt_foot(A,B,C):
  m = B-C
  n = np.matmul(omat,m) 
  N=np.vstack((m,n))
  p = np.zeros(2)
  p[0] = m@A 
  p[1] = n@B
  #Intersection
  P=np.linalg.inv(N.T)@p
  return P

#Radius and centre of the circumcircle
#of triangle ABC
def ccircle(A,B,C):
  p = np.zeros(2)
  n1 = dir_vec(B,A)
  p[0] = 0.5*(np.linalg.norm(A)**2-np.linalg.norm(B)**2)
  n2 = dir_vec(C,B)
  p[1] = 0.5*(np.linalg.norm(B)**2-np.linalg.norm(C)**2)
  #Intersection
  N=np.vstack((n1,n2))
  O=np.linalg.inv(N)@p
  r = np.linalg.norm(A -O)
  return O,r

#Inradius of triangle ABC
def tri_iradius(a,b,c):
  s = (a+b+c)/2
  area = tri_hero(a,b,c)
  r = area/s
  return r

#Circumradius of triangle ABC
def tri_cradius(a,b,c):
  s = (a+b+c)/2
  area = tri_hero(a,b,c)
  R = a*b*c/(4*area)
  return R

#Circumcentre of triangle ABC
def tri_ccentre(A,B,C):
  p = np.zeros(2)
  n1 = dir_vec(B,A)
  p[0] = 0.5*(np.linalg.norm(A)**2-np.linalg.norm(B)**2)
  n2 = dir_vec(C,B)
  p[1] = 0.5*(np.linalg.norm(B)**2-np.linalg.norm(C)**2)
  #Intersection
  N=np.vstack((n1,n2))
  O=np.linalg.inv(N)@p
  return O

#vertices of the intriangle
def intri_vert(a,b,c):
#Finding D, E, F
  tan_mat = np.array([[1,1, 0],[0,1,1],[1,0,1]]) 
  cvec = np.array([a,b,c])
  x = np.linalg.inv(tan_mat)@cvec
  A,B,C= tri_vert(a,b,c)
  D = tri_section(B,C,x[1]/x[0])
  E = tri_section(C,A,x[2]/x[1])
  F = tri_section(A,B,x[0]/x[2])
  return D,E,F



#Radius and centre of the incircle
#of triangle ABC
def icentre(A,B,C,k1,k2):
  p = np.zeros(2)
  t = norm_vec(B,C)
  n1 = t/np.linalg.norm(t)
  t = norm_vec(C,A)
  n2 = t/np.linalg.norm(t)
  t = norm_vec(A,B)
  n3 = t/np.linalg.norm(t)
  p[0] = n1@B- k1*n2@C
  p[1] = n2@C- k2*n3@A
  N=np.vstack((n1-k1*n2,n2-k2*n3))
  I=np.matmul(np.linalg.inv(N),p)
  r = n1@(I-B)
  #Intersection
  return I,r

dvec = np.array([-1,1]) 
#Orthogonal matrix
omat = np.array([[0,1],[-1,0]]) 

#Code by GVV Sharma
#March 6, 2019
#released under GNU GPL

#This program plots the incircle of Triangle ABC
import numpy as np
import matplotlib.pyplot as plt

#from incentre import icentre
#if using termux
import subprocess
import shlex
#end if


A,B,C=tri_vert(10,11,7)
len = 100

#A,B,C=np.loadtxt('./codes/vert.dat', dtype='double')
p = np.zeros(2)
k1 = 1
k2 = 1
I,r = icentre(A,B,C,k1,k2)
#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

#Generating circle
theta = np.linspace(0,2*np.pi,len)
x_circ = np.zeros((2,len))
x_circ[0,:] = r*np.cos(theta)
x_circ[1,:] = r*np.sin(theta)
x_circ = (x_circ.T + I).T

#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')

#Plotting circle
plt.plot(x_circ[0,:],x_circ[1,:],label='$incircle$')

P = np.array([4.5,4.5]) 

plt.plot(A[0], A[1], 'o')
plt.text(A[0] * (1 + 0.1), A[1] * (1 - 0.1) , 'A')
plt.plot(B[0], B[1], 'o')
plt.text(B[0] * (1 - 0.2), B[1] * (1) , 'B')
plt.plot(C[0], C[1], 'o')
plt.text(C[0] * (1 + 0.03), C[1] * (1 - 0.1) , 'C')
plt.plot(I[0], I[1], 'o')
plt.text(I[0] * (1 + 0.1), I[1] * (1 - 0.1) , 'I')
plt.plot(I[0], I[1], 'o')
plt.text(I[0] * (1 + 0.1), I[1] * (1 - 0.1) , 'I')

P = np.array([4.5,4.5])
plt.plot(P[0], P[1], 'o')
plt.text(4.5,4.8, 'P')

Q = np.array([0.6,3])
plt.plot(Q[0], Q[1], 'o')
plt.text(0,3, 'Q')

R = np.array([3,0])
plt.plot(R[0], R[1], 'o')
plt.text(3,0.2, 'R')



plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')


#else
#plt.show()

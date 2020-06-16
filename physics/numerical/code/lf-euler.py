'''
Solving the mass spring problem using Euler integration
The time values are stored in a numpy array. Two arrays of the size are
created to store the displacement and velocity.
x1 is displacement and x2 is velocity.
'''
    
import numpy as np
import matplotlib.pyplot as plt

E = np.array([0, 0, 0]) # Electric field components Ex,Ey & Ez
B = np.array([0, 0,  5]) # Magnetic field 

q = 5.0	# spring constant
m = 25.0
dt = .01
tmax = 7
ta = np.arange(0, tmax, dt)		# array of time
N =len(ta)

def f1(x1,x2,t):		# dx1/dt is x2
	dx1dt = x2
	return dx1dt

def f2(x1,x2,t):
	dx2dt = q/m * (E + np.cross(x2,B)) 
	return dx2dt

def euler(X1,X2, f1, f2, t):	# X1 and X2 have three components
	X1 += dt * f1(X1, X2,  t)
	X2 += dt * f2(X1, X2,  t)
	#print (X1, X2)
	return X1,X2

X1a = np.zeros([N,3])   # array of N positions, 3 coordinates
X1a[0] = [0,0,0]        # inital position. x, y and z
X2a = np.zeros([N,3])
X2a[0] = [0,1,0]        # initial velocity. vx, vy and vz

for i in range(N-1):
	X1a[i+1], X2a[i+1] = euler(X1a[i], X2a[i], f1, f2, dt)	
  
plt.plot(X1a[:,0], X1a[:,1])
plt.show()  


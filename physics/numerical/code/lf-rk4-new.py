'''
Solving the mass spring problem using 4th order Runge-Kutta method
The time values are stored in a numpy array. Two arrays of the size are
created to store the displacement and velocity.
x1 is displacement and x2 is velocity.
'''
    
import numpy as np
import matplotlib.pyplot as plt

E = np.array([0, 0, .1]) # Electric field components Ex,Ey & Ez
B = np.array([0, 0,  5]) # Magnetic field 

q = 5.0	# spring constant
m = 25.0
dt = .1
tmax = 50
ta = np.arange(0, tmax, dt)		# array of time
N =len(ta)

def f1(x1,x2,t):		# dx1/dt is x2
	dx1dt = x2
	return dx1dt

def f2(x1,x2,t):
	dx2dt = q/m * (E + np.cross(x2,B)) 
	return dx2dt

def rk4(x1,x2, f1, f2, t):
	k11 = dt*f1(x1,x2,t);
	k21 = dt*f2(x1,x2,t);
	k12 = dt*f1(x1 + k11/2, x2 + k21/2, t + dt/2);
	k22 = dt*f2(x1 + k11/2, x2 + k21/2, t + dt/2);
	k13 = dt*f1(x1 + k12/2, x2 + k22, t + dt/2);
	k23 = dt*f2(x1 + k12/2, x2 + k22/2, t + dt/2);
	k14 = dt*f1(x1 + k13, x2 + k23, t + dt);
	k24 = dt*f2(x1 + k13, x2 + k23, t+dt);
	x1 += (k11+2*k12+2*k13+k14)/6;
	x2 += (k21+2*k22+2*k23+k24)/6;
	return x1,x2

X1a = np.zeros([N,3])   # array of N positions, 3 coordinates
X1a[0] = [0,0,0]     # inital position. x, y and z
X2a = np.zeros([N,3])
X2a[0] = [0,1,0]    # initial velocity. vx, vy and vz

for i in range(N-1):
	X1a[i+1], X2a[i+1] = rk4(X1a[i], X2a[i], f1, f2, dt)	
 
print(X1a[:,0][-1], ta[-1])
 
from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(plt.figure())
ax.plot(X1a[:,0], X1a[:,1], X1a[:,2])		# 3d plot of x, y and z
ax.set_zlabel('Z axis')
plt.show()
 
  


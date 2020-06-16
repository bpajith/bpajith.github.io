'''
Solving the mass spring problem using 4th order Runge-Kutta method
The time values are stored in a numpy array. Two arrays of the size are
created to store the displacement and velocity.
x1 is displacement and x2 is velocity.
'''
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

E = np.array([0, 0, .1]) # Electric field components Ex,Ey & Ez
B = np.array([0, 0,  5]) # Magnetic field 
q = 5.0	# spring constant
m = 25.0
dt = .1
tmax = 50

def solver(X, t0): # X is a six element array, t0 required for the solver
	v = np.array([X[3], X[4], X[5]])  # make the velocity vector
	a = q * (E + np.cross(v,B)) / m   # F = q v x B ; a = F/m
	return [X[3], X[4], X[5], a[0], a[1], a[2]]	

ta = np.arange(0, tmax, dt)		# array of time
x0 = [0, 0, 0, 0, 1, 0]			# initual position & velocity, at t = 0
res = integrate.odeint(solver, x0, ta)	# integrate for position and velocity

print(res[:,0][-1], ta[-1])

from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(plt.figure())
ax.plot(res[:,0], res[:,1], res[:,2])		# 3d plot of x, y and z
ax.set_zlabel('Z axis')
plt.show()
 
  


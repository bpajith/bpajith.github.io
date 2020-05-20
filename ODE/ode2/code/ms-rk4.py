'''
Solving the mass spring problem using Euler integration
The time values are stored in a numpy array. Two arrays of the size are
created to store the displacement and velocity.
x1 is displacement and x2 is velocity.
'''
    
import numpy as np
import matplotlib.pyplot as plt

k = 10.0	# spring constant
m = 1.0
dt = .1
tmax = 5
ta = np.arange(0, tmax, dt)		# array of time
N =len(ta)


def f1(x1,x2,t):		# dx1/dt is x2
	dx1dt = x2
	return dx1dt

def f2(x1,x2,t):
	dx2dt = -k * x1 / m    
	return dx2dt

def rk4(x1,x2, f1, f2, t):
	k11 = dt*f1(x1,x2,t);
	k21 = dt*f2(x1,x2,t);
	k12 = dt*f1(x1+0.5*k11,x2+0.5*k21,t+0.5*dt);
	k22 = dt*f2(x1+0.5*k11,x2+0.5*k21,t+0.5*dt);
	k13 = dt*f1(x1+0.5*k12,x2+0.5*k22,t+0.5*dt);
	k23 = dt*f2(x1+0.5*k12,x2+0.5*k22,t+0.5*dt);
	k14 = dt*f1(x1+k13,x2+k23,t+dt);
	k24 = dt*f2(x1+k13,x2+k23,t+dt);
	x1 += (k11+2*k12+2*k13+k14)/6;
	x2 += (k21+2*k22+2*k23+k24)/6;
	return x1,x2


x1a = np.zeros(N)		# array to store the computed displacements
x2a = np.zeros(N)		# and velocities
x1 = 2
x2 = 0
x1a[0] = x1				# initial value of displacement
x2a[0] = x2				# and position are filled

	
for i in range(N-1):
	x1a[i+1], x2a[i+1] = rk4(x1a[i], x2a[i], f1, f2, dt)	
  
plt.plot(ta, x1a)
plt.plot(ta, x2a)
plt.show()  


def sineFunc(x, a1, a2, a3, a4):
    return a4 + a1* np.sin(abs(a2*(2*np.pi))*x + a3)

import scipy.optimize as optimize

fr = np.sqrt(k/m)/2/np.pi
par = [x1a[0], fr, 0.0, 0.0] 	# Amp, freq, phase , offset
par, pcov = optimize.curve_fit(sineFunc, ta, x1a, par)
print('Frequencies ', fr, par[1])

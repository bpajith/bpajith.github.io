'''
Initial Value Problem with Ordinary Differential Equations. 
We choose a known functions so that the result can be cross checked. 

dF/dx = cos(x) and F(0) = 0  are given.
Using this evaluate F(x) from 0 to pi 
'''

import numpy as np
import matplotlib.pyplot as plt

def f1(x,y):		# derivative of the function to be evaluated
	return np.cos(x)

N = 4				# number of points
h = np. pi/(N-1)	# stepsize. Note that the Number of steps = Number of points -1 

xa = np.zeros(N)   	# numpy array of store the variable values
ya = np.zeros(N)	# numpy array to store results

xa[0] = 0.0
ya[0] = 0.0

def euler(x, y, fxy, h):
	return y + h * fxy(x,y)   # Euler method
	
for i in range(N-1):
	ya[i+1] = euler(xa[i], ya[i], f1, h)
	xa[i+1] = xa[i] + h

print (xa[-1], ya[-1],  'Err = ', ya[-1] - np.sin(xa[-1]))   # the last point and error

plt.plot(xa, ya, 'x-')
plt.plot(xa, np.sin(xa), 'o-')
plt.show()

   

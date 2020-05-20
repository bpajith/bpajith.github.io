'''
Initial Value Problem with Ordinary Differential Equations. 

dF/dx = constant, and F(0) = 0  are given.
'''

import numpy as np
import matplotlib.pyplot as plt


def f1(x,y):		# derivative of the function to be evaluated
	return 4 

N = 4				# number of points
h = np.pi/(N-1)   		# stepsize. Note that the Number of steps = Number of points -1 

xa = np.zeros(N)   	# numpy array of store the variable values
ya = np.zeros(N)		# numpy array to store results

xa[0] = 0.0
ya[0] = 0.0

def euler(x, y, fxy, h):
	return y + h * fxy(x,y)   # Euler method
	
for i in range(N-1):
	ya[i+1] = euler(xa[i], ya[i], f1, h)
	xa[i+1] = xa[i] + h

print (xa[-1], ya[-1],  'Error = ', ya[-1] - 4*xa[-1])
plt.plot(xa, ya, 'x')
plt.plot(xa, 4*xa, '+')
plt.show()
   

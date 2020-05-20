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
h = np.pi/(N-1)   	# stepsize. Note that the Number of steps = Number of points -1 

xa = np.zeros(N)   	# numpy array of store the variable values
ya = np.zeros(N)	# numpy array to store results

xa[0] = 0.0
ya[0] = 0.0

def rk4(x, y, fxy, h):   # x, y , f(x,y)
	k1 = h * fxy(x, y)
	k2 = h * fxy(x + h/2.0, y+k1/2)
	k3 = h * fxy(x + h/2.0, y+k2/2)
	k4 = h * fxy(x + h, y+k3)
	return  y + ( k1/6 + k2/3 + k3/3 + k4/6 )

for i in range(N-1):
	ya[i+1] = rk4(xa[i], ya[i], f1, h)
	xa[i+1] = xa[i] + h

print (xa[-1], ya[-1],  'Err = ', ya[-1] - np.sin(xa[-1]))

plt.plot(xa, ya, 'x-')
plt.plot(xa, np.sin(xa), 'o-')
plt.show()
   

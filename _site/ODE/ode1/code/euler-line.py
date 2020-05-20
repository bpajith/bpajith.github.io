'''
Initial Value Problem with Ordinary Differential Equations. 
We choose a known functions so that the result can be cross checked. 

dy/dx = cos(x) and y(0) = 0  are given.
'''

import numpy as np
import matplotlib.pyplot as plt

def f1(x,y):		# derivative of the function to be evaluated
	return 4

dx = .5	# step size
xmin = 0	# initial value, where 
xmax = np.pi	# calculate up to this only

xa = np.arange(xmin, xmax, dx)	# array of the independent variable
N =len(xa)

ya = np.empty(N)	# numpy array to store results

xa[0] = 0.0	# initial values
ya[0] = 0.0

def euler(x, y, fxy, h):
	return y + h * fxy(x,y)   # Euler method
	
for i in range(N-1):
	ya[i+1] = euler(xa[i], ya[i], f1, dx)

print (xa[-1], ya[-1],  'Err = ', ya[-1] - np.sin(xa[-1]))   # the last point and error

plt.plot(xa, ya, 'x-')
plt.plot(xa, 4*xa, 'o-')
plt.show()

   

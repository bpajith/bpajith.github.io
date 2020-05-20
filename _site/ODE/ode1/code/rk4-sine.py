'''
Initial Value Problem with Ordinary Differential Equations. 
We choose a known functions so that the result can be cross checked. 

dy/dx = cos(x) and y(0) = 0  are given.
'''

import numpy as np
import matplotlib.pyplot as plt

def f1(x,y):		# derivative of the function to be evaluated
	return np.cos(x)

dx = .5	# step size
xmin = 0	# initial value, where 
xmax = np.pi	# calculate up to this only

xa = np.arange(xmin, xmax, dx)	# array of the independent variable
N =len(xa)

ya = np.empty(N)	# numpy array to store results
ya[0] = 0.0	# initial value 

def rk4(x, y, fxy, h):   # x, y , f(x,y)
	k1 = h * fxy(x, y)
	k2 = h * fxy(x + h/2.0, y+k1/2)
	k3 = h * fxy(x + h/2.0, y+k2/2)
	k4 = h * fxy(x + h, y+k3)
	return  y + ( k1/6 + k2/3 + k3/3 + k4/6 )

	
for i in range(N-1):
	ya[i+1] = rk4(xa[i], ya[i], f1, dx)

print (xa[-1], ya[-1],  'Err = ', ya[-1] - np.sin(xa[-1]))   # the last point and error

plt.plot(xa, ya, 'x-')
plt.plot(xa, np.sin(xa), 'o-')
plt.show()

   

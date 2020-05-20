'''
Initial Value Problem with Ordinary Differential Equations. 
We choose a known functions so that the result can be cross checked. 

dy/dx = -Ly  
x : time
y : number ot radio-active atoms 
'''

import numpy as np
import matplotlib.pyplot as plt

L1 = 0.5    # decay constant

def f1(x,y):		# derivative of the function to be evaluated
	return -L1 * y

dx = .1	# step size
xmin = 0	# initial value, where 
xmax = 10	# calculate for 10 seconds

xa = np.arange(xmin, xmax, dx)	# array of the independent variable
N =len(xa)

ya = np.empty(N)	# numpy array to store results
ya[0] = 1000.0		# initial value 

def rk4(x, y, fxy, h):   # x, y , f(x,y)
	k1 = h * fxy(x, y)
	k2 = h * fxy(x + h/2.0, y+k1/2)
	k3 = h * fxy(x + h/2.0, y+k2/2)
	k4 = h * fxy(x + h, y+k3)
	return  y + ( k1/6 + k2/3 + k3/3 + k4/6 )

	
for i in range(N-1):
	ya[i+1] = rk4(xa[i], ya[i], f1, dx)

y1 = ya[0] * np.exp(-L1 * xa)  # from analytical solution   

print ('Final value of y :', ya[-1] , y1[-1])   # numerical and analytical results

plt.plot(xa, ya)
#plt.plot(xa, ya[0] * np.exp(-L1 * xa))
plt.show()

   

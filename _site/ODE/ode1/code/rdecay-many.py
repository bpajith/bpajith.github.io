'''
dy/dx = -Ly  
x : time
y : number ot radio-active atoms 
'''

import numpy as np
import matplotlib.pyplot as plt

NS = 3		   			# number of different siotopes
DC = np.array([0.3, 0.1, 0.1])  	# decay constant for each isotope

def f1(x,y, lam, dx):		# derivative of the function to be evaluated
	return -lam * y * dx

dx = 2	# step size
xmin = 0	# initial value, where 
xmax = 10	# calculate for 10 seconds

xa = np.arange(xmin, xmax, dx)	# array of the independent variable
N =len(xa)

ya = np.zeros([N,NS])	# numpy array to store results
ya[0][0] = 1000.0	# initial value for the first isotope
dya = np.zeros([N,NS])	# numpy array to store decay counts

def rk4(x, y, fxy, lam, dx):   # x, y , f(x,y)
	k1 = dx * fxy(x, y, lam,dx)
	k2 = dx * fxy(x + dx/2.0, y+k1/2, lam,dx)
	k3 = dx * fxy(x + dx/2.0, y+k2/2, lam,dx)
	k4 = dx * fxy(x + dx, y+k3, lam,dx)
	return  y + ( k1/6 + k2/3 + k3/3 + k4/6 )

	
for i in range(N-1):
	for k in range(NS):
		ya[i+1][k] = rk4(xa[i], ya[i][k], f1, DC[k], dx)
		dya[i][k] = ya[i][k] - ya[i+1][k]
		if k < NS-1:		# the decayed count adds up to the daughter count
			ya[i][k+1] += dya[i][k]
	print(ya[i])

#y1 = ya[0] * np.exp(-L1 * xa)  # from analytical solution   

#print ('Final value of y :', ya[-1] , y1[-1])   # numerical and analytical results

plt.plot(xa, ya[:,0])
plt.plot(xa, ya[:,1])
plt.plot(xa, ya[:,2])
#plt.plot(xa, ya[0] * np.exp(-L1 * xa))
plt.show()

   

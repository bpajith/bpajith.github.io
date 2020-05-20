#Solving the mass spring problem using Euler integration

k = 10.0	# spring constant
m = 1.0
dt = 0.1
t = 0
x = 2.0		# initial displacement
v = 0.0		# initial velocity

ta = [t]   # list for time, add start time
xa = [x]   # list for displacement, x = 2 at t= 0
va = [v]   # list for velocity, v=0 at t=0

while t < 5:
	f = -k * x
	v =  v + (f/m) * dt	  # a = F/m;  a = dv/dt
	x =  x + v * dt         # v = dx/dt
	t = t + dt
	ta.append(t)
	xa.append(x)
	va.append(v)

import matplotlib.pyplot as plt
plt.plot(ta, xa)
plt.plot(ta, va)
plt.show()  

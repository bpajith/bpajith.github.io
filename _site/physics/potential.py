from pylab import *          # to install, sudo apt install python3-matplotlib

q1 = 1.60217662e-19
q2 = q1
e0 = 8.85418782e-12         	# permittivity of free space
hcross = 6.62607004e-34 / (2*pi)
m1 = 1.6726219e-27      		# proton mass
m2 = 9.10938356e-31     		# electron mass
mu = m1*m2/(m1+m2)				# reduced mass


rmin = 0.2e-10
rmax = 5.0e-10
ra = linspace(rmin, rmax, 200)  # range of radius for calculation

l = 1

Vl = l * (l + 1) * hcross ** 2 / (2 * mu * ra **2)

Vr = -q1 * q2 / (4*pi*e0*ra)


xlim(0, rmax)
ylim(-2e-18,2e-18)      # adjusted to view the most relevant part of the graph

xlabel('radius')
ylabel('potential')

plot([0, rmax], [0,0])    # x-axis

plot(ra, Vr)
plot(ra, Vl)
plot(ra, Vr + Vl)

show()


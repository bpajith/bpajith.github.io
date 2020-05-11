from pylab import *          # to install, sudo apt install python3-matplotlib

q1 = 1.60217662e-19
q2 = q1
e0 = 8.85418782e-12
hcross = 6.62607004e-34/(2*pi)
m1 = 1.6726219e-27      # proton
m2 = 9.10938356e-31     # electron
mu = m1*m2/(m1+m2)


ra = linspace(0.2, 5, 200) * 1e-10   # r from 10 to 100 Fermi

l = 1

Vl = l * (l + 1) * hcross ** 2 / (2 * mu * ra **2)

Vr = -q1 * q2 / (4*pi*e0*ra)


xlim(0, 5e-10)
ylim(-2e-18,2e-18)
xlabel('radius')
ylabel('potential')

plot([0,5e-10], [0,0])    # x-axis

plot(ra, Vr)
plot(ra, Vl)
plot(ra, Vr + Vl)

show()


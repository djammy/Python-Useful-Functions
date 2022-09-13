import numpy as np
import math
def Lorentzian(parameters,m_eff,epsilon,permitivitty,q,Frequency):
	g, Neq, Npd, gamma, A = parameters
	Neq = Neq*1e17
	Npd = Npd*1e22
	gamma = gamma*1e13
	A = 1
	omega = 2*math.pi*Frequency
	plasma_freq1 = math.sqrt((g*Neq*q**2)/(m_eff*epsilon*permitivitty))
	plasma_freq2 = math.sqrt((g*(Neq+Npd)*q**2)/(m_eff*epsilon*permitivitty))
	fit = A*((1j*(Neq + Npd)*(q**2)*omega)/(m_eff*(omega**2 - plasma_freq2**2 + 1j*omega*gamma))) - \
		A*((1j*(Neq)*(q**2)*omega)/(m_eff*(omega**2 - plasma_freq1**2 + 1j*omega*gamma)))
	return(fit)

def Modified_Surface_Plasmon(parameters,m_eff,m_hol,epsilon,permitivitty,q,Frequency,approach):
	g, Ne_eq, Nh_eq, Nphoto, gamma_e, gamma_h = parameters
	Ne_eq = Ne_eq*1e22
	Nh_eq = Nh_eq*1e22
	Nphoto = Nphoto*1e21
	gamma_e = gamma_e*1e13
	gamma_h = gamma_h*1e12

	# Angular Frequency
	omega = 2*math.pi*Frequency

	# Function Handles
	w_resonance = lambda g,N,q,m,epsilon,permitivitty: math.sqrt((g*N*q**2)/(m*epsilon*permitivitty))
	lorentzian = lambda N,q,omega,m,w_0,gamma: (1j*N*(q**2)*omega)/(m*(omega**2 - w_0**2 + 1j*omega*gamma))
	
	if approach == 'Surface Plasmon':
		fit = lorentzian(Ne_eq+Nphoto,q,omega,m_eff,w_resonance(g,Ne_eq+Nphoto,q,m_eff,epsilon,permitivitty),gamma_e) - \
			lorentzian(Ne_eq,q,omega,m_eff,(w_resonance(g,Ne_eq,q,m_eff,epsilon,permitivitty)),gamma_e)
	elif approach == 'Modified Surface Plasmon':
		fit = lorentzian(Ne_eq+Nphoto,q,omega,m_eff,w_resonance(g,Ne_eq+Nphoto,q,m_eff,epsilon,permitivitty),gamma_e) + \
			lorentzian(Nh_eq+Nphoto,q,omega,m_hol,w_resonance(g,Nh_eq+Nphoto,q,m_hol,epsilon,permitivitty),gamma_h) - \
			lorentzian(Ne_eq,q,omega,m_eff,(w_resonance(g,Ne_eq,q,m_eff,epsilon,permitivitty)),gamma_e) - \
			lorentzian(Nh_eq,q,omega,m_hol,(w_resonance(g,Nh_eq,q,m_hol,epsilon,permitivitty)),gamma_h)
	return(fit)

approach = 'Modified Surface Plasmon'
Frequency = np.linspace(0,3,500)
Frequency = Frequency*1e12
m_eff = 0.02*9.10938215E-31
m_hol = 0.12*9.10938215E-31
epsilon = 12
permittivity = 8.854187817620E-12
q = 1.60217646E-19
vary = [1, 10, 30, 40, 50, 60, 70, 100]
n = 8
y_real = np.zeros((500,n))
y_imag = np.zeros((500,n))

import matplotlib.pyplot as plt
import matplotlib.patches as patch
from matplotlib.pyplot import cm
import pandas as pd

fig = plt.figure()
fig.patch.set_facecolor('white')

for index in range(n):
	parameters = [0.0113, 2.348, vary[index], 10.136, 2.17256, 7.31699]
	y = Modified_Surface_Plasmon(parameters,m_eff,m_hol,epsilon,permittivity,q,Frequency,approach)
	y_real[:,index] = y.real/100
	y_imag[:,index] = y.imag/100
	plt.plot(Frequency*1e-12, y_real[:,index] + 4*index, Frequency*1e-12, y_imag[:,index] - 2*index)
	plt.axhline(y=0, color='k')
# # plt.grid(True, which='both')
	plt.xlim(0,2)
# np.savetxt("Frequency.csv", Frequency*1e-12, delimiter=",")
# np.savetxt("YReal.csv", y_real, delimiter=",")
# np.savetxt("YImag.csv", y_imag, delimiter=",")

# np.savetxt("ConductivityDataReal.csv",conductivityfull[:,:].real/100, delimiter=",")
# np.savetxt("ConductivityDataImag.csv",conductivityfull[:,:].imag/100, delimiter=",")
# np.savetxt("FrequencyFit.csv", Frequency_fit*1e-12, delimiter=",")
# np.savetxt("FitReal.csv",y_fit[:,:].real/100, delimiter=",")
# np.savetxt("FitImag.csv",y_fit[:,:].imag/100, delimiter=",")
# breakpoint()

# df = pd.DataFrame(y_real, index=Frequency*1e-12)
# df.plot.area()



# import matplotlib.pyplot as plt
# import matplotlib.patches as patch
# from matplotlib.pyplot import cm

# fig = plt.figure()
# fig.patch.set_facecolor('white')
# plt.plot(Frequency*1e-12, y.real/100, Frequency*1e-12, y.imag/100)
# plt.axhline(y=0, color='k')
# # plt.grid(True, which='both')
# plt.xlim(0,2)
# # plt.ylim(-25 ,60)
plt.show()
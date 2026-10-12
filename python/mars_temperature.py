import numpy as np

# http://www.bartonlevenson.com/NewPlanetTemps.html
stefan = 5.6704e-8
emiss = 1.0
S = 589.1 # solar constant for mars W/m**2
A = 0.250 # bond albedo (no units)
F = 0.25*S*(1.-A) # climate flux W/m**2
Te = (F/(emiss*stefan))**0.25
Psurf=636. # surface pressure, Pa
Psurf=2e5
f_CO2 = 0.9532
f_H2O = 0. #0.000210
PCO2 = f_CO2 * Psurf 
PH2O = f_H2O * Psurf 
kCO2 = 0.029
kH2O = 0.087
tauCO2 = kCO2*PCO2**0.5
tauH2O = kH2O*PH2O**0.5

tau = tauCO2 + tauH2O
T0 = Te*(1.0+0.75*tau)**0.25

tauVIS = 0.36*(tau - 0.723)**0.411
Labs = F*(1.-np.exp(-tauVIS))
Fsi = F - Labs
F0 = stefan*0.95*T0**4
Fabs = (1.-A)*Fsi + 0.95*(F0-F)
Fc=0.369*Fabs*tau/(-0.6+2*tau)
Fs=F0-Labs-Fc
Ts=(Fs/(0.95*stefan))**0.25

def co2_vap_press(t):
# 	Rv = 8.314/44e-3
# 	Ls = 32.3e3/44e-3
# 	t0=216.8  # k
# 	p0=5.13e5 # pa
# 	co2_vp = p0*np.exp(Ls/Rv*(1./t0-1./t))
	
	co2_vp = np.exp(np.log(760/101.325)-24.03761*np.log(t) - \
		7062.404/(t)+166.3861+3.368548e-5*(t)**2)*133.322
	return co2_vp
	
if __name__ == "__main__":
	print(Ts)
	import matplotlib.pyplot as plt
	plt.ion()
	
	t = np.linspace(120,300,1000)
	plt.plot(t-273.15,co2_vap_press(t)/1e5)
	x1=-125
	x2=-70
	logy1=np.log10(0.006)
	logy2=np.log10(5)
	m=(logy2-logy1)/(x2-x1)
	c=logy2-m*x2
	
	plt.plot(t-273.15,10**(m*(t-273.15)+c),'--b')
	plt.xlabel(r'T ($^\circ$C)')
	plt.ylabel('$P_{CO_2}$ (atm)')
	plt.yscale('log')
	plt.ylim((1e-5,1e1))
	plt.xlim((-160,-50))
	plt.grid()
	plt.show()
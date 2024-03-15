# This is one of the component of the MO cycling loss calculation codes
# This files stores the parameters needed for calculation, including unchange parameters and variables.

import math
import numpy as np
import scipy
import matplotlib.cm as cm
from pylab import *
from scipy.integrate import ode
import time
import os.path
#import astropy.unit as u
from MOparams import *


hostfile = 'BHAC15_0.09Msun.txt' # corresponds to M8 host star
# WATER FUGACITY
# Need to calculate non-dimensional water fugacity at each step.
# Formula original from Li et al. (2008), based on experimental data.

def lnf_w(x):
	'''
	This is a calculation requires one variable x.
	In the following eta() calculation, x = 1.29e16/(f_M*M_E)
	Do not depend on any parameters.
	'''
	# NOTE: Requires non-dimensionalized mantle water mass fraction, and converts it (see Moore & Cowan 2020).
	c0 = -7.9859
	c1 = 4.3559
	c2 = -0.5742
	c3 = 0.0337
	B = 2.0e6
	mu_oliv = 153.31
	mu_wat = 18.01528

 #   M_E, R_E, d_hE, rho_c, rho_m, f_degasE, G, P_E, year, omega_E, x_E, E_a, R_g, T_liq_dry, T_sol_dry = params1
 #   x_h, chi, d_melt, f_M, f_b, f_btwid, alpha, Ra_c, kappa, T_ref, k, beta, c_p, Q_0, T_s, theta, \
 #       T_serp, K_cnst, gamma, r_fug, d_b, f_w_min = params2
 #   omega_0, num_oceans, sigma, mu, t_loss, loss_factor, chi_d, chi_r, dt_nom = params3

	ln_term = np.log((B*x*(mu_oliv/mu_wat))/(1.-(x*(mu_oliv/mu_wat))))

	return c0 + c1*ln_term + c2*ln_term**2. + c3*ln_term**3.


def f_w(x):
	'''
	This is a calculation based on Inf_w(), requires the same input x.
	'''
	return np.exp(lnf_w(x))


# MANTLE VISCOSITY

# Need a function to calculate the viscosity.
# x and T are non-dimensionalized here as well.
def eta(temp, x):
	'''
	This function requires input x for calculations in f_w().
	and requires temp to calculate minimum value.
	parameters: r_fug, M_E, f_M, E_a, R_g
	requires mass to calculate the desiccation limit.
	'''
	eta_scale = (np.exp(lnf_w(5.8e-4)))**(-r_fug)
	eta_0 = 1.0e21/eta_scale

	# Need a minimum fugacity, to avoid f_w --> 0 and eta --> infinity.
	# Let's take the 'desiccation limit' from Moore & Cowan (2020) as the minimum value to be 'trapped' in the mantle.
	if x <= 1.29e16/(f_M*M_E):
		return eta_0*(f_w(1.29e16/(f_M*M_E))**(-r_fug)) * np.exp((E_a/R_g)*((1./temp) - (1./T_ref)))
	else:
		return eta_0*(f_w(x)**(-r_fug)) * np.exp((E_a/R_g)*((1./temp) - (1./T_ref)))


# Planetary parameter scaling laws (for Rp, Rc, h, g, A, V) from Valencia et al. (2006).



# Planetary radius [m]

def Rp(M):
	'''
	This is a calculation of Planetary radius
	Requires radius and mass of the planet for import
	'''
	M_E = 5.972e24
	R_E = 6.371e6
	return R_E*(M/M_E)**0.27

# Core radius [m]

def Rc(M):
	'''
	This is a calculation of core radius, requires mass of radius of the planet.
	'''
	M_E = 5.972e24
	R_E = 6.371e6
	return 0.547*R_E*(M/M_E)**0.25



# Mantle thickness [m]

def h(M):
	'''
	Mantle thickness, implicitly depends on mass.
	'''
	return Rp(M) - Rc(M)


# Surface gravity [m/s^2]

def g(M):
	''' 
	Surface gravity, depends on mass.
	'''
	G = 6.67e-11
	return G*M/(Rp(M)**2.)


# Surface area [m^2]

def A(M):
	'''
	Surface area, depends on planetary radius.
	'''
	return 4.*np.pi*(Rp(M)**2.)


# Mantle volume [m^3]

def V(M):
	'''
	Mantle volume, depends on mass.
	'''
	return (4.*np.pi/3.)*((Rp(M)**3.) - (Rc(M)**3.))




# Below relations based on Earth.

# Mid-ocean ridge length [m]

def L_MOR(M):
	'''
	mid-ocean ridge length, depends on mass.
	Needs to make sure if the scale changes for a heavier planet.
	'''
	return 1.5*2.*np.pi*Rp(M)


# Spreading rate [m/s]

def S(t, temp, x, T_surf, M):
	'''
	Spreading rate, depends on temperature, x, surface temperature and mass.
	t is not used in the computation, bus used as a place holder in case of integration.
	'''
	return 10.76*(kappa**(1.-2.*beta))*(h(M)**(6.*beta-1.))*((alpha*rho_m*g(M) *
															  (temp-T_surf)/(eta(temp, x)*Ra_c))**(2.*beta))



# Rayleigh number

def Ra(temp, x, T_surf, M):
	'''
	Rayleigh number, depends on temperature, x, surface temperature, and mass.
	'''

	return (alpha*rho_m*g(M)*(temp - T_surf)*(h(M)**3))/(eta(temp,x)*kappa)



# Nusselt number, for looking at heat flux -- UNUSED

def Nu(temp, x, T_surf, M):
	'''
	Nusselt numeber, depends on radius.
	Not used
	'''

	return (Ra(temp, x, T_surf, M)/Ra_c)**beta



# Start to run the simulation here



# Function to read in all relevant data from corresponding host star file.


def hoststar(file):

	data = np.loadtxt(file, skiprows=3)
	log_age = data[:, 1]  # log(yr)
	T_eff_star = data[:, 2]  # [K]
	Lbol_Ls = data[:, 3]  # log luminosity
	R_Rs = data[:, 5]

	Ls = 3.839e33  # [erg/s]; solar bolometric luminosity
	Lbol = np.zeros(len(Lbol_Ls))
	Lbol_Ls_star = np.zeros(len(Lbol_Ls))
	for idx in range(0, len(Lbol_Ls)):
		Lbol[idx] = (10.**Lbol_Ls[idx])*Ls/1.0e7  # [W]
		Lbol_Ls_star[idx] = 10.**Lbol_Ls[idx]

	Rs = 6.96e10  # [cm]; solar radius
	Rstar = (R_Rs*Rs)/100.0  # [m]

	star_age = np.zeros(len(log_age))  # [s], for later calculations

	# Convert stellar age to time, in seconds, for compatibility with simulations.
	for idx in range(0, len(log_age)):
		star_age[idx] = ((10.0**log_age[idx])*year)

	return star_age, Lbol, Lbol_Ls_star, T_eff_star


star_age, Lbol, Lbol_Ls_star, T_eff_star = hoststar(hostfile)


# Calculate the top-of-atmosphere bolometric flux, using BHAC15 stellar track.

def S_0(t, a_orb):
	'''
	calculate top-of-atmosphere bolometric flux, using BHAC15 stellar track.
	Requires the hostfile as an import.
	Doesn't depend on mass.
	'''

	t_int = t + (5.0e6*year)  # add offset of 5 Myr, to account for formation
	return np.interp(t_int, star_age, Lbol)/(4.*np.pi*(a_orb**2.))


# Calculate XUV luminosity at time, t. (Ribas et al., 2005; Luger & Barnes 2015)

def L_XUV(t):
	'''
	Calculate XUV luminosity at time t.
	Doesn't depend on mass.
	'''

	t_int = t + (5.0e6*year)  # add offset of 5 Myr, to account for formation

	f_sat = 1.0e-3  # saturation fraction
	# [s], saturation timescale (same as Luger & Barnes 2015)
	t_sat = 1.0e9*year
	beta_XUV = -1.23

	if t <= t_sat:
		return f_sat*np.interp(t_int, star_age, Lbol)
	else:  # t > t_sat
		return f_sat*((t_int/t_sat)**beta_XUV)*np.interp(t_int, star_age, Lbol)


# Calculate top-of-atmosphere XUV flux at time, t. [W/m^2]

def F_XUV(t, a_orb):
	'''
	This is a calculation depends on L_XUV(t).
	Doesn't depend on mass, requires time and a_orb as import.
	'''
	return L_XUV(t)/(4.*np.pi*(a_orb**2.))


# Calculate the effective temperature of the planet at time, t. [K]

def f_T_eff(t, a_orb):
	'''
	Calculate the effective temperature of the planet at time t.
	Depends on S_0(), doesn't depends on mass.
	'''
	return ((S_0(t, a_orb)*(1.0-alb))/(4.0*sigma_sb))**(0.25)



# Calculate the skin temperature (i.e., the isothermal stratosphere temperature). [K]

def T_strat(t, a_orb):
	'''
	Calculate the isothermal stratosphere temperature (skin temp).
	Depends on S_0, doesn't depend on mass.
	'''

	return ((1./2.)**0.25)*((S_0(t, a_orb)*(1.0-alb))/(4.0*sigma_sb))**(0.25)





# Parameterize depth of a solidifying MO as a decreasing exponential.

def d_MO(t, tau_MO, M):
	'''
	Parameterize depth of a solidifying MO as a decreasing exponential.
	Depth of MO at time, t, for a given MO solidification timescale (free parameter).
	Based on Fig. 3 of Barth et al. 2020, Fig. 4 of Schaefer et al. 2016
	Implicitly depends on mass because the calculation requires mantle thickness h(M)FT
	'''
	d_MO_0 = h(M)/1000.  # ENTIRE MANTLE BEGINS MOLTEN [km]
	beta_MO = 1./(np.exp(1.)-1)

	return (beta_MO*d_MO_0)*(np.exp((-t/tau_MO) + 1.) - 1)



# Temperature of MO as a function of MO depth, simply using liquidus temperature (Eqn. 1 of Elkins-Tanton 2008). [K]

def T_MO(t, tau_MO, M):
	'''
	Calculate temperature of MO as a function of MO depth.
	Depends on mass implicitly, utilized Rp(M) function to calculate r.
	Note the unit is in kilometer.
	'''
	r = (Rp(M) - d_MO(t, tau_MO, M))/1000.  # [km]

	return (-1.16e-7*(r**3.)) + (0.0014*(r**2.)) - (6.382*r) + 1.444e4




# Parameterization of water in each inventory during MO solidification.
# Details in Moore et al. (2022); simplified from Boukare et al. (2019).
# RECALL r = R_p - d_MO(t)


def C_l(C_0, r, M):
	'''
	Concentration of water in liquid phase.
	Depends on C_0 and radius, so implicitly on mass.
	'''
	R_p = Rp(M)/1000
	R_c = Rc(M)/1000
	return C_0*((R_p**3.-R_c**3.)/(R_p**3.-r**3.))**(1.-D)


# Concentration of water in solid phase
def C_s(C_0, r):
	'''
	Concentration of water in solid phase.
	Depends on C_l, hence on mass.
	'''
	return D*C_l(C_0, r, M)


# Saturation radius [km]

def R_sat(C_0, M):
	'''
	Saturation radius.
	If C_0 > C_sat, use C_0 = C_sat.
	Depends on C_0 and radius of the planet, so implicitly on mass.
	Units in kilometer.
	'''
	R_p = Rp(M)/1000
	R_c = Rc(M)/1000

	return (R_p**3. - ((R_p**3. - R_c**3.)/((C_sat/C_0)**(1./(1.-D)))))**(1./3.)


# Mass of water in unsaturated MO [kg]

def M_MO_unsat(C_0, r, M):
	'''
	Mass of water in unsaturated MO.
	Depends on r and radius of the planet, so implicitly depends on mass.
	'''
	R_p = Rp(M)/1000
	R_c = Rc(M)/1000
	rho_MO = 3.0e3*1.0e9 #density of MO [kg/km^3]; assumed to be same for solid mantle & MO for mass balance purposes

	return C_l(C_0, r, M)*(4.*np.pi/3.)*rho_MO*(R_p**3. - r**3.)


# Mass of water in unsaturated solid mantle [kg]

def M_SM_unsat(C_0, r, M):
	'''
	Mass of water in unsaturated solid mantle.
	Depends on r and radius of the planet, so implicitly depends on mass.
	'''
	R_p = Rp(M)/1000
	R_c = Rc(M)/1000
	rho_MO = 3.0e3*1.0e9 #density of MO [kg/km^3]; assumed to be same for solid mantle & MO for mass balance purposes
	return C_0*(4.*np.pi/3.)*(rho_MO)*((R_p**3. - R_c**3.)**(1.-D)) *\
		(((R_p**3. - R_c**3.)**D) - ((R_p**3. - r**3.)**D))


# Mass of water in saturated MO [kg]

def M_MO_sat(r, M):
	'''
	Mass of water in saturated MO.
	Similar to M_SM_unsat, depends on r and radius of the planet, so the mass.
	'''
	R_p = Rp(M)/1000
	R_c = Rc(M)/1000
	rho_MO = 3.0e3*1.0e9 #density of MO [kg/km^3]; assumed to be same for solid mantle & MO for mass balance purposes

	return C_sat*(4.*np.pi/3.)*rho_MO*(R_p**3. - r**3.)


# Mass of water in saturated solid mantle [kg]

def M_SM_sat(C_0, r, M):
	'''
	Mass of water in saturated solid mantle.
	Depends on r and radius of the planet, so depends on mass.
	'''
	rho_MO = 3.0e3*1.0e9 #density of MO [kg/km^3]; assumed to be same for solid mantle & MO for mass balance purposes

	return M_SM_unsat(C_0, R_sat(C_0, M), M) + D*C_sat*(4.*np.pi/3.)*(rho_MO)*(r**3.-R_sat(C_0, M)**3.)





# Loss rate equations from Luger & Barnes (2015)

# Energy-limited loss rate [kg/s]
def loss_rate_EL_MO(F_XUV, M):
	'''
	Energy-limited loss rate.
	Depends on mass.
	'''
	# R_XUV = R_p for simplicity
	# M = planetary mass

#	G = 6.67e-11  # [m^3/(kg s^2)]; gravitational constant, already available in parameter module.
	eps_XUV = 0.1  # XUV efficiency factor; varies within literature

	return eps_XUV*np.pi*F_XUV*(Rp(M)**3.)/(G*M)


# Diffusion-limited loss rate [kg/s]

def loss_rate_DL_MO(M):
	'''
	Diffusion-limited loss rate [kg/s]
	Depends on mass for radius and local gravity in calculation.
	'''
	g_p = g(M)*100.  # [cm s^-2], to match the units of b
	k_B = 1.38e-23  # [m^2 kg s^-2 K^-1]; Boltzmann constant

	# Mole fractions of H & O (assuming all H2O is dissociated & pure H2O atmosphere)
	X_H = 2./3.
	X_O = 1./3.

	m_H = 1.66e-27  # [kg]; mass of H atom
	m_O = 16.*m_H  # mass of O atom

	# Just take the thermospheric temperature as 400 K at all times, for now.
	b = 4.8e17*(400.**0.75)  # [cm^-1 s^-1]

	return m_H*(np.pi*(Rp(M)**2.)*b*g_p*(m_O-m_H))/(k_B*400.*(1.+(X_O/X_H)))



# Loss rate; EL during RG, minimum of EL/DL otherwise

def f_loss_MO(t, M, a_orb):
	'''
	Calculate loss rate. During runaway greengas stage the function returns energy-limited loss, and EL/DL otherwise.
	Depends on mass and also requires the hostfile as an input.
	'''
	# Note that the "check" for not losing more water than on the surface is in the MO loop itself

	# Runaway greenhouse limit; divide by 4 because planet is a sphere
	if S_0(t, a_orb)/4. > (325./(1.-alb)):

		# Energy-limited loss; set R_XUV = R_p, for simplicity
		d_loss = loss_rate_EL_MO(F_XUV(t, a_orb), M)  # [kg/s]

	else:

		# Take the lower value of EL and DL loss rates.
		d_loss = np.minimum(loss_rate_EL_MO(F_XUV(t, a_orb), M),
							loss_rate_DL_MO(M))  # [kg/s]

	return d_loss




# Function to determine orbital distances

def f_a_orb(t, HZ_data_file, hostfile):
	'''
	This function extract orbital distances from Kopparapu et al. (2013) HZ data. 
	Also requires Baraffe et al. 2015 (BHAC15) startrack data file as an import.
	This function also requires a time in years to find orbital distances at the time. 
	For the original MO calculation, we take t = 4.5Gyr, roughly the age of earth.
	'''
	# Read in Kopparapu et al. (2013) HZ data, for different stellar hosts.
	data = np.loadtxt(HZ_data_file, skiprows=2)
	T_eff_kopp = data[:, 0]  # [K]
	# All below values give S_eff -- see Eqn.(2) of Kopparapu+(2013)
	RV_kopp = data[:, 1]  # Recent Venus
	RG_kopp = data[:, 2]  # Runaway Greenhouse
	MG_kopp = data[:, 3]  # Maximum Greenhouse
	EM_kopp = data[:, 4]  # Early Mars

	# Read in data from the hostfile(BHAC15)
	star_age, Lbol, Lbol_Ls_star, T_eff_star = hoststar(hostfile)
	# Roughly find limits at t = 4.5 Gyr (~age of Earth)
	T_eff_star_t = np.interp(4.5e9*year, star_age, T_eff_star)
	Lbol_Ls_star_t = np.interp(4.5e9*year, star_age, Lbol_Ls_star)

	# Effective flux at given orbital distance
	S_eff_RV = np.interp(T_eff_star_t, T_eff_kopp, RV_kopp)  # recent Venus
	S_eff_RG = np.interp(T_eff_star_t, T_eff_kopp, RG_kopp)  # runaway greenhouse
	S_eff_MG = np.interp(T_eff_star_t, T_eff_kopp, MG_kopp)  # maximum greenhouse
	S_eff_EM = np.interp(T_eff_star_t, T_eff_kopp, EM_kopp)  # early Mars

	# Orbital distance
	a_RV = (Lbol_Ls_star_t/S_eff_RV)**0.5  # [AU]
	a_RG = (Lbol_Ls_star_t/S_eff_RG)**0.5  # [AU]
	a_MG = (Lbol_Ls_star_t/S_eff_MG)**0.5  # [AU]
	a_EM = (Lbol_Ls_star_t/S_eff_EM)**0.5  # [AU]

	return a_RV, a_RG, a_MG, a_EM


# Seafloor pressure
def P(s, M):
	'''
	Calculates seafloor pressure.
	Requires mass as an input.
	'''

	return M*g(M)*s/(f_b*4*np.pi*(Rp(M)**2))


# Hydration depth
def d_h(temp, T_surf, x, M):
	'''
	Calculates hydration depth, requires temperature, surface temperature, x, and mass.
	'''
	return (h(M)**(1.0-3.0*beta))*((temp-T_surf)**(-(1.0+beta))) *\
		(T_serp-T_surf)*((eta(temp, x) * kappa*Ra_c/(alpha*rho_m*g(M)))**beta)



# Radionuclide heating, using 4 most common Earth mantle species (i.e., nominal bulk silicate Earth, 21 ppb U)
# (Schaefer & Sasselov 2015, Eqn. (4), with constants as given in Laura's MATLAB code (see refs there))

def Q_sum(t):
	'''
	Radionuclide heating, using 4 most common Earth mantle species.
	Requires time as input, not mass depenedent.
	(Schaefer & Sasselov 2015, Eqn. (4), with constants as given in Laura's MATLAB code (see refs there))
	'''
	# Mantle concentration of element by mass (Schubert et al. 2001, Ch. 4)
	Uran = 21e-9  # [U] = 21 ppb
	C_238U = 0.9927 * Uran
	C_235U = 0.0072 * Uran
	C_40K = 1.28 * Uran
	C_232Th = 4.01 * Uran

	# Heat production per unit mass (Schubert et al. 2001, Ch. 4)
	H_238U = 9.37e-5
	H_235U = 5.69e-4
	H_232Th = 2.69e-5
	H_40K = 2.79e-5

	# Decay constants
	lam_238U = 0.155e-9
	lam_235U = 0.985e-9
	lam_232Th = 0.0495e-9
	lam_40K = 0.555e-9

	s_in_yr = 365.25*24*60*60
	t /= s_in_yr


	return rho_m*(C_238U*H_238U*np.exp(lam_238U*(4.6e9-t)) + C_235U*H_235U*np.exp(lam_235U*(4.6e9-t)) +
				  C_232Th*H_232Th*np.exp(lam_232Th*(4.6e9-t)) + C_40K*H_40K*np.exp(lam_40K*(4.6e9-t)))


# Saturation water vapour pressure, from Nakajima et al. (1992)
def P_star(T):
	'''
	Calculate saturation water vapour pressure, from Nakajima et al. (1992)
	Depends on temperature
	'''
	return P_0_star*np.exp(-l_c_NJ/(R_gas_NJ*T))



# Calculate surface temperature based on Fig. 2b of Goldblatt et al. (2013).
# If above Simpson-Nakajima limit, set to same RG value as Barth et al. (2020): T_surf = 1800 K.

def T_surf_OLR(t, alb, a_orb):
	'''
	Calculate surface temperature based water condensation limit of Turbet et al. (2021).
    The surface temperature is fixed at 1800 K when T_surf > threshold, and 293.15 otherwise.
	'''

	if S_0(t, a_orb)/4. > (325./(1.-alb)):  # runaway greenhouse

		return 1800.  # ~maximum T_surf corresponding to given OLR

	else:  # not in runaway greenhouse

		return 293.15 #20 C


# Energy-limited escape rate
def loss_rate_EL(F_XUV, R_XUV, M):  # [kg/s]
	'''
	Energy-limited escape rate. Similar as in MO stage.
	Depends on mass, the rate is a constant for all t.
	'''
	eps_XUV = 0.1

	return eps_XUV*np.pi*F_XUV*Rp(M)*(R_XUV**2.)/(G*M)


# Diffusion-limited escape rate
def loss_rate_DL(M, T = 400.):
	'''
	Diffusion-limited escape rate.
	Depends on mass.
	!!!Note!!! This function assumes thermospheric temperature at 400K for all time.
	'''
	g_p = g(M)*100.  # [cm s^-2], to match the units of b

	# Just take the thermospheric temperature as 400 K at all times, for now.
	b = 4.8e17*(T**0.75)  # [cm^-1 s^-1]

	return m_H*(np.pi*(Rp(M)**2.)*b*g_p*(m_O-m_H))/(k_B*400.*(1.+(X_O/X_H)))


# Thermal evolution of mantle temperature [K]
def f_delta_temp(t, temp, W_m, W_s, T_surf, M):
	'''
	Calculate thermal evolution of mantle temperature in K.
	Depends on mass.
	'''
	x = W_m/(f_M*M)  # mantle water mass fraction
	s = W_s/M  # surface water mass fraction

	return (1./(rho_m*c_p))*((Q_sum(t)) -
							 (k*A(M)*(temp-T_surf)/(h(M)*V(M)))*(Ra(temp, x, T_surf, M)/Ra_c)**beta)


# Degassing rate
def f_degas(t, temp, W_m, W_s, T_surf, M):
	'''
	Calculates degassing rate.
	Depends on mass.
	t is not used in the calculation, but it is a place holder to track indices of input array.
	'''
	x = W_m/(f_M*M)
	s = W_s/M

	# Mantle temperature check (if below solidus, no melting)
	if temp > T_sol_dry-K_cnst*x**gamma:

		# Piecewise limit from Cowan & Abbot (2014)
		if x > 0. and s > 0. and f_degasE*(P(s, M)/P_E)**(-mu) < 1.:
			return L_MOR(M)*S(t, temp, x, T_surf, M)*(x*rho_m*d_melt *\
													(f_degasE*(P(s, M)/P_E)**(-mu)) *\
													((temp-(T_sol_dry-K_cnst*x**gamma))/(T_liq_dry-T_sol_dry))**theta)

		elif x > 0. and s > 0. and f_degasE*(P(s, M)/P_E)**(-mu) >= 1.:
			return L_MOR(M)*S(t, temp, x, T_surf, M)*(x*rho_m*d_melt *\
													(1.)*((temp-(T_sol_dry-K_cnst*x**gamma))/(T_liq_dry-T_sol_dry))**theta)

		elif x > 0. and s <= 0.:
			return L_MOR(M)*S(t, temp, x, T_surf, M)*(x*rho_m*d_melt *\
												 (1.)*((temp-(T_sol_dry-K_cnst*x**gamma))/(T_liq_dry-T_sol_dry))**theta)

		else:  # x <= 0.; no degassing
			return 0.

	else:  # temp <= T_sol_dry-K_cnst*x**gamma; no degassing
		return 0.



# Regassing rate (adding regassing check of Cowan & Abbot 2014 to T-dependent regassing of Schaefer & Sasselov 2015)

def f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb):
	'''
	Calculate Regassing rate.
	Depends on mass. This function requires hostfile to be defined in previous codes.
K_cnst	'''
	x = W_m/(f_M*M)
	s = W_s/M
	if S_0(t, a_orb)/4. > (325./(1.-alb)):

		# Runaway greenhouse; all water in atmosphere so NO regassing!!!
		return 0.

	else:

		# Piecewise limit from Cowan & Abbot (2014)
		if x > 0. and s > 0. and d_h(temp, T_surf, x, M)*(P(s, M)/P_E)**sigma < d_b:
			return L_MOR(M)*S(t, temp, x, T_surf, M)*(x_h*rho_c*chi_r *
												 d_h(temp, T_surf, x, M)*(P(s, M)/P_E)**sigma)

		elif x > 0. and s > 0. and d_h(temp, T_surf, x, M)*(P(s, M)/P_E)**sigma >= d_b:
			return L_MOR(M)*S(t, temp, x, T_surf, M)*(x_h*rho_c*chi_r*d_b)

		elif x <= 0. and s > 0.:  # constant spreading rate assumed if mantle dry
			return L_MOR(M)*(0.1/year)*(x_h*rho_c*chi_r*d_b)

		elif s <= 0.:  # no regassing
			return 0.


# Water loss rate; either energy- or diffusion-limited
def f_loss(t, W_m, W_s, alb, dt, M, a_orb):

	'''
	Calculates water loss rate, either energy or diffusion-limited.
	Depends on mass. This function requires hostfile to be defined in previous codes.
	'''
	x = W_m/(f_M*M)
	s = W_s/M

	# Runaway greenhouse limit; divide by 4 because planet is a sphere
	if S_0(t, a_orb)/4. > (325./(1.-alb)):

		# Energy-limited loss; set R_XUV = R_p, for simplicity
		M_loss = loss_rate_EL(F_XUV(t, a_orb), Rp(M), M)  # [kg/s]

	else:

		# Take the lower value of EL and DL loss rates.
		M_loss = np.minimum(loss_rate_EL(F_XUV(t, a_orb), Rp(M), M),
							loss_rate_DL(M))  # [kg/s]

	# Our timesteps are comparatively long, replenishment time above exobase is ~instantaneous; therefore, we can take
	# water directly from the surface reservoir.
	# Piecewise definition, so you don't lose more than on the surface
	if W_s == 0.:
		return 0.
	elif M_loss < W_s/dt:
		return M_loss
	else:
		return W_s/dt

	# return 0. #no loss


# Change in mantle water with time
def f_delta_W_m(t, temp, W_m, W_s, T_surf, alb, dt, M, a_orb, W_BMO_init, tau_BMO, tau_RG):
	'''
	Calculate change in mantle water with time.
	Depends on mass. This function requires hostfile to be defined in previous codes.
	'''
	x = W_m/(f_M*M)
	s = W_s/M

	# Rate of degassing from basal MO to solid mantle, PER SECOND, for lifetime of BMO
	if t <= tau_BMO:
		
		W_BMO_rate = W_BMO_init/(tau_BMO-tau_RG)
		
	else: #t > tau_BMO
		
		W_BMO_rate = 0.
	
	d_W_m = f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
		f_degas(t, temp, W_m, W_s, T_surf, M)+\
			W_BMO_rate #slow but constant degassing of water from BMO to mantle, PER SECOND
	
	# Mantle water capacity check
	if W_m + (d_W_m*dt) <= (M/M_E*12.0)*1.4e21:
		return d_W_m
	else:
		return (((M/M_E*12.0)*1.4e21) - W_m)/dt


# Change in surface water with time

def f_delta_W_s(t, temp, W_m, W_s, T_surf, alb, dt, M, a_orb):
	'''
	Calculate change in surface water with time.
	Depends on mass. This function requires hostfile to be defined in previous codes.
	'''
	x = W_m/(f_M*M)
	s = W_s/M

	d_W_m = f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
		f_degas(t, temp, W_m, W_s, T_surf, M)

	# Mantle water capacity check
	if W_m + (d_W_m*dt) <= (M/M_E*12.0)*1.4e21:
		return f_degas(t, temp, W_m, W_s, T_surf, M) - \
			f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
			f_loss(t, W_m, W_s, alb, dt, M, a_orb)
	else:
		return -((((M/M_E*12.0)*1.4e21) - W_m)/dt) - \
			f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
			f_loss(t, W_m, W_s, alb, dt, M, a_orb)



# Cycling equation, to be integrated during the loop

def f_cycling(t, z, T_surf, alb, dt, M, a_orb, W_BMO_rate, tau_BMO, tau_RG):
	'''
	Defines the cycling equation to be integrated into the loop.
	'''

	temp = z[0] #mantle temperature
	W_m = z[1] #mantle water mass
	W_s = z[2] #surface water mass
	
	return [f_delta_temp(t, temp, W_m, W_s, T_surf, M),\
			f_delta_W_m(t, temp, W_m, W_s, T_surf, alb, dt, M, a_orb, W_BMO_rate, tau_BMO, tau_RG),\
			f_delta_W_s(t, temp, W_m, W_s, T_surf, alb, dt, M, a_orb)]



# The following equations are only used in MO cycling stage.
###############################################################################################
# Cycling equation, in the magma ocean stage --- # UNUSED; DIFFERENT MECHANISMS AND FUNCTIONS USED, SEE MAIN NOTEBOOK

# Change in mantle water with time
def f_delta_W_m_MO(t, temp, W_m, W_s, T_surf, alb, dt, M, a_orb):

	x = W_m/(f_M*M)
	s = W_s/M

	f_loss(t, W_m, W_s, alb, dt, M, a_orb)
	d_W_m = f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
			f_degas(t, temp, W_m, W_s, T_surf, M)
    
    # Mantle water capacity check; maximum of 12 Earth Oceans
	if W_m + (d_W_m*dt) <= 12.0*1.4e21:
		return d_W_m
	else:
		return ((12.0*1.4e21) - W_m)/dt

# Change in surface water with time
def f_delta_W_s_MO(t, temp, W_m, W_s, T_surf, T_strat, alb, dt, M, a_orb):
   
	x = W_m/(f_M*M)
	s = W_s/M
    
	d_W_m = f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
			f_degas(t, temp, W_m, W_s, T_surf, M)
    
	# Mantle water capacity check
	if W_m + (d_W_m*dt) <= 12.0*1.4e21:
		return f_degas(t, temp, W_m, W_s, T_surf, M) - \
			f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
			f_loss(t, W_m, W_s, alb, dt, M, a_orb)
	else:
		return -(((12.0*1.4e21) - W_m)/dt) - \
			f_regas(t, temp, W_m, W_s, T_surf, alb, M, a_orb) - \
			f_loss(t, W_m, W_s, alb, dt, M, a_orb)


# Cycling equation, to be integrated during the loop
def f_cycling_MO(t, z, T_surf, T_strat, alb, dt, M, a_orb):
    
	temp = z[0] #mantle temperature
	W_m = z[1] #mantle water mass
	W_s = z[2] #surface water mass
    
	return [f_delta_temp(t, temp, W_m, W_s, T_surf, M),\
			f_delta_W_m_MO(t, temp, W_m, W_s, T_surf, alb, dt, M, a_orb),\
			f_delta_W_s_MO(t, temp, W_m, W_s, T_surf, T_strat, alb, dt, M, a_orb)]
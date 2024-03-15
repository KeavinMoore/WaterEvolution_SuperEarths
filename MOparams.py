# This is one of the component of the MO cycling loss calculation codes
# This files stores the parameters needed for calculation, including unchange parameters and variables.
import math
import numpy as np
import scipy
import matplotlib.cm as cm
from pylab import *
from scipy.integrate import ode
#from mpl_toolkits.axes_grid.inset_locator import inset_axes
import time
import os.path

#import astropy.unit as u


# 1) Unchanging constants
# [M_E, R_E, d_hE, rho_c, rho_m, f_degasE, G, P_E, year, omega_E, x_E, E_a, R_g, T_liq_dry, T_sol_dry]
# One can change the numbers here so the file "const_params_input.txt" is no longer necessary

M_E = 5.972e24
M = 1.*M_E # CHANGE FOR DIFFERENT MASS PLANETS
D = 0.001 #solid-liquid H2O distribution coefficient, =0.2 to mimic basal magma ocean behaviour
NUM_EARTH_OCEANS = np.array([2., 4., 6., 8.])*(1.4e21/M_E) # water mass fraction of initial water inventory of the planet -- scaled to 2, 4, 6, 8 Earth Oceans for 1 Earth Mass
NUM_EARTH_OCEANS_FILE = (M/M_E)*np.array([2, 4, 6, 8]) # used for saving files for writing and reading
HABITABLE_ZONE = 'ALL' #the region of the habitable zone the planet is in; Inner, Middle, or Outer; or ALL to test all three.
R_E = 6.371e6 
d_hE = 3.0e3
rho_c = 3.0e3
rho_m = 3.3e3
f_degasE = 0.9
G = 6.67e-11
P_E = 4.0e7
year = 3.154e7
omega_E = 6.2e-4
x_E = 5.8e-4
E_a = 335.0e3
R_g = 8.314
T_liq_dry = 1498.15
T_sol_dry = 1248.15

# M_E, R_E, d_hE, rho_c, rho_m, f_degasE, G, P_E, year, omega_E, x_E, E_a, R_g, T_liq_dry, T_sol_dry = params1

# 2) Constants that could be varied
# [x_h, chi, d_melt, f_M, f_b, f_btwid, alpha, Ra_c, kappa, T_ref, k, beta, c_p, Q_0, T_s, theta, \
#  T_serp, K_cnst, gamma, r_fug, d_b]

x_h = 0.05
chi = 0.23
d_melt = 60.0e3
f_M = 0.68
f_b = 0.9
f_btwid = 1.3
alpha = 2.0e-5
Ra_c = 1100.0
kappa = 1.0e-6
T_ref = 1600.0
k = 4.2
beta = 0.3
c_p = 1200.0
Q_0 = 5.0e8
T_s = 280.0
theta = 1.5
T_serp = 973.0
K_cnst = 43.0
gamma = 0.75
r_fug = 1.0
d_b = 6.0e3
f_w_min = 0.09

# x_h, chi, d_melt, f_M, f_b, f_btwid, alpha, Ra_c, kappa, T_ref, k, beta, c_p, Q_0, T_s, theta, \
#    T_serp, K_cnst, gamma, r_fug, d_b, f_w_min = params2

# 3) Parameters that could be different between, e.g., different planets (currently based on Earth)
# [omega_0, num_oceans, sigma, mu, t_loss, loss_factor, chi_d, chi_r, dt]


omega_0 = 2.3e-4
num_oceans = 1.0
sigma = 1.0
mu = 1.0
t_loss = 1.0e8
loss_factor = 7.5
chi_d = 0.02
chi_r = 0.03
dt_nom = 1.0e-5
# omega_0, num_oceans, sigma, mu, t_loss, loss_factor, chi_d, chi_r, dt_nom = params3

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

# NEW constants, for updated atmosphere & loss calculations




R_p = Rp(M)/1000. #[km] #6400 km
R_c = Rc(M)/1000. #[km] #3000 km
rho = 3.0e3*1.0e9 #[kg/km^3] 
rho_MO = rho #density of MO [kg/km^3]; assumed to be same for solid mantle & MO for mass balance purposes

# Saturation concentration of water in MO
C_sat = 0.01 #H2O saturation limit of MO



sigma_sb = 5.67e-8 #[W/m^2/K^4]; Stefan-Boltzmann constant
m_H = 1.66e-27 #[kg]; mass of H atom
m_O= 16.*m_H
k_B = 1.38e-23 #[m^2 kg s^-2 K^-1]; Boltzmann constant
m_H2O = 18.02*m_H #molecular mass of water
m_air = 28.7*m_H #molecular mass of air on Earth
rho_water = 997. #[kg/m^3]; density of water
cp_water_vapour = 1.996e3 #[J/(kg K)]; specific heat capacity of water vapour (steam)
alb = 0.3 #albedo; value roughly similar to Earth
#alb = 0.75 #albedo; upper limit to roughly test "steam atmosphere"

# Other parameters used in functions.
k_B = 1.38e-23  # [m^2 kg s^-2 K^-1]; Boltzmann constant
# Mole fractions of H & O (assuming all H2O is dissociated)
X_H = 2./3.
X_O = 1./3.

def available_var():
	'''
	show all available variables defined in this module
	'''
	print("available variables are:")
	print("M_E, R_E, d_hE, rho_c, rho_m, f_degasE, G, P_E, year, omega_E, x_E, E_a, R_g, T_liq_dry, T_sol_dry")
	print("x_h, chi, d_melt, f_M, f_b, f_btwid, alpha, Ra_c, kappa, T_ref, k, beta, c_p, Q_0, T_s, theta, T_serp, K_cnst, gamma, r_fug, d_b, f_w_min")
	print("omega_0, num_oceans, sigma, mu, t_loss, loss_factor, chi_d, chi_r, dt_nom")
	print("sigma_sb, m_H, k_B, m_H2O, m_air, rho_water, cp_water_vapour, alb")

# Loop to determine whether or not MO is saturated & degas an atmosphere if it is
dt_MO = 2.0e3*year #2000 yr per timestep

# Variables/functions, required for the adiabats and water saturation pressure.
P_0_star = 1.4e11 #[Pa]
l_c_NJ = 43655 #[J/mol] 
R_gas_NJ = 8.314 #[J/mol/K]


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


# Things that need to be calculated with constants and above functions:
# NOTE: If dimensions are being used anyway, many of these are redundant/won't be used in the code.
# HOWEVER, the way the functions are written depend on params4 as an input, so leave for now.
# Details can mostly be found in Komacek & Abbot (2016) and Moore & Cowan (2020).
omg_Etwid = num_oceans*omega_E/omega_0
omegatwid = omg_Etwid/(f_btwid)
T_stwid = T_s/T_ref
T_mtwid = T_ref*R_g/E_a
T_liqtwid = T_liq_dry/T_ref
T_soltwid = T_sol_dry/T_ref
Ktwid = K_cnst/T_ref
Xtwid = x_h*rho_c*chi*d_hE*f_M/(rho_m*d_melt*f_degasE*omega_0*f_btwid)
Pi = (rho_m*d_melt*chi_d*(omega_0*f_btwid/f_M)*((T_liqtwid-T_soltwid)**-theta))
lambdatwid = Ktwid*(omega_0*f_btwid/f_M)**(gamma)
eta_scale = (np.exp(lnf_w(5.8e-4)))**(-r_fug)
eta_0 = 1.0e21/eta_scale
#f_w = np.exp(lnf_w(1.,params1,params2,params3))
f_wtwid_min = 1.0e-5 #CONSERVATIVELY CHOSEN FOR NOW
E = rho_m*d_melt*f_degasE*omega_0*f_btwid/f_M
#tau_heat = Q_0/(rho_m*c_p*T_ref)
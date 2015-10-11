'''MSSM Back Propagation: Plots and Analysis
Jonny Arnold'''

import models as MO
import methods as ME
import plots as P

import numpy as N
import scipy.optimize as SO

# See 'sm-mssm.py' for a simple back-propagation graph (use plot_coupling)

def mc_min_mssm(log_m_susy_range, log_m_gut_range, a_inv_gut_range, **kwargs):
	'''Uses a Monte-Carlo Levenberg-Marquardt (MCLM) algorithm to determine the minimum of chi^2 for the 2-loop MSSM'''
	# NOTE: Uncaught errors can happen if ranges are too extreme!
	
	# Minimising function for LM algorithm
	def f(x):
		# Set up model
		log_m_susy, log_m_gut, a_inv_gut = x[0], x[1], x[2]
		model = MO.MSSMModel(loops=2, log_m_susy=log_m_susy)
		model.set_init(log_m_gut, N.array([a_inv_gut, a_inv_gut, a_inv_gut]), model.a_inv_error_from_exp())
		# Propagate model back to initial values, return chi^2
		log_mass, a_inv = model.propagate(0, return_masses=True)
		return model.get_chi_mz(log_mass, a_inv)
	
	# Minimising function for MCLM algorithm
	def g(x):
		result = SO.leastsq(f,x)
		if result[1] in (1,2,3,4): return N.sum(f(result[0])**2)
		else: return N.inf
	
	# Do the minimisation!
	ME.mc_minimise(g, [log_m_susy_range, log_m_gut_range, a_inv_gut_range], **kwargs)

# Minimise the MSSM		
#mc_min_mssm(log_m_susy_range=[3.3,3.6], log_m_gut_range=[16.0,16.2], a_inv_gut_range=[25.4,25.6])

def f(x):
		# Set up model
		log_m_susy, log_m_gut, a_inv_gut = x[0], x[1], x[2]
		model = MO.MSSMModel(loops=2, log_m_susy=log_m_susy)
		model.set_init(log_m_gut, N.array([a_inv_gut, a_inv_gut, a_inv_gut]), model.a_inv_error_from_exp())
		# Propagate model back to initial values, return chi^2
		log_mass, a_inv = model.propagate(0, return_masses=True, prop_method=MO.INT)
		return model.get_chi2_mz(log_mass, a_inv)
	
def simplex_min_mssm(initial_guesses): 
	result = SO.fmin(f, initial_guesses, full_output=1)
	print result
	return result

def mc_mssm_simplex(guess_range, **kwargs):
	return ME.mc_minimise(simplex_min_mssm, guess_range, num_trials=100, **kwargs)

#mc_mssm_simplex([[2.5,4.3], [15.5, 16.1], [24.4, 28.2]])

#[3.36825754, 16.07515093, 25.59902942]
#print f([3.37, 16.07, 25.60])
#print simplex_min_mssm([3.39232668, 16.09896961, 25.59925585])
#print simplex_min_mssm([3.39, 16.10, 25.60])

# Chi^2 + 1 Determination
#def g(x): return abs(f(x)-1)
#print [3.36983676, 16.07484974, 25.60088492]
#vary_msusy = ME.minimise_function(g, [N.arange(3.340,3.400,0.005), 16.07484974, 25.60088492])
#vary_mgut = ME.minimise_function(g, [3.36983676, N.arange(16.055, 16.100, 0.002), 25.60088492])
#vary_ainv = ME.minimise_function(g, [3.36983676, 16.07484974, N.arange(25.56, 25.64, 0.005)])

def h(x): return N.log10(f([x, 16.07515093, 25.59902942]))
#P.plot_function(h, N.arange(3.35,3.39,0.00001))
#P.plot_function(h, N.arange(3.3682,3.3683,0.0000001))

def mssm_3d_chi2(log_m_susy, log_m_gut, a_inv_gut, **kwargs):
	def g(x): return N.log10(f(x))
	P.parameter_3d(g, [log_m_susy, log_m_gut, a_inv_gut], **kwargs)

#mssm_3d_chi2(N.arange(3.34,3.38,0.002), N.arange(16.06, 16.10, 0.002), 25.59, labels=[r"$\log(M_\mathrm{SUSY})$ [GeV]", r"$\log(M_\mathrm{GUT})$ [GeV]", r"$\log(\chi^2)$"])

def mssm_contour_chi2(log_m_susy, log_m_gut, a_inv_gut, **kwargs):
	def g(x): return N.log10(f(x))
	P.parameter_contour(g, [log_m_susy, log_m_gut, a_inv_gut], **kwargs)

#mssm_contour_chi2(N.arange(3.340,3.400,0.005), N.arange(16.055, 16.100, 0.002), 25.60088492, labels=[r"$\log(M_\mathrm{SUSY})$ [log(GeV)]", r"$\log(M_\mathrm{GUT})$ [log(GeV)]", r"$\log(\chi^2)$"], contours=[0])
#mssm_contour_chi2(N.arange(3.25,3.50,0.01), 16.07484974, N.arange(25.4, 25.8, 0.02), labels=[r"$\log(M_\mathrm{SUSY})$ [log(GeV)]", r"$\alpha_\mathrm{GUT}^{-1}$", r"$\log(\chi^2)$"], contours=[0])
mssm_contour_chi2(3.36983676, N.arange(16.050, 16.100, 0.005), N.arange(25.56, 25.64, 0.005), labels=[r"$\log(M_\mathrm{GUT})$ [log(GeV)]", r"$\alpha_\mathrm{GUT}^{-1}$", r"$\log(\chi^2)$"], contours=[0])

#mssm_contour_chi2(N.arange(3.3,3.5,0.02), N.arange(15.90, 16.10, 0.02), 25.59, labels=[r"$\log(M_\mathrm{SUSY})$ [GeV]", r"$\log(M_\mathrm{GUT})$ [GeV]", r"$\log(\chi^2)$"])

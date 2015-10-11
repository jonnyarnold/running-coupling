'''Extra Dimension plots and analysis
Jonny Arnold'''

import models as MO
import plots as P
import methods as ME

import numpy as N
import scipy.optimize as SO

### COUPLING PLOTS ###

def ed_coupling(model, mass_to=20, plot_crossing=True):
	log_mass, a_inv = model.propagate(mass_to, return_masses=True)
	err_lo, err_hi = model.get_error_lines(mass_to)
	plot = P.CouplingPlot(log_mass, a_inv, err_hi, err_lo)
	if isinstance(model, MO.MSSMExtraDim): plot.plotVLine(model.log_m_susy, color='k', ls=':')
	plot.plotVLine(model.log_kkor, color='k', ls=':')
	if plot_crossing: plot.addCrossingZoom()
	plot.show()

#ed_coupling(MO.MSSMExtraDim(loops=1, extra_dims=1, log_m_susy=N.log10(MO.mz), log_kkor=5))

#ed_coupling(MO.MSSMExtraDim(loops=1, extra_dims=1, log_m_susy=N.log10(MO.mz), log_kkor=15.85))
#ed_coupling(MO.MSSMExtraDim(loops=1, extra_dims=1, log_m_susy=3.38, log_kkor=12.53))

### OTHER DDG PLOTS ###

def dim_plot(dim_range, log_kkor=5):
	m_guts, a_invs = [], []
	for d in dim_range:
		# Determine M_GUT and a_inv_gut
		model = MO.MSSMExtraDim(loops=1, log_m_susy=N.log10(MO.mz), extra_dims=d, log_kkor=log_kkor)
		log_mass, a_inv = model.propagate(20, return_masses=True, prop_method=MO.INT)
		ci, cm, cc = model.findUnification(log_mass, a_inv)
		a_inv_gut = sum(cc)/3
		m_gut = sum(cm)/3
		m_guts.append(m_gut)
		a_invs.append(a_inv_gut)
	plot = P.Plot()
	plot.plotData(m_guts, a_invs)
	plot.show()

#dim_plot(range(1,11))

#colors = [(1,0,0), (0,1,0), (0,0,1), (1,0,1), (0,1,1), (1,0.5,0), (0,1,0.5), (0.5,0,1), (1,0,0.5), (0.5,1,0), (0,0.5,1)]
def dimension_dependence_plot(model, dim_range, log_kkor_range):
	plot = P.Plot()
	
	for d in dim_range:
		print d
		mguts, ainvs, log_kkors = [], [], log_kkor_range
		min_mgut, min_kkor, min_ainv, min_chi2 = float(N.inf), float(N.inf), float(N.inf), float(N.inf)
		for k in log_kkor_range:
			# Get the couplings
			model = MO.MSSMExtraDim(loops=1, log_m_susy=N.log10(MO.mz), extra_dims=d, log_kkor=k)
			log_mass, a_inv = model.propagate(20, steps=10000, return_masses=True)#, prop_method=MO.INT)
			# Determine unification parameters
			ci, cm, cc = model.findUnification(log_mass, a_inv)
			a_inv_gut = sum(cc)/3
			m_gut = sum(cm)/3
			mguts.append(m_gut)
			ainvs.append(a_inv_gut)
			# Check if it's the minimum
			chi2 = model.get_chi2_uni(log_mass, a_inv)
			if chi2 < min_chi2:
				min_mgut, min_kkor, min_ainv, min_chi2 = m_gut, k, a_inv_gut, chi2
		#plot.plotData(mguts, log_kkors, color=colors[d-1])
		#plot.plotData([min_mgut], [min_kkor], marker='o', label=str(d), color=colors[d-1])
		plot.plotData(mguts, ainvs) #, color=colors[d-1])
		#plot.plotData([min_mgut], [min_ainv], marker='o', label=str(d), color=colors[d-1])
		
	#plot.setLabels(r"$\log(M_\mathrm{GUT})$ [log(GeV)]", r"$\log(\mu_0)$ [log(GeV)]")
	plot.setLabels(r"$\log(M_\mathrm{GUT})$ [log(GeV)]", r"$\alpha_\mathrm{GUT}^{-1}$")
	#plot.axes.legend(loc=1)
	plot.show()
	
#dimension_dependence_plot(MO.MSSMExtraDim(loops=1, log_m_susy=MO.mz), range(1,6), N.arange(15.7,15.9, 0.005))		
dimension_dependence_plot(MO.MSSMExtraDim(loops=1, log_m_susy=N.log10(MO.mz)), [4], N.arange(3,17, 1))		
#dimension_dependence_plot(MO.MSSMExtraDim(loops=1, log_m_susy=MO.mz), range(1,16), N.arang(15.5,15.9, 0.02))		

def dim_dep_dist(dim_range, log_kkor=5):
	dims, mguts = dim_range, []
	
	for d in dims:
		model = MO.MSSMExtraDim(loops=1, log_m_susy=N.log10(MO.mz), extra_dims=d, log_kkor=log_kkor)
		log_mass, a_inv = model.propagate(20, return_masses=True, prop_method=MO.INT)
		# Determine unification parameters
		ci, cm, cc = model.findUnification(log_mass, a_inv)
		m_gut = sum(cm)/3
		mguts.append(m_gut)
	
	print dims
	print mguts
	
	plot = P.Plot()
	plot.plotData(dims, mguts)
	plot.setLabels(r"$\delta$", r"$\log(M_\mathrm{GUT})$ [log(GeV)]")
	plot.show()

#dim_dep_dist(range(1,11))

### CHI^2 PLOTS ###

def plot_ed_chi2(model, log_kkor_range):
	'''Plot chi^2 against log_kkor for a model'''
	# Define Chi^2 function
	def chi2(x):
		model.log_kkor = x
		log_mass, a_inv = model.propagate(20, return_masses=True)
		return model.get_chi2_uni(log_mass, a_inv)
	
	P.plot_function(chi2, log_kkor_range)

#plot_ed_chi2(MO.MSSMExtraDim(loops=1, extra_dims=1, log_m_susy = N.log10(MO.mz)), log_kkor_range = N.arange(15.5,16.0,0.01))



### CHI^2 ANALYSIS TOOLS ###

def ed_mc_min(model, log_m_susy_range, log_kkor_range, log_m_gut_range, a_inv_gut_range, **kwargs):
	# Minimising function for MC algorithm
	def f(x):
		# Set up model
		log_m_susy, log_kkor, log_m_gut, a_inv_gut = x[0], x[1], x[2], x[3]
		model.log_m_susy, model.log_kkor = log_m_susy, log_kkor
		model.set_init(log_m_gut, N.array([a_inv_gut, a_inv_gut, a_inv_gut]), model.a_inv_error_from_exp())
		# Propagate model back to initial values, return chi^2
		log_mass, a_inv = model.propagate(0, return_masses=True)
		return model.get_chi2_mz(log_mass, a_inv)
	
	# Do the minimisation!
	ME.mc_minimise(f, [log_m_susy_range, log_kkor_range, log_m_gut_range, a_inv_gut_range], **kwargs)

# Minimise the extra dimensions model!
#ed_mc_min(MO.MSSMExtraDim(extra_dims=1), [3.3, 3.5], [12.4,12.6], [12.0,14.0], [32,34])

def ed_simplex_min(model, init_guesses):
	def f(x):
		# Set up model
		log_kkor, log_m_gut, a_inv_gut = x[0], x[1], x[2]
		model.log_kkor = log_kkor
		model.set_init(log_m_gut, N.array([a_inv_gut, a_inv_gut, a_inv_gut]), model.a_inv_error_from_exp())
		# Propagate model back to initial values, return chi^2
		log_mass, a_inv = model.propagate(0, return_masses=True)
		return model.get_chi2_mz(log_mass, a_inv)
	return SO.fmin_powell(f,init_guesses)

#print ed_simplex_min(MO.MSSMExtraDim(extra_dims=1), [12.5, 13.0, 33.6])

def ed_min_chi2(model, log_m_susy, log_kkor, plot=False, **kwargs):
	'''Finds the minimum of Chi^2 for a range/value of M_SUSY and R^-1'''
	
	def chi2(x):
		model.log_m_susy, model.log_kkor = x[0], x[1]
		log_mass, a_inv = model.propagate(20, return_masses=True)
		return model.get_chi2_uni(log_mass, a_inv)
	
	# Minimise!
	min_chi2, min_x = ME.minimise_function(chi2, [log_m_susy, log_kkor], **kwargs)
	print "X^2 = %s \t(%s)" % (min_chi2, min_x)
	
	# Plot, if required
	if plot:
		model.log_m_susy, model.log_kkor = min_x[0], min_x[1]
		ed_coupling(model)

#ed_min_chi2(MO.SMExtraDim(loops=1), [0], N.arange(5,20,1), plot=True)

def ed_chi2_err(model, log_m_susy, log_kkor, **kwargs):
	'''Find the minimum of |Chi^2 - 1| to determine error in parameters'''
	
	def chi2_min1(x):
		model.log_m_susy, model.log_kkor = x[0], x[1]
		log_mass, a_inv = model.propagate(20, return_masses=True)
		return abs(model.get_chi2_uni(log_mass, a_inv) - 1)
	
	# Minimise!
	min_chi2, min_x = ME.minimise_function(chi2_min1, [log_m_susy, log_kkor], **kwargs)
	print "X^2 = %s \t(%s)" % (min_chi2, min_x)

# Determine the optimum value for log_kkor and its error
#ed_min_chi2(MO.MSSMExtraDim(loops=1, extra_dims=2), log_m_susy = N.arange(3.35,3.45,0.01), log_kkor = N.arange(12.5,12.56,0.01), verbose=True, plot=True)
#ed_chi2_err(MO.MSSMExtraDim(loops=1, extra_dims=2), log_m_susy = N.arange(3.3,3.5,0.01), log_kkor = 12.53)

### 3D PLOTS ###

def ed_mssm_chi2_surface(model, log_m_susy_range, log_kkor_range):
	'''Plot a surface showing chi^2 as it varies for log_m_susy and log_kkor'''
	def chi2(x):
		model.log_m_susy, model.log_kkor = x[0], x[1]
		log_mass, a_inv = model.propagate(20, return_masses=True)
		return model.get_chi2_uni(log_mass, a_inv)
	
	P.parameter_3d(chi2, [log_m_susy_range, log_kkor_range], [r"$\log(M_\mathrm{SUSY})$ [log(GeV)]", r"$\log(\mu_0)$ [log(GeV)]", r"$\chi^2 (M_Z)$"])

# Plot how chi^2 varies in the ED-MSSM
#ed_mssm_chi2_surface(MO.MSSMExtraDim(extra_dims=1), N.arange(2.5, 5.0, 0.1), N.arange(8,20,0.5))

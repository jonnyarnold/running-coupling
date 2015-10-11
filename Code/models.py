'''MODELS.PY: Models used for determining running coupling
Jonny Arnold

INPUT PARAMETERS:
* M_Z: Mass of Z boson [Experimental]
* a^-1(M_Z): Inverse EM coupling strength at M_Z [Experimental]
* sin^2 W (M_Z): Weak mixing angle at M_Z [Experimental]
* a_s(M_Z): Strong coupling constant at M_Z [Experimental]
* # flavours and # Higgs doublets in model [Model]
* Beta coefficient formulae for model [Model]
* b_ij coefficients [Model]

Data from S. Eidelman et al. [Particle Data Group], Phys. Lett. B 592, 1 (2004).
'''

from __future__ import division
from math import *
from scipy.integrate import odeint
from scipy.special import gamma
import numpy as N

### EXPERIMENTAL DATA ###
mz = 91.19 # GeV
a_inv_em_mz, a_inv_em_error = 127.906, 0.019 # EM inverse coupling strength at M_Z
sin2_w_mz, sin2_w_error = 0.2312, 0.0002 # Weak mixing angle at M_Z
a_s_mz, a_s_error = 0.1187, 0.002 # Strong coupling strength at M_Z
#########################

### CONSTANTS ###
# Propagation methods
TAYLOR, TAYLOR_2DIFF, INT = 1,2,3
#################

def sign(x):
	'''Determine the sign of a number, taking 0 to be positive'''
	if x >= 0: return 1
	else: return -1

### MODULE CLASSES ###

class UnificationModel:
	'''Base class for unification models'''
	
	def __init__(self, loops, log_M_init=None, a_inv_init=None, a_inv_init_error=None, prop_method = TAYLOR):
		'''Two ways of declaring:
		UnificationModel(loops): Defines a^-1 at M_Z using default experimental parameters, defined above
		UnificationModel(loops, M_init, a_inv_init [, a_inv_init_error]): Define a^-1 at M_init directly (a_inv_init is an array)'''
		if log_M_init == None or a_inv_init == None:
			self.init_mass = log10(mz)
			self.init_a_inv = self.a_inv_from_exp(a_inv_em_mz, sin2_w_mz, a_s_mz)
			self.init_error = self.a_inv_error_from_exp((a_inv_em_mz, sin2_w_mz, a_s_mz), (a_inv_em_error, sin2_w_error, a_s_error))
		else:
			self.init_error = []
			self.set_init(log_M_init, a_inv_init, a_inv_init_error)
		
		# Set propagation method
		self.prop_method = prop_method
		
		# Set beta function
		self.beta_function = self.normal_beta_function
		
		# Calculate two-loop coefficients, if necessary
		self.loops = loops
		if loops != 1 and loops != 2:
			raise Exception('Beta function not defined for %s loops' % loops)
	
	def set_init(self, log_M_init, a_inv_init, a_inv_init_error = None):
		self.init_mass, self.init_a_inv = log_M_init, N.array(a_inv_init, N.float)
		if self.init_error == None or a_inv_init_error != None: self.init_error =  N.array(a_inv_init_error, N.float)
	
	def a_inv_from_exp(self, a_em = a_inv_em_mz, sin2_w = sin2_w_mz, a_s = a_s_mz):
		'''Set the initial values of a^-1(M_Z) based on experimental parameters (Dimensional Reduction scheme)'''
		init_a_inv = N.zeros(3)
		init_a_inv[0] = (3/5)*((cos(asin(sqrt(sin2_w))))**2)*a_em
		init_a_inv[1] = sin2_w*a_em - 2/(12*pi)
		init_a_inv[2] = 1/a_s - 3/(12*pi)
		return N.array(init_a_inv, N.float)
	
	def a_inv_error_from_exp(self, exp_params = [a_inv_em_mz, sin2_w_mz, a_s_mz], exp_error = [a_inv_em_error, sin2_w_error, a_s_error]):
		'''Set error on a^-1(M_Z) based of experimental parameters'''
		init_error = N.zeros(3, N.float)
		a_em, sin2_w, a_s = exp_params[0], exp_params[1], exp_params[2]
		a_em_error, sin2_error, a_s_error = exp_error[0], exp_error[1], exp_error[2]
		init_a_inv = self.a_inv_from_exp(a_em, sin2_w, a_s)
		# Error formulae found in appendices
		init_error[0] = sqrt(((3*a_em*sin2_error)/(5*sqrt(sin2_w)))**2 + ((init_a_inv[0]*a_em_error)/a_em)**2)
		init_error[1] = init_a_inv[1] * sqrt((sin2_error/sin2_w)**2 + (a_em_error/a_em)**2)
		init_error[2] = init_a_inv[2] * a_s_error / a_s
		return N.array(init_error, N.float)
	
	def normal_beta_function(self, a_inv, ln_mass):
		if self.loops == 1: return -(self.b_i)/(2*pi)
		else:
			# Build summation matrices
			b_ij_sum = [0.0,0.0,0.0]
			for i in (0,1,2):
				b_ij_partsum = 0.0
				for j in (0,1,2):
					b_ij_partsum += self.b_ij[i,j] / a_inv[j]
				b_ij_sum[i] = b_ij_partsum
			b_ij_sum = N.array(b_ij_sum)
			return -(1/(2*pi))*(self.b_i + (1/(4*pi))*b_ij_sum)
	
	def get_mass_range(self, mass_to, steps=1000):
		stepsize = (mass_to - self.init_mass)/steps
		return N.arange(self.init_mass, mass_to, stepsize)
	
	def propagate(self, mass_to, steps=10000, prop_method = None, return_masses = False):
		'''Propagate from log(initial mass) to a different log(mass) using beta-functions'''
		if prop_method == None: prop_method = self.prop_method
		# Get range
		log_masses = self.get_mass_range(mass_to, steps)
		ln_masses = N.log(10**log_masses)
		if prop_method == TAYLOR:
			# Taylor method
			a_inv = [self.init_a_inv]
			for m in range(0, len(ln_masses)-1):
				stepsize = ln_masses[m+1] - ln_masses[m]
				a_inv_next = a_inv[m] + self.beta_function(N.array(a_inv[m]),ln_masses[m])*stepsize
				a_inv.append(a_inv_next)
			a_inv = N.array(a_inv)
		elif prop_method == INT:
			# Integrator method
			a_inv = odeint(self.beta_function, self.init_a_inv, ln_masses, rtol=1e-10, atol=1e-10, mxhnil=1, mxstep=32768)
		else:
			raise Exception('Propagation type not understood')
		if return_masses: return log_masses, a_inv
		else: return a_inv
	
	def findUnification(self, mass, couplings):
		cross_index = []
		cross_mass = []
		cross_coupling = []
		m = 0
		mass_steps = len(mass)
		while len(cross_index) < 3 and m < mass_steps-1:
			m += 1
			if sign(abs(couplings[m,0]) - abs(couplings[m,1])) != sign(abs(couplings[m-1,0]) - abs(couplings[m-1,1])): 
				cross_index.append(m)
				cross_mass.append((mass[m] + mass[m-1])/2)
				cross_coupling.append((couplings[m,0] + couplings[m-1,0])/2)
			if sign(abs(couplings[m,0]) - abs(couplings[m,2])) != sign(abs(couplings[m-1,0]) - abs(couplings[m-1,2])): 
				cross_index.append(m)
				cross_mass.append((mass[m] + mass[m-1])/2)
				cross_coupling.append((couplings[m,0] + couplings[m-1,0])/2)
			if sign(abs(couplings[m,1]) - abs(couplings[m,2])) != sign(abs(couplings[m-1,1]) - abs(couplings[m-1,2])): 
				cross_index.append(m)
				cross_mass.append((mass[m] + mass[m-1])/2)
				cross_coupling.append((couplings[m,1] + couplings[m-1,1])/2)
		
		if len(cross_index) != 3: return [float('nan')], [float('nan')], [float('nan')]
		else: return cross_index, cross_mass, cross_coupling
	
	def get_chi2_uni(self, mass, couplings):
		ci, cm, cc = self.findUnification(mass, couplings)
		if len(ci) != 3: return N.inf
		else:
			cc_avg = sum(cc)/3
			chi2_min = N.inf
			for m in range(min(ci), max(ci)+1):
				chi2 = N.sum(((couplings[m] - cc_avg)/self.init_error)**2)
				if chi2 < chi2_min: chi2_min = chi2
			return chi2
			
	def get_error_lines(self, mass_to, steps=10000, prop_type=None):
		'''Get line corresponding to the upper and lower errors in a^-1'''
		self.init_a_inv = self.init_a_inv + self.init_error
		a_inv_hi = self.propagate(mass_to, steps, prop_type)
		self.init_a_inv = self.init_a_inv - 2*self.init_error
		a_inv_lo = self.propagate(mass_to, steps, prop_type)
		self.init_a_inv = self.init_a_inv + self.init_error
		return a_inv_lo, a_inv_hi
	
	def get_chi_mz(self, log_mass, a_inv):
		'''Return a list of Chi values at M_Z, in comparison with experimental data'''		
		# Get the experimental values at M_Z
		a_exp, a_exp_error = self.a_inv_from_exp(), self.a_inv_error_from_exp()
		# Find M_Z, interpolate if necessary	
		log_mz = log10(mz)
		a_inv_mz = []
		for i in range(0, len(log_mass)-1):
			if log_mass[i] == log_mz: 
				a_inv_mz = a_inv[i]
				break
			elif log_mass[i+1] == log_mz:
				a_inv_mz = a_inv[i+1]
				break
			elif (log_mass[i] > log_mz and log_mass[i+1] < log_mz) or (log_mass[i] < log_mz and log_mass[i+1] > log_mz):
				# Interpolate!
				a_inv_mz = ((log_mz - log_mass[i])*a_inv[i+1] + (log_mass[i+1] - log_mz)*a_inv[i])*(1/(log_mass[i+1] - log_mass[i]))
				break
		if a_inv_mz == []:
			print "Warning! Chi^2 not determinable"
			a_inv_mz = N.zeros(len(a_inv))
			a_inv_mz.fill(N.inf)
			return a_inv_mz
		else: return abs((a_inv_mz - a_exp) / (a_exp_error))
	
	def get_chi2_mz(self, log_mass, a_inv):
		'''Return the sum of Chi^2 at M_Z, in comparison with experimental data'''
		return N.sum(self.get_chi_mz(log_mass, a_inv)**2)

class StandardModel(UnificationModel):
	'''Unification model with Standard Model parameters and beta coefficients'''
	
	def __init__(self, loops = 2, flavours = 3, higgs_doublets = 1, M_init=None, a_inv_init=None, a_inv_init_error=None, prop_method = TAYLOR):
		UnificationModel.__init__(self, loops, M_init, a_inv_init, a_inv_init_error, prop_method)
		if loops == 1:
			self.title = 'Standard Model (1 loop)'
		else:
			self.title = 'Standard Model, (2 loops)'
		self.set_coefficients(loops, flavours, higgs_doublets)
	
	def set_coefficients(self, loops, f = 3, h = 1):
		'''Set one loop (and two-loop, if required) coefficients'''
		# Calculate one-loop coefficients
		self.b_i = N.array([0.0,0.0,0.0])
		self.b_i[0] = (4/3)*f + (1/10)*h
		self.b_i[1] = (-22/3) + (4/3)*f + (1/6)*h
		self.b_i[2] = -11 + (4/3)*f
		
		# Calculate two-loop coefficients, if necessary
		if loops == 2:
			self.b_ij = N.array([[(19/15)*f+(9/50)*h, (3/5)*f+(9/10)*h, (44/15)*f], [(1/5)*f+(3/10)*h, -(136/3)+(49/3)*f+(13/6)*h, 4*f], [(11/30)*f, (3/2)*f, -102+(76/3)*f]], N.float)

class MSSMModel(UnificationModel):
	'''Unification model with MSSM parameters and beta coefficients'''
	
	def __init__(self, loops=2, log_m_susy = 0.0, flavours = 3, higgs_doublets = 2, log_m_init=None, a_inv_init=None, a_inv_init_error=None, prop_method = TAYLOR):
		UnificationModel.__init__(self, loops, log_m_init, a_inv_init, a_inv_init_error, prop_method)
		self.set_log_m_susy(log_m_susy)
		self.loops = loops
		if loops == 1:
			self.title = 'MSSM (1 loop'
		else:
			self.title = 'MSSM (2 loops'
		if log_m_susy == 0:
			self.title += ')'
		else:
			self.title += ', log(M_SUSY) = %s)' % self.log_m_susy
		
		self.set_coefficients(loops, flavours, higgs_doublets)
		
		self.sm = StandardModel(loops)
		self.beta_function = self.scaled_beta_function
	
	def set_log_m_susy(self, log_m_susy):
		self.log_m_susy = log_m_susy
		self.ln_m_susy = N.log(10**log_m_susy)
	
	def set_coefficients(self, loops, f = 3, h = 2):
		# Calculate one-loop coefficients
		self.b_i = N.array([0.0,0.0,0.0])
		self.b_i[0] = 2*f + (3/10)*h
		self.b_i[1] = -6 + 2*f + (1/2)*h
		self.b_i[2] = -9 + 2*f
		
		# Calculate two-loop coefficients, if necessary
		if loops == 2:
			self.b_ij = N.array([[(38/15)*f+(9/50)*h, (6/5)*f+(9/10)*h, (88/15)*f], [(2/5)*f+(3/10)*h, -24+14*f+(7/2)*h, 8*f], [(11/15)*f, 3*f, -54+(68/3)*f]], N.float)
	
	def scaled_beta_function(self, a_inv, ln_mass):
		if ln_mass < self.ln_m_susy:
			# Standard model
			return self.sm.beta_function(a_inv, ln_mass)
		else:
			# MSSM
			return self.normal_beta_function(a_inv, ln_mass)

class ExtraDimModel(UnificationModel):
	'''Base class for extra-dimensional models'''
	def __init__(self, model, log_kkor=5, extra_dims=0, kk_fermions=0):
		# kkor = Kaluza-Klein Orbifold Radius (mu_0)
		# extra_dims = # of extra dimensions (delta)
		# kk_fermions = # of fermion generations that admit Kaluza-Klein excitations
		self.model = model
		self.log_kkor, self.extra_dims, self.kk_fermions = log_kkor, extra_dims, kk_fermions
		#self.beta_function = self.extra_dim_beta_function
		
		# New beta coefficients
		self.set_beta_coefficients()
	
	def __getattr__(self, name): return getattr(self.model, name)
	
	def set_beta_coefficients(self, kk_fermions=None):
		pass
	
	def X(self, extra_dims): return (N.pi**(extra_dims/2.0))/(gamma(1 + extra_dims/2.0))
	
	def propagate(self, mass_to, steps=10000, prop_method = None, return_masses = False):
		'''Propagate from log(initial mass) to a different log(mass)'''
		# Get range
		log_masses = self.get_mass_range(mass_to, steps)
		ln_masses = N.log(10**log_masses)
		a_inv = []
		for m in ln_masses: a_inv.append(self.get_a_inv(m))
		a_inv = N.array(a_inv)
		
		if return_masses: return log_masses, a_inv
		else: return a_inv
	
	def get_a_inv(self, ln_mass):
		mass = N.exp(ln_mass)
		if N.log10(N.exp(ln_mass)) < self.log_kkor: 
			return self.init_a_inv - self.b_i / (2 * N.pi) * N.log(mass / 10**self.init_mass)
		else:
			result = self.init_a_inv - self.b_i / (2 * N.pi) * N.log(mass / 10**self.init_mass) + self.b_i_tilde / (2* N.pi) * N.log(mass / 10**self.log_kkor) 
			if self.extra_dims > 0: result -= (self.b_i_tilde * self.X(self.extra_dims)) / (2 * N.pi * self.extra_dims) * ((mass/10**self.log_kkor)**self.extra_dims - 1)
			return result
	
	#def extra_dim_beta_function(self, a_inv, ln_mass):
	#	mass = N.exp(ln_mass)
	#	if N.log10(mass) < self.log_kkor: 
	#		return self.model.beta_function(a_inv, ln_mass)
	#	else:
	#		return (self.b_i_tilde - self.b_i) / (2 * N.pi) - (self.b_i_tilde * self.X(self.extra_dims) * self.extra_dims) / (2 * N.pi) * N.exp(self.extra_dims * (ln_mass - N.log(10**self.log_kkor)))

class MSSMExtraDim(ExtraDimModel):
	'''One-loop extra-dimensional model using MSSM'''
	
	def __init__(self, log_kkor = 5, log_m_susy = 3, extra_dims = 0, kk_fermions = 0, **kwargs):
		ExtraDimModel.__init__(self, MSSMModel(log_m_susy=log_m_susy, **kwargs), log_kkor, extra_dims, kk_fermions)
		self.log_m_susy = log_m_susy
	
	def set_beta_coefficients(self, kk_fermions=None):
		if kk_fermions == None: kk_fermions = self.kk_fermions
		self.b_i_sm = self.model.sm.b_i
		self.b_i_mssm = self.model.b_i
		
		self.b_i_tilde = N.array([0.0,0.0,0.0])
		self.b_i_tilde[0] = 0.6 + 4*kk_fermions
		self.b_i_tilde[1] = -3 + 4*kk_fermions
		self.b_i_tilde[2] = -6 + 4*kk_fermions
	
	def get_a_inv(self, ln_mass):
		mass = N.exp(ln_mass)
		log_mass = N.log10(mass)
		if log_mass <= self.log_m_susy:
			return self.init_a_inv - self.b_i_sm / (2 * N.pi) * N.log(mass / 10**self.init_mass)
		elif log_mass <= self.log_kkor: 
			result = self.init_a_inv - self.b_i_sm / (2 * N.pi) * N.log(10**self.log_m_susy / 10**self.init_mass)
			result -= self.b_i_mssm / (2 * N.pi) * N.log(mass / 10**self.log_m_susy)
			return result
		else:
			result = self.init_a_inv - self.b_i_sm / (2 * N.pi) * N.log(10**self.log_m_susy / 10**self.init_mass)
			if self.extra_dims > 0: 
				result -= self.b_i_mssm / (2 * N.pi) * N.log(10**self.log_kkor / 10**self.log_m_susy)
				result += self.b_i_tilde / (2* N.pi) * N.log(mass / 10**self.log_kkor) - (self.b_i_tilde * self.X(self.extra_dims)) / (2 * N.pi * self.extra_dims) * ((mass/10**self.log_kkor)**self.extra_dims - 1)
			else: 
				result -= self.b_i_mssm / (2 * N.pi) * N.log(mass / 10**self.log_m_susy)
			return result

######################

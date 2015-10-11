'''METHODS: Useful routines for gathering data about models
Jonny Arnold'''

import plots as P

import copy
import random as R

import numpy as N

def mc_minimise(f, x, num_trials=N.inf, verbose=True):
	'''Use a Monte-Carlo method to minimise f, subject to the ranges in x'''
	
	if(num_trials==N.inf and verbose==False): exit(msg="Fatal Error: No output specified for mc_minimise")
	
	# Set up minimum values
	f_min, x_min = N.inf, []
	
	# Start trials
	trial = 0
	while trial < num_trials:
		trial += 1
		
		# Set up initial values
		x_i = []
		for j in range(len(x)):
			x_i.append(R.uniform(min(x[j]), max(x[j])))
		
		# Evaluate function at this point
		f_i = f(x_i)
		
		# Determine if new minimum has been reached
		if f_i < f_min:
			f_min, x_min = f_i, x_i
			if verbose: print "#%s: f=%s\tx=%s" % (trial, f_min, x_min)
	
	# Iterations over: return results
	return f_min, x_min
			

def minimise_function(f, x, verbose=False):
	'''Minimise the function f, given the parameters x (one or more can be a list). Returns the minimum value of f and the value of x at which it occurs.'''
	# Determine the lists
	list_ids = P.get_lists(x)
	if list_ids == []: return f(x), x
	else:
		# Set up the indices used during iteration
		list_indices = []
		for i in range(len(list_ids)): list_indices.append(0)
		
		# Set up minimum tracker
		min_f, min_x = N.inf, []
		
		# Start iterating
		done_iterating = False
		while not done_iterating:
			x_i = copy.copy(x)
			# Set list parameters
			for i in range(len(list_ids)): x_i[list_ids[i]] = x[list_ids[i]][list_indices[i]]
			# Evaluate function
			f_i = f(x_i)
			if verbose: print "%s\t(%s)" % (f_i, x_i)
			# Check for minimum; replace if lower than current minimum
			if f_i < min_f:
				min_f, min_x = f_i, x_i
		
			# Increment list indices
			j = len(list_indices) - 1
			done_incrementing = False
			while j >= 0 and not done_incrementing:
				list_indices[j] += 1
				if list_indices[j] >= len(x[list_ids[j]]):
					list_indices[j] = 0
					j -= 1
				else: done_incrementing = True
			if j < 0: done_iterating = True
		
		# Finished iteration: output minimum value
		return min_f, min_x
	


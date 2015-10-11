'''PLOTS: Graph plotting class
Jonny Arnold'''

import numpy as N
import pylab as P
from matplotlib import rc, cm
rc('mathtext', fontset='stixsans')
rc('font', size=18)
from mpl_toolkits.mplot3d import axes3d, Axes3D

from math import *
import copy

def sign(x):
	if x >= 0: return 1
	else: return -1
	
def get_lists(x):
	'''Returns the indices corresponding to lists in the given list'''
	list_indices = []
	for i in range(0, len(x)):
		if isinstance(x[i], list) or (type(x[i]) == type(N.array([]))):
			list_indices.append(i)
	return list_indices

### 2D PLOTS ###

class Plot:
	
	def __init__(self):
		self.fig = P.figure()
		self.axes = self.fig.gca()
	
	def show(self): P.show()
	
	# Plotting
	def plotData(self, x, y, **kwargs): self.axes.plot(x, y, **kwargs)
	def plotErrors(self, x, yerr_hi, yerr_lo, **kwargs): self.axes.fill_between(x, yerr_hi, yerr_lo, **kwargs)
	
	def plotVLine(self, x, label=None, **kwargs): 
		self.axes.axvline(x=x, **kwargs)
		if label != None: self.axes.text(x, 0, label)
	
	def plotHLine(self, y, label=None, **kwargs):
		self.axes.axvline(y=y, **kwargs)
		if label != None: self.axes.text(0, y, label)
	
	# Titles, labels, ticks and legends
	def setTitle(self, title): self.axes.set_title(title)
	def setXLabel(self, xlabel): self.axes.set_xlabel(xlabel)
	def setYLabel(self, ylabel): self.axes.set_ylabel(ylabel)
	def setLabels(self, xlabel, ylabel):
		self.axes.set_xlabel(xlabel)
		self.axes.set_ylabel(ylabel)
	def setTicks(self, xticks, yticks):
		self.axes.set_xticks(xticks)
		self.axes.set_yticks(yticks)
		self.axes.set_xlim(min(xticks), max(xticks))
		self.axes.set_ylim(min(yticks), max(yticks))
	def setLegend(self, showlegend=True, labels=None, location=0):
		pass
	
	# Subfigure
	def addSubfig(self, pos, size): 
		self.axes = self.fig.add_axes([pos[0], pos[1], size[0], size[1]])
	
class CouplingPlot(Plot):
	# Class constants
	cols = [(1,0,0,0), (0,1,0,0), (0,0,1,0)] # Graphing colours
	err_cols = [(1,0.75,0.75,0.5), (0.75,1,0.75,0.5), (0.75,0.75,1,0.5)] # Error bar colours
	zoom_pos, zoom_size = [0.375,0.62], [0.25,0.25]
	
	def __init__(self, mass, couplings, errors_hi=None, errors_lo=None):
		Plot.__init__(self)
		self.setLabels(r"$\log(\mu_0)$ [log(GeV)]", r"$\alpha^{-1}$")
		self.mass, self.couplings, self.errors = mass, couplings, [errors_hi, errors_lo]
		self.plotCoupling()
	
	def plotCoupling(self, x_range = None):
		if x_range == None: 
			ci, cm, cc = self.getCrossing()
			if len(ci) < 3: 
				# Find the least ridiculous values
				cutoff = 1
				ok = False
				while not ok and cutoff < len(self.mass)-1:
					if N.any(self.couplings[cutoff,:] > 100): ok = True
					cutoff += 1
				x_range = [0, cutoff]
			else: x_range = [0, min(len(self.mass), max(ci)+floor(0.05*max(ci)))]
		xmin, xmax = x_range[0], x_range[1]
		
		P.hold(True)
		if any(i != None for i in self.errors):
			for c in range(3): self.plotErrors(self.mass[xmin:xmax], self.errors[0][xmin:xmax,c], self.errors[1][xmin:xmax,c], color=self.err_cols[c])
		for c in range(3): self.plotData(self.mass[xmin:xmax], self.couplings[xmin:xmax,c], color=self.cols[c])
		self.axes.set_xlim(self.mass[xmin], self.mass[xmax-1])
		P.hold(False)
	
	def getCrossing(self):
		cross_index = []
		cross_mass = []
		cross_coupling = []
		m = 0
		mass_steps = len(self.mass)
		while len(cross_index) < 3 and m < mass_steps-1:
			m += 1
			if sign(abs(self.couplings[m,0]) - abs(self.couplings[m,1])) != sign(abs(self.couplings[m-1,0]) - abs(self.couplings[m-1,1])): 
				cross_index.append(m)
				cross_mass.append((self.mass[m] + self.mass[m-1])/2)
				cross_coupling.append((self.couplings[m,0] + self.couplings[m-1,0])/2)
			if sign(abs(self.couplings[m,0]) - abs(self.couplings[m,2])) != sign(abs(self.couplings[m-1,0]) - abs(self.couplings[m-1,2])): 
				cross_index.append(m)
				cross_mass.append((self.mass[m] + self.mass[m-1])/2)
				cross_coupling.append((self.couplings[m,0] + self.couplings[m-1,0])/2)
			if sign(abs(self.couplings[m,1]) - abs(self.couplings[m,2])) != sign(abs(self.couplings[m-1,1]) - abs(self.couplings[m-1,2])): 
				cross_index.append(m)
				cross_mass.append((self.mass[m] + self.mass[m-1])/2)
				cross_coupling.append((self.couplings[m,1] + self.couplings[m-1,1])/2)
		return cross_index, cross_mass, cross_coupling
	
	def addCrossingZoom(self):
		ci, cm, cc = self.getCrossing()
		
		P.hold(True)
		self.addSubfig(self.zoom_pos, self.zoom_size)
		if (max(ci) - min(ci)) > 10: width = (2*(max(ci) - min(ci)))
		else: width = 10
		#width = int(0.01*len(self.mass))
		subfig_xmin, subfig_xmax = max(0, min(ci)-width), min(len(self.mass)-1, max(ci)+width)
		
		self.plotCoupling([subfig_xmin, subfig_xmax])
		
		xticks = [ceil(self.mass[subfig_xmin]*100)/100, round(sum(cm)/3, 2), floor(self.mass[subfig_xmax-1]*100)/100]
		yticks = [ceil(min(self.couplings[subfig_xmin])*100)/100, round(sum(cc)/3, 2), floor(max(self.couplings[subfig_xmax-1])*100)/100]
		
		self.setTicks(xticks, yticks)
		P.hold(False)

def plot_function(f, x, verbose=False):
	'''Plot a function f, dependent on one range of parameters x'''
	y = []
	for x_i in x:
		y.append(f(x_i))
	plot = Plot()
	plot.plotData(x,y)
	P.show()

### CONTOUR PLOTS ###

def parameter_contour(f, x, labels = [], levels = 1000, contours=[], verbose = False):
	'''Returns a contour plot of a pair of varying parameters against a function f. Takes a list of variables, two of which should be lists of values'''
	li = get_lists(x)
	if len(li) != 2: raise Exception('Two list arguments required in parameters of parameter_graph')
	
	X, Y, Z = [], [], []
	for xi in x[li[0]]:
		templist_x = []
		templist_y = []
		templist_z = []
		for yi in x[li[1]]:
			x_now = copy.copy(x)
			x_now[li[0]] = xi
			x_now[li[1]] = yi
			f_now = f(x_now)
			if verbose: print "x = %s => f(x) = %s" % (x_now, f_now)
			templist_x.append(xi)
			templist_y.append(yi)
			templist_z.append(f_now)
		X.append(templist_x)
		Y.append(templist_y)
		Z.append(templist_z)
	X, Y, Z = N.array(X), N.array(Y), N.array(Z)
	
	#if swap_xy == False: X, Y, Z = x[li[0]], x[li[1]], []
	#else: X, Y, Z = x[li[1]], x[li[0]], []
	#for xi in X:
	#	templist = []
	#	for yi in Y:
	#		x_now = copy(x)
	#		x_now[li[0]] = xi
	#		x_now[li[1]] = yi
	#		templist.append(f(x_now))
	#	Z.append(templist)
	#Z = N.transpose(N.array(Z))
	
	fig = P.figure()
	P.hold(True)
	fig = P.contourf(X,Y,Z, 1000, antialias=True)
	if labels != []:
		if len(labels) >= 1: P.xlabel(labels[0])
		if len(labels) >= 2: P.ylabel(labels[1])
		if len(labels) >= 3:
			cb = P.colorbar()
			cb.set_label(labels[2])
	if contours:
		P.contour(X,Y,Z, contours, colours='k')
	P.hold(False)
	#P.xlim([min(X), max(X)])
	#P.ylim([min(Y), max(Y)])
	
	P.show()

### 3D PLOTS ###

def parameter_3d(f, x, labels=[], levels = 1000, verbose = False):
	'''Returns a 3D plot of a pair of varying parameters against a function f. Takes a list of variables, two of which should be lists of values'''
	li = get_lists(x)
	if len(li) != 2: raise Exception('Two list arguments required in parameters of parameter_3d')
	
	X, Y, Z = [], [], []
	for xi in x[li[0]]:
		templist_x = []
		templist_y = []
		templist_z = []
		for yi in x[li[1]]:
			x_now = copy.copy(x)
			x_now[li[0]] = xi
			x_now[li[1]] = yi
			f_now = f(x_now)
			if verbose: print "x = %s => f(x) = %s" % (x_now, f_now)
			templist_x.append(xi)
			templist_y.append(yi)
			templist_z.append(f_now)
		X.append(templist_x)
		Y.append(templist_y)
		Z.append(templist_z)
	X, Y, Z = N.array(X), N.array(Y), N.array(Z)
	fig = P.figure()
	axes = Axes3D(fig)
	axes.plot_surface(X,Y,Z,rstride=1,cstride=1,linewidth=0,cmap=cm.spectral)
	if labels != []:
		if len(labels) >= 1: P.xlabel(labels[0])
		if len(labels) >= 2: P.ylabel(labels[1])
		if len(labels) >= 3: axes.set_zlabel(labels[2])
	P.show()

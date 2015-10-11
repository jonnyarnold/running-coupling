'''SM/MSSM plots
Jonny Arnold'''

import models as MO
import methods as ME
import plots as P

import numpy as N

### COUPLING PLOTS ###

def plot_coupling(model, mass_to=20, plot_crossing=True):
	'''Show coupling plot for a MSSM Model'''
	log_mass, a_inv = model.propagate(mass_to, return_masses=True)
	err_lo, err_hi = model.get_error_lines(mass_to)
	plot = P.CouplingPlot(log_mass, a_inv, err_hi, err_lo)
	if isinstance(model, MO.MSSMModel): plot.plotVLine(model.log_m_susy, color='k', ls=':')
	if plot_crossing: plot.addCrossingZoom()
	plot.show()

# SM, Forward Propagation
#plot_coupling(MO.StandardModel(loops=1), plot_crossing=False)
#plot_coupling(MO.StandardModel(loops=2), plot_crossing=False)

# MSSM, Forward propagation
#plot_coupling(MO.MSSMModel(loops=1, log_m_susy=N.log10(MO.mz)))
#plot_coupling(MO.MSSMModel(loops=2, log_m_susy=3.4))

# MSSM, Back propagation
#model = MO.MSSMModel(loops=2, log_m_susy=3.39)
#model.set_init(log_M_init=16.09, a_inv_init=N.array([25.59,25.59,25.59]))
#plot_coupling(model, mass_to=2, plot_crossing=False)

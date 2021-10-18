from __future__ import print_function
from astropy.io import fits
import numpy as np
import math as math
import matplotlib.pyplot as plt
from pylab import *
import scipy.integrate as si
import os.path

# =================================================================================
# ================================== INPUT ========================================
# =================================================================================


directories = ['']
output_directory = ''

Tstart_all = []
Tend_all = []
highest_energy_photons = []
OPT_plot_limit_vE_highest = True #Option to choose true or false

# define model components
models = ['']

# eac true or false
eac = False

E1 = 1.  # keV   # Energy range  minimum for plotting SED
E2 = 7.  # in keV = 10 GeV   # Energy range meximum for plotting SED

color_list = [''] # define colors
component_linestyle = ':'  # define component linestyle
component_linewidth = 2. # define component linewidth

model_linestyle = '-' # define model linestyle
model_linewidth = 1. # define model linestyle

totalmodel_linewidth = 3. # define total model linewidth

OPT_LEGEND = 'model'  # choose between interval or model
OPT_SHADE = False  # shade confidence interval

# Constant value
keVtoErg = 1.60218e-9

# define plot limits:
xlim_min = 10.
xlim_max = 1.e6
ylim_min = 9.e-10
ylim_max = 1.0e-5

# change figure size
fig = plt.figure(1, figsize=(7., 4.))

if len(models) > 1:
    str_allmodels = '__'.join(models)
else:
    str_allmodels = models[0]

for im, model_complete in enumerate(models):
    print(im, color_list[im])
    color10 = color_list[im]
    directory = directories[im]
    if OPT_plot_limit_vE_highest:
        Emax_highest = highest_energy_photons[im]

    Tstart = Tstart_all[im]  # Fill value of start time of interval used for fitting.
    Tend = Tend_all[im]  # Fill value of end time of interval used for fitting.
    Tint = Tend - Tstart  # integration time. This will be time interval over which the fitting was done.
    # model_complete='CPL_PL' #String with main model name. The first model should be the base model. This string will be used to obtain the full file name. If EAC is applied and as a result is a part of the filename, it should not be added in this string. Instead the value of the next parameter eac should be set to True. This will add _EAC to the fits file name. If eac is not applied, eac should be set to False.
    str_modelname = model_complete.replace('_', '+')
    str_timeinterval = str(Tstart) + ' - ' + str(Tend) + ' s'
    if OPT_LEGEND == 'interval':
        str_label = str_timeinterval
    elif OPT_LEGEND == 'model':
        str_label = str_modelname
    with open('multiple_SEDs_loop.py') as infile:
        exec (infile.read())
    print(im)

# change legend location, x and y axes labels
plt.legend(loc='upper right')
plt.xlabel('Energy (keV)')
plt.ylabel(r'$\nu F\nu$ (erg cm$^{-2}$ s$^{-1}$)')
plt.xlim([xlim_min,xlim_max])
plt.ylim([ylim_min,ylim_max])

# change output figure name
plt.savefig(output_directory + 'SED_' + str_allmodels + '.png', bbox_inches="tight")
plt.savefig(output_directory + 'SED_' + str_allmodels + '.pdf', bbox_inches="tight")
plt.show()



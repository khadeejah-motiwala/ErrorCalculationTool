from __future__ import print_function
from astropy.io import fits
import numpy as np
import math as math
import matplotlib.pyplot as plt
from pylab import *
import scipy.integrate as si
import os.path
#from numpy.core.defchararray import add



#=================================================================================
#================================== INPUT ========================================
#=================================================================================


directories=['/home/khadeejah/error_prop/']
output_directory='/home/khadeejah/error_prop/'

Tstart_all=[3.84, 8.2, 10.8]
Tend_all=[8.2, 10.8,  17.408]
highest_energy_photons=[ 408.11203e3,  277.63062e3, 152.10425e3]
plot_limit_vE_highest=True


models=['BAND_BB','BAND_BB','BAND_BB']

colors=['navy','C0', 'mediumturquoise']
component_linestyle=':'


if len(models)>1:
    str_allmodels='__'.join(models)
else:
    str_allmodels=models[0]

eac=True

E1= 1. #keV   # Energy range  minimum for plotting SED
E2=7. # in keV = 10 GeV   # Energy range meximum for plotting SED

OPT_LEGEND='interval' # choose between interval or model
OPT_SHADE='True'    # shade confidence interval

#Constant values
keVtoErg=1.60218e-9



#Plot parameters

#yscale

xlim_min=10.
xlim_max=1.e6
xlim_max=1.e5
ylim_min=9.e-10
ylim_min=9.e-8
ylim_max=1.0e-5

fig = plt.figure(1,figsize=(7.,4.))
print (models)
for im, model_complete in enumerate(models):
    print (im, colors[im])
    color10=colors[im]
    directory=directories[im]
    Emax_highest=highest_energy_photons[im]

    Tstart=Tstart_all[im] # Fill value of start time of interval used for fitting.
    Tend=Tend_all[im] # Fill value of end time of interval used for fitting.
    Tint=Tend-Tstart #integration time. This will be time interval over which the fitting was done. 
    # model_complete='CPL_PL' #String with main model name. The first model should be the base model. This string will be used to obtain the full file name. If EAC is applied and as a result is a part of the filename, it should not be added in this string. Instead the value of the next parameter eac should be set to True. This will add _EAC to the fits file name. If eac is not applied, eac should be set to False.
    str_modelname=model_complete.replace('_','+')
    str_timeinterval=str(Tstart)+' - '+str(Tend)+' s'
    if OPT_LEGEND=='interval':
        str_label=str_timeinterval
    elif OPT_LEGEND=='model':
        str_label=str_modelname
    with open('multiple_SEDs_loop.py') as infile:
                exec(infile.read())
    print (im)
plt.legend(loc='upper right')
plt.xlabel('Energy (keV)')
plt.ylabel(r'$\nu F\nu$ (erg cm$^{-2}$ s$^{-1}$)')
plt.xlim([xlim_min,xlim_max])
plt.ylim([ylim_min,ylim_max])
plt.savefig(output_directory+'SED_'+str_allmodels+'.png',bbox_inches="tight")
plt.savefig(output_directory+'SED_'+str_allmodels+'.pdf',bbox_inches="tight")
plt.show()

            

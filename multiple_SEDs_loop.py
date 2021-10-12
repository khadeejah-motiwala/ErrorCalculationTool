
#!/usr/bin/env python


#Code prepared by Khadeejah Motiwala and Saeeda Sajjad
#2021/07/12

from __future__ import print_function
from astropy.io import fits
import numpy as np
import math as math
import matplotlib.pyplot as plt
from pylab import *
from scipy.special import comb
import scipy.integrate as si
#from numpy.core.defchararray import add




#This code is run in the following way.
#Input: fits file from Rmfit containing the results of spectral fits
#Output: SED butterfly plot with complete error propagation


#Note: This code actually only requires the values of parameters obtained from a spectral fit, the error on the parameters and the covariance matrix. These values are stored in the arrays param_mean, param_err and covariance_matrix respectively. The parameters are ordered according to the Rmfit default order (see comments below). This means that this code can be easily adapted to plot SEDs if the fit parameters, their errors and the covariance matrix is available in some other form.

#Note: To see how the fits file can be obtained from Rmfit, see comments below (one can search for Rmfit).

#Note: The Rmfit file should be named in the following way: BASEMODEL_ADDITIONALCOMPONENT1_ADDITIONALCOMPONENT2_EAC.fit. This means if the spectrum is fitted with the Band model alone the file name should be BAND.fit. If it is fitted with the Band model and an EAC factor is applied, it should be named BAND_EAC.fit. Other examples can be: BAND_BB.fit, BAND_PL.fit, BAND_PL_BB.fit, BAND_PL_BB_EAC.fit

#Note: The base models currently used are CPL, BAND, SBPL. The additional models currently used are PL, BB. More functions can be added if needed.

#Method of use:
#=============
#In the INPUT section, give the input parameters, file names etc. See comments in input section.



#List of models and their parameters in Rmfit
'''
For details look at the file photon_models.ps among rmfit files
SBPL = Smoothly Broken Power Law
6 parameters - 4 free - 2 fixed
1----Amplitude - A - vary photon/(s cm^2 keV)
2----Pivot E - Epiv - fix at 100 keV
3----Index1 - lambda1 - vary
4----Break E - Eb - vary keV
5----Break scale - Delta -  fix at 0.3 decades E
6----Index2 - lambda2 - vary

Band
4 parameters
1----Amplitude - A - vary photon/(s cm^2 keV)
2----Epeak - vary keV
3----alpha - vary
4----beta - vary

Compton
4 parameters - 3 free - 1 fixed
1----Amplitude - A - vary photon/(s cm^2 keV)
2----Epeak - vary keV
3----Index - vary
4----Pivot energy - Epiv - fix at 100 keV

Power Law
3 parameters - 2 free - 1 fixed
1----Amplitude - A - vary photon/(s cm^2 keV)
2----Pivot energy - Epiv - fix at 100 keV
3----Index - vary

Black Body
2 parameters
1----Amplitude - A - vary photon/(s cm^2 keV)
2----kT - vary (electron energy in keV)

Order in rmfit
PL SBPL Band Comptonized BlackBody

'''

#The fits file used  is the one obtained from rmfit by going to the Fit Display window then Fit Results -> Write results to file
#The code can be adapted to reading a covariance matrix and values of parameters and their errors from some other format too.

#allmodel_names=np.asarray(['SBPL','BAND','CPL','PL','BB'])
#parnumbers=[6,4,4,3,2]

#Parameter names
#For SBPL
#parnames=['A_sbpl','Epiv_sbpl','Alpha','Ebreak','Break Scale','Beta']
#Use symbol definition from Ep.png (from Feraol). It is equivalent to the definition from photon_model.ps (rmfit).

#For Band
#parnames=['A_band','Epeak_band','alpha','beta']

#For Compton
#parnames=['A_cpl','Epeak_cpl','Index_cpl','Epiv_cpl']

#For Power Law
#parnames=['A_pl','Epiv_pl','Index_pl']

#For Black Body
#parnames=['A_bb','kT']


#=================================================================================
#=========================== FUNCTION DEFINITION =================================
#=================================================================================




#Functions can be definied in a standalone file which can be called through the following two lines. These should be uncommented if needed. However, in this standalone version, all the functions are defined below, from this point, till the INPUT section. These lines would be present in the function_definitions.py file, if it was being used. The file would have to be in the same folder as this code in order to run the code.
#with open('function_definitions.py') as infile:
#    exec(infile.read())





def get_model_details(model_name):
    #order of models in rmfit
    rmfitorder=np.array(['PL', 'SBPL', 'BAND', 'CPL', 'BB', 'EAC'])
    #number of parameters in each of these models. I am putting 0 for EAC for now, but I will update this later when the fits file is read.
    rmfit_par_numbers=np.array([3, 6, 4, 4, 2, 0])
    #names of all parameters (model-wise)
    rmfit_parnames=np.array([np.array(['A_pl','Epiv_pl','Index_pl']),
                             np.array(['A_sbpl','Epiv_sbpl','alpha_sbpl','Ebreak_sbpl','Delta_sbpl','beta_sbpl']),
                             np.array(['A_band','Epeak_band','alpha_band','beta_band']),
                             np.array(['A_cpl','Epeak_cpl','Index_cpl','Epiv_cpl']),
                             np.array(['A_bb','kT_bb']),
                             np.array(['EAC'])],dtype=object)
    rmfit_parfix=np.asarray([np.asarray(['v','f','v']),
                             np.asarray(['v','f','v','v','f','v']),
                             np.asarray(['v','v','v','v']),
                             np.asarray(['v','v','v','f']),
                             np.asarray(['v','v'])],dtype=object)
    model_array = np.array(model_name.split('_')) #Get array of models in the original order
    models=np.zeros(model_array.shape[0],dtype='str') #make new empty array of the same length for rmfit order
    #print model_name
    #print 'model_array',model_array
    indices=np.where(np.in1d(rmfitorder,model_array))[0] #Get indices of the models in the original order according to rmfit order
    #print indices
    models=rmfitorder[indices] # new model array
    Nparameters=rmfit_par_numbers[indices] #Array for number of parameters according to rmfit order
    parnames=np.concatenate(rmfit_parnames[indices]) #Array of parameter names according to rmfit order.
    parfix=np.concatenate(rmfit_parfix[indices]) #Array of parameter fix or vary according to rmfit order.
    #print rmfit_parnames[indices]
    #print 'nparnames',parnames
    return models, Nparameters,parnames,parfix



###################################################################################
#FUNCTIONS WITH MODEL DEFINITIONS, ERROR CALCULATION, OTHERS
###################################################################################






def fSBPL(E,A,Epiv,alpha,Ebreak,delta,beta):
    Epiv=100.
    delta=0.3
    r=math.log10(E/Ebreak)/delta
    rp=math.log10(Epiv/Ebreak)/delta
    a=1./2.*delta*(beta-alpha)*np.log((np.exp(r)+np.exp(-r))/2.)
    ap=1./2.*delta*(beta-alpha)*np.log((np.exp(rp)+np.exp(-rp))/2.)
    N=A*(E/Epiv)**((alpha+beta)/2.)*10.**(a-ap)
    return N



def fBAND(E,Amplitude,Epeak, alpha, beta):
    Ec = (alpha-beta)*Epeak/(2.0+alpha)   #
    if E<=Ec:
        N=Amplitude*(E/100)**alpha * np.exp(-(alpha+2.)*E/Epeak)
    else:
        N=Amplitude*(E/100)**beta*np.exp(beta-alpha)*(((alpha-beta)*Epeak)/(100.0*(alpha+2.)))**(alpha-beta)
        
    return N


def fCPL(E,A,Epeak,index,Epiv):
    Epiv=100.
    N=A*(E/Epiv)**index*math.exp(-E*(2.+index)/Epeak)
    return N


def fPL(E,A,Epiv,index):
    #Epiv=100.
    N=A*(E/Epiv)**index
    return N

def fBB(E,A,kT):
    # np.exp(E/kT) is very large i. e. greater than 1.e304, then the denominator is infinite. Which means that the numerator is practically zero (~1.e-304), therefore do not  calculate N. Just put it equal to zero.
    #exp(E/kT) becomes greater than ~1.e304 when E/kT > 700.
    if  E/kT >700:
        N=0.
    else:
        N=A*(E**2/(np.exp(E/kT)-1.))
    return N


#================================================================================
#Flux functions
#================================================================================



def ENE_SBPL(E,A,Epiv,alpha,Ebreak,delta,beta):
    Epiv=100.
    delta=0.3
    r=math.log10(E/Ebreak)/delta
    rp=math.log10(Epiv/Ebreak)/delta
    a=1./2.*delta*(beta-alpha)*log((np.exp(r)+np.exp(-r))/2.)
    ap=1./2.*delta*(beta-alpha)*log((np.exp(rp)+np.exp(-rp))/2.)
    N=A*(E/Epiv)**((alpha+beta)/2.)*10**(a-ap)
    return E*N


def ENE_BAND(E,Amplitude,Epeak, alpha, beta):        
    if E<=(((alpha-beta)/(alpha+2))*Epeak):
        N=Amplitude*(E/100)**alpha * np.exp(-(alpha+2)*E/Epeak)
    else:
        N=Amplitude*(E/100)**beta*np.exp(beta-alpha)*(((alpha-beta)*Epeak)/(100*(alpha+2)))**(alpha-beta)        
    #E_iso = E*N
    return E*N
    #return E_iso*Tint*(4.*np.pi*D**2)/(1+z)*1.60217657e-9 # Tint is the time interval you used to perform the the joint GBM-LAT spectral analysis 


def ENE_CPL(E,A,Epeak,index,Epiv=100.):
    N=A*(E/Epiv)**index*math.exp(-E*(2+index)/Epeak)
    return E*N


def ENE_PL(E,A,Epiv,index):
    N=A*(E/Epiv)**index
    return E*N

def ENE_BB(E,A,kT):
    #N=A*(E**2/(np.exp(E/kT)-1.))
    if  E/kT >700:
        N=0.
    else:
        N=A*(E**2/(np.exp(E/kT)-1.))
    return E*N






#================================================================================
#Define integrals for various models. Integrals are calculated between the energy limits E1 and E2.
#================================================================================

def integral_SBPL(A,Epiv,alpha,Ebreak,delta,beta,E1,E2):
    Epiv=100.
    delta=0.3
    #print 'inside integral_SBPL z=',z
    #print A,Epiv,alpha,Ebreak,delta,beta
    return si.quad(ENE_SBPL, E1, E2,args=(A,Epiv,alpha,Ebreak,delta,beta))
    


def integral_BAND(Amplitude,Epeak,alpha,beta,E1,E2):
    #print 'inside integral_band z=',z
    #print Amplitude,Epeak,alpha,beta
    #test=si.quad(fBand, 1./(1.+z), 10000./(1.+z),args=(Amplitude,Epeak,alpha,beta))
    #print test
    return si.quad(ENE_BAND, E1, E2,args=(Amplitude,Epeak,alpha,beta))
    


def integral_CPL(A,Epeak,index,Epiv,E1,E2):
    Epiv=100.
    #print 'inside integral_cpl z=',z
    #print A,Epeak,index,Epiv
    return si.quad(ENE_CPL, E1, E2,args=(A,Epeak,index,Epiv))


def integral_PL(A,Epiv,index,E1,E2):
    Epiv=100.
    #print 'inside integral_PL z=',z
    #print A,Epiv,index
    return si.quad(ENE_PL, E1, E2,args=(A,Epiv,index))

def integral_BB(A,kT,E1,E2,FNsteps):
    # Steps need to be logarithmic
    Nsteps=int(FNsteps)
    E1log=math.log10(E1)
    E2log=math.log10(E2)
    
    vE=np.logspace(E1log,E2log,Nsteps)
    dE=np.zeros(Nsteps)
    for i in range(Nsteps-1):
        dE[i+1]=vE[i+1]-vE[i]
        
    integral=0.
    for i in range (Nsteps-1):
        Emid=vE[i]+dE[i+1]/2.
        ENE_mid=ENE_BB(Emid, A ,kT)
        integral=integral+ENE_mid*dE[i+1]

    return integral
    #return si.quad(ENE_BB, E1, E2,args=(A,kT))

def fPL_err2(E, A, Epiv, index, s2A, s2Epiv, s2index, sAindex):
    f = fPL(E, A, Epiv, index)
    Epiv=100.
    partialA_pl=f/A
    partialEpiv_pl=0.
    partialIndex_pl=f*np.log(E/Epiv)
    err2 = partialA_pl**2*s2A+partialIndex_pl**2*s2index+2.*partialA_pl*partialIndex_pl*sAindex
    return err2, partialA_pl, partialIndex_pl,partialEpiv_pl



def fBAND_err2(E, A_band, Epeak, alpha, beta, s2band_A, s2band_Epeak, s2band_alpha, s2band_beta, s_bandAEpeak, s_bandAalpha, s_bandAbeta, s_bandEpeakalpha, s_bandEpeakbeta, s_bandalphabeta):
    Ec = (alpha-beta)*Epeak/(2.0+alpha)   #
    f=fBAND(E, A_band, Epeak, alpha, beta)
    partialA_band = f / A_band
    if E <= Ec:
        partialEpeak_band = f * E * (2. + alpha) / Epeak ** 2
        partialalpha_band = f * (np.log(E / 100.) - E / Epeak)
        partialbeta_band = 0.
    elif E > Ec:
        partialEpeak_band = f * (alpha - beta) / Epeak
        partialalpha_band = f * (np.log(((alpha - beta) * Epeak) / (100. * (2. + alpha))) - (alpha - beta) / (2. + alpha))
        partialbeta_band = f * (np.log(E) - np.log(((alpha - beta) * Epeak) / (2. + alpha)))

    err2 = (partialA_band)**2*s2band_A+(partialEpeak_band)**2*s2band_Epeak+(partialalpha_band)**2*s2band_alpha+(partialbeta_band)**2*s2band_beta+2.*(partialA_band*partialEpeak_band*s_bandAEpeak+partialA_band*partialalpha_band*s_bandAalpha+partialA_band*partialbeta_band*s_bandAbeta+partialEpeak_band*partialalpha_band*s_bandEpeakalpha+partialEpeak_band*partialbeta_band*s_bandEpeakbeta+partialalpha_band*partialbeta_band*s_bandalphabeta)
    return err2, partialA_band, partialbeta_band, partialalpha_band, partialEpeak_band


def fCPL_err2(E, A, Epeak, index, Epiv, s2A, s2Epeak, s2index, s2Epiv,s_AEpeak, s_Aindex, s_Epeakindex):
    Epiv=100.
    f = fCPL(E, A, Epeak, index, Epiv)

    partialA_cpl = f / A
    partialEpeak_cpl = f * (2. + index) * E / Epeak ** 2
    partialIndex_cpl = f * (math.log(E / Epiv) - E / Epeak)
    partialEpiv_cpl=0.

    err2 = partialA_cpl ** 2 * s2A + partialEpeak_cpl ** 2 * s2Epeak + partialIndex_cpl ** 2 * s2index + 2.*(partialA_cpl*partialEpeak_cpl*s_AEpeak+partialA_cpl*partialIndex_cpl*s_Aindex+partialEpeak_cpl*partialIndex_cpl*s_Epeakindex)

    return err2, partialA_cpl, partialIndex_cpl, partialEpeak_cpl,partialEpiv_cpl


def fBB_err2(E, A, kT, s2A, s2kT, s_AkT):
    f = fBB(E, A, kT)
    # For very small f (less than 1.e304) do not calculate errors. The function is practically zero and so will be its errors. This avoids warnings and nan in the output.
    if f < 1.e304:
        partialA_bb = 0.
        partialkT_bb = 0.
    else:
        partialA_bb = f / A
        partialkT_bb = f * np.exp(E / kT) * E / kT ** 2 * (1. / (np.exp(E / kT) - 1.))

    err2 = (partialA_bb) ** 2 * s2A + (partialkT_bb) ** 2 * s2kT + 2.*partialA_bb*partialkT_bb*s_AkT

    return err2, partialA_bb, partialkT_bb


def fSBPL_err2(E, A, Epiv, alpha, Ebreak, delta, beta, s2A, s2Epiv, s2alpha, s2Ebreak, s2delta, s2beta, s_alphabeta,s_AEbreak,s_Aalpha, s_Abeta,s_Ebreakalpha, s_Ebreakbeta ):
    f = fSBPL(E, A, Epiv, alpha, Ebreak, delta, beta)
    r = math.log10(E / Ebreak) / delta
    rp = math.log10(Epiv / Ebreak) / delta
    a = 1. / 2. * delta * (beta - alpha) * np.log((np.exp(r) + np.exp(-r)) / 2.)
    ap = 1. / 2. * delta * (beta - alpha) * np.log((np.exp(rp) + np.exp(-rp)) / 2.)

    # partial derivative w. r. t. A
    partialA_sbpl = f / A
    partialEbreak_sbpl = f * (beta - alpha)/(2.*Ebreak)*((np.exp(rp)-np.exp(-rp))/(np.exp(rp)+np.exp(-rp))-((np.exp(r)-np.exp(-r))/(np.exp(r)+np.exp(-r))))
    partialalpha_sbpl = f / 2. * (np.log(E/100.)+delta*np.log(10)*(np.log(np.exp(rp)+np.exp(-rp))-np.log(np.exp(r)+np.exp(-r))))
    partialbeta_sbpl = f / 2. * (np.log(E/100.)+delta*np.log(10)*(np.log(np.exp(r)+np.exp(-r))-np.log(np.exp(rp)+np.exp(-rp))))
    partialEpiv_sbpl=0.
    partialdelta_sbpl=0.
    partialDelta_sbpl= partialdelta_sbpl

    err2 = partialA_sbpl ** 2 * s2A + partialEbreak_sbpl ** 2 * s2Ebreak + partialalpha_sbpl ** 2 * s2alpha + partialbeta_sbpl ** 2 * s2beta + 2.*(partialA_sbpl*partialEbreak_sbpl*s_AEbreak+partialA_sbpl*partialalpha_sbpl*s_Aalpha+partialA_sbpl*partialbeta_sbpl*s_Abeta+partialEbreak_sbpl*partialalpha_sbpl*s_Ebreakalpha+partialEbreak_sbpl*partialbeta_sbpl*s_Ebreakbeta+partialalpha_sbpl*partialbeta_sbpl*s_alphabeta)

    return err2, partialA_sbpl, partialbeta_sbpl, partialalpha_sbpl, partialEbreak_sbpl,partialEpiv_sbpl,partialDelta_sbpl




#=================================================================================
#========================= END FUNCTION DEFINITIONS  =============================
#=================================================================================













#=================================================================================
#================================== MAIN ========================================
#=================================================================================



#FILE NAMES
#==========

str_EAC=''
if eac:
    str_EAC='_EAC'
str_MODEL=model_complete+str_EAC
str_file=directory+str_MODEL+'.fit'  # 020_GRB130427324 BAND_PL_BB_EAC.fit file


str_modelname=model_complete.replace('_','+')
#open file
hdu=fits.open(str_file)   # open fits file

#MODEL, COMPONENTS
#=================
model_components=np.array(model_complete.split('_'))
basemodel=model_components[0]
print (model_components,basemodel)
# Get array of models, array of parameter numbers, and array of parameter names in rmfit order.
models,Nparameters,parnames,parfix=get_model_details(model_complete)
indexmain=np.where(models == basemodel)[0][0] #Find position of main model after reordering models according to rmfit order.
Nparameters_model_only=np.sum(Nparameters) #Total number of parameters in the model.
#Get number of model components (excluding EAC)
Ncomponents=models.shape[0]
if models[Ncomponents-1]=='EAC':
    Ncomponents=Ncomponents-1

#Find indices of fixed parameters in this model
indfix=np.where(parfix=='f')[0]
#Find indices of variable parameters in this model
indvary=np.where(parfix=='v')[0]

#Initialise parameter names
fullparnames=np.full(Nparameters_model_only,' '*40)
fullparnames_err=np.full(Nparameters_model_only,' '*40)


# GET COVARIANCE MATRIX AND MEAN VALUES OF ALL PARAMETERS
#========================================================
tmpcovariance_matrix=hdu[2].data.field('COVARMAT') # with EAC
if tmpcovariance_matrix[0][0].shape[0] != Nparameters_model_only:
    print('Incompatible covariance matrix size - correct for EAC')
covariance_matrix=tmpcovariance_matrix[0][0:Nparameters_model_only,0:Nparameters_model_only] # Copy covariance matrix for model parameters only (ignore EAC)

#Replace columns for fixed parameters with zero.
for j in indfix:
   covariance_matrix[:,j]=0.
#Replace rows for fixed parameters with zero.
for i in indfix:
    covariance_matrix[i,:]=0.                
    covariance_matrix[np.abs(covariance_matrix) < 1.e-18] = 0.


#Read parameters values and their erros    
param_mean=[]
param_err=[]

for ip in range(Nparameters_model_only):
    param_mean.append(hdu[2].data.field('PARAM'+str(ip))[0][0])
    param_err.append(hdu[2].data.field('PARAM'+str(ip))[0][1])
    # strout=parnames[ip].ljust(18)+'= '+str("%0.3e"%param_mean[ip]).rjust(12)+'    +/- '+str("%0.3e"%param_err[ip]).rjust(12)+'\n'
    # fileout.write(strout)    
    #print(parnames[ip],'= ',param_mean[ip],'+/-',param_err[ip])
param_mean=np.asarray(param_mean)

# Calculate flux, get Epeak array
#================================
# empty array for the flux, fluence for each of the simulated model.
nvect=1
flux_array=np.zeros(nvect)
fluence_array=np.zeros(nvect)
flux_BB_array=np.zeros(nvect)
flux_bolometric_array=np.zeros(nvect)
photon_flux_E=np.zeros(nvect) # photon flux at a given energy Energy_flux_fluence

Epeak_array=np.zeros(nvect)



vE=np.logspace(E1,E2,10000)   #Vector for energies keV

if plot_limit_vE_highest:
    ie_low=0
    ie_high= np.argwhere(vE >= Emax_highest)[0][0]-1
    
    #print(np.argwhere(vE >= Emax_highest))
    #print (ie_low, ie_high)
else:
    ie_low=0
    ie_high= vE.shape[0]
#print (ie_low, ie_high)
    
#exit()

# Initialisations for model
SED_components=[] # Create list for model component vectors for SED for later plots. This list will be saved in a npy file.
SED_components.append(vE)  #Add energy vector to list to make plots later
SED_sum=np.zeros(vE.shape)   #Create empy numpy array for SED for sum of models

isim=0
#Loop on models
ip=0   # number of parameters already read

for ic, component in enumerate(models):
    base=False
    if ic==indexmain:
        base=True
        
    if component=='PL':
        val_A=param_mean[ip]
        val_Epiv=param_mean[ip+1]
        val_Epiv=100.
        val_index=param_mean[ip+2]
        # Calculate flux
        # Call function/calculations
        vfPL = np.vectorize(fPL)
        vvfPL =vfPL(vE,val_A,val_Epiv, val_index)
        SED_PL=vE**2*vvfPL*keVtoErg
        SED_components.append(SED_PL) #Add PL SED to SED components list
        SED_sum=SED_sum + SED_PL
        plt.loglog(vE[ie_low:ie_high],SED_PL[ie_low:ie_high],component_linestyle, linewidth=2.,c=color10,alpha=0.8)

        val_s2A =covariance_matrix[ip,ip]
        val_s2Epiv=covariance_matrix[ip+1,ip+1]
        val_s2index=covariance_matrix[ip+3,ip+3]
        val_sAindex=covariance_matrix[ip, ip+3]

        vfPL_err2 = np.vectorize(fPL_err2)
        vvfPL_err2, partialA_pl, partialIndex_pl,partialEpiv_pl = vfPL_err2(vE, val_A, val_Epiv, val_index, val_s2A, val_s2Epiv, val_s2index, val_sAindex)

        SED_PL_ERR = vE ** 2 * np.sqrt(vvfPL_err2) * keVtoErg
        SED_PL_high = SED_PL + SED_PL_ERR
        SED_PL_low = SED_PL - SED_PL_ERR
        #plt.fill_between(vE, SED_PL_high, SED_PL_low, color=color10, alpha=.4)

    if component=='SBPL':
        val_A=param_mean[ip]
        val_Epiv=param_mean[ip+1]
        val_Epiv=100.
        val_alpha=param_mean[ip+2]
        val_Ebreak=param_mean[ip+3]
        val_delta=param_mean[ip+4]
        val_delta=0.3
        val_beta=param_mean[ip+5]
        #Epeak
        if base: #get rid of this?
            Epeak_array[isim]=val_Ebreak*10**(1./(2.)*val_delta*np.log((val_alpha+2.)/(-val_beta-2.)))
        vfSBPL = np.vectorize(fSBPL)
        vvfSBPL=vfSBPL(vE,val_A,val_Epiv,val_alpha, val_Ebreak, val_delta, val_beta )
        SED_SBPL=vE**2*vvfSBPL*keVtoErg
        SED_components.append(SED_SBPL) #Add PL SED to SED components list
        SED_sum=SED_sum + SED_SBPL
        # plot
        plt.loglog(vE[ie_low:ie_high],SED_SBPL[ie_low:ie_high],component_linestyle, linewidth=2.,c=color10,alpha=0.8)

        val_s2A =covariance_matrix[ip,ip]
        val_s2Epiv =covariance_matrix[ip+1,ip+1]
        val_s2alpha =covariance_matrix[ip+2,ip+2]
        val_s2Ebreak =covariance_matrix[ip+3,ip+3]
        val_s2delta=covariance_matrix[ip+4,ip+4]
        val_s2beta =covariance_matrix[ip+5,ip+5]
        val_s_alphabeta=covariance_matrix[ip+2,ip+5]
        val_s_AEbreak =covariance_matrix[ip,ip+3]
        val_s_Aalpha =covariance_matrix[ip,ip+2]
        val_s_Abeta =covariance_matrix[ip,ip+5]
        val_s_Ebreakalpha=covariance_matrix[ip+2,ip+3]
        val_s_Ebreakbeta=covariance_matrix[ip+3,ip+5]
        vfSBPL_err2 = np.vectorize(fSBPL_err2)
        vvfSBPL_err2, partialA_sbpl, partialbeta_sbpl, partialalpha_sbpl, partialEbreak_sbpl,partialEpiv_sbpl,partialDelta_sbpl = vfSBPL_err2(vE, val_A, val_Epiv, val_alpha, val_Ebreak, val_delta, val_beta, val_s2A, val_s2Epiv, val_s2alpha, val_s2Ebreak, val_s2delta,
                                   val_s2beta, val_s_alphabeta, val_s_AEbreak, val_s_Aalpha, val_s_Abeta, val_s_Ebreakalpha, val_s_Ebreakbeta)

        SED_SBPL_ERR = vE ** 2 * np.sqrt(vvfSBPL_err2) * keVtoErg
        SED_SBPL_high = SED_SBPL + SED_SBPL_ERR
        SED_SBPL_low = SED_SBPL - SED_SBPL_ERR
        #plt.fill_between(vE, SED_SBPL_high, SED_SBPL_low, color=color10, alpha=.4)

    if component=='BAND':
        val_A=param_mean[ip]
        val_Epeak=param_mean[ip+1]
        val_alpha=param_mean[ip+2]
        val_beta=param_mean[ip+3]
        
        #Epeak
        if base:
            Epeak_array[isim]=val_Epeak
        vfBAND = np.vectorize(fBAND)
        vvfBAND=vfBAND(vE,val_A,val_Epeak,val_alpha, val_beta )
        SED_BAND=vE**2*vvfBAND*keVtoErg
        SED_components.append(SED_BAND) #Add PL SED to SED components list
        SED_sum=SED_sum + SED_BAND
        # plot
        plt.loglog(vE[ie_low:ie_high],SED_BAND[ie_low:ie_high],component_linestyle, linewidth=2.,c=color10,alpha=0.8)

        val_s2band_A=covariance_matrix[ip,ip]
        val_s2band_Epeak=covariance_matrix[ip+1,ip+1]
        val_s2band_alpha=covariance_matrix[ip+2,ip+2]
        val_s2band_beta=covariance_matrix[ip+3,ip+3]
        val_s_bandAEpeak=covariance_matrix[ip,ip+1]
        val_s_bandAalpha=covariance_matrix[ip,ip+2]
        val_s_bandAbeta=covariance_matrix[ip,ip+3]
        val_s_bandEpeakalpha=covariance_matrix[ip+1,ip+2]
        val_s_bandEpeakbeta=covariance_matrix[ip+1,ip+3]
        val_s_bandalphabeta=covariance_matrix[ip+2,ip+3]

        print (val_Epeak, math.sqrt(val_s2band_Epeak))
        vfBAND_err2 = np.vectorize(fBAND_err2)
        vvfBAND_err2, partialA_band, partialbeta_band, partialalpha_band, partialEpeak_band = vfBAND_err2(vE, val_A, val_Epeak, val_alpha, val_beta, val_s2band_A, val_s2band_Epeak, val_s2band_alpha, val_s2band_beta,val_s_bandAEpeak, val_s_bandAalpha, val_s_bandAbeta, val_s_bandEpeakalpha, val_s_bandEpeakbeta, val_s_bandalphabeta)

        SED_BAND_ERR = vE ** 2. * np.sqrt(vvfBAND_err2) * keVtoErg
        SED_BAND_high = SED_BAND + SED_BAND_ERR
        SED_BAND_low = SED_BAND - SED_BAND_ERR
        #plt.fill_between(vE, SED_BAND_high, SED_BAND_low, color=color10, alpha=.4)

    if component=='CPL':
        val_A=param_mean[ip]
        val_Epeak=param_mean[ip+1]
        val_index=param_mean[ip+2]
        val_Epiv=param_mean[ip+3]

        val_Epiv=100.
        #Epeak
        if base:
            Epeak_array[isim]=val_Epeak
        vfCPL = np.vectorize(fCPL)
        #vvfCPL, partialA_cpl, partialIndex_cpl, partialEpeak_cpl =vfCPL(vE,val_A,val_Epeak,val_index, val_Epiv )
        vvfCPL =vfCPL(vE,val_A,val_Epeak,val_index, val_Epiv )
        SED_CPL=vE**2*vvfCPL*keVtoErg
        SED_components.append(SED_CPL) #Add PL SED to SED components list
        SED_sum=SED_sum + SED_CPL
        #plot
        plt.loglog(vE[ie_low:ie_high],SED_CPL[ie_low:ie_high],component_linestyle, linewidth=2.,c=color10,alpha=0.8)

        val_s2A=covariance_matrix[ip,ip]
        val_s2Epeak=covariance_matrix[ip+1,ip+1]
        val_s2index=covariance_matrix[ip+2,ip+2]
        val_s2Epiv=covariance_matrix[ip+3,ip+3]
        val_s_AEpeak=covariance_matrix[ip,ip+1]
        val_s_Aindex=covariance_matrix[ip,ip+2]
        val_s_Epeakindex=covariance_matrix[ip+1,ip+2]
        vfCPL_err2 = np.vectorize(fCPL_err2)
        vvfCPL_err2,partialA_cpl, partialIndex_cpl, partialEpeak_cpl,partialEpiv_cpl = vfCPL_err2(vE, val_A, val_Epeak, val_index, val_Epiv, val_s2A, val_s2Epeak, val_s2index, val_s2Epiv, val_s_AEpeak, val_s_Aindex,
                                 val_s_Epeakindex)

        SED_CPL_ERR = vE ** 2. * np.sqrt(vvfCPL_err2) * keVtoErg
        SED_CPL_high = SED_CPL + SED_CPL_ERR
        SED_CPL_low = SED_CPL - SED_CPL_ERR
        #plt.fill_between(vE, SED_CPL_high, SED_CPL_low, color=color10, alpha=.4)

    if component=='BB':
        val_A=param_mean[ip]
        val_kT=param_mean[ip+1]
        vfBB = np.vectorize(fBB)
        vvfBB=vfBB(vE,val_A,val_kT)
        SED_BB=vE**2*vvfBB*keVtoErg
        SED_components.append(SED_BB) #Add PL SED to SED components list
        SED_sum=SED_sum + SED_BB
        # plot
        plt.loglog(vE[ie_low:ie_high],SED_BB[ie_low:ie_high],component_linestyle, linewidth=2.,c=color10,alpha=0.8)

        val_s2A=covariance_matrix[ip,ip]
        val_s2kT=covariance_matrix[ip+1,ip+1]
        val_s_AkT=covariance_matrix[ip,ip+1]
        vfBB_err2 = np.vectorize(fBB_err2)
        vvfBB_err2,partialA_bb, partialkT_bb = vfBB_err2(vE, val_A, val_kT, val_s2A, val_s2kT, val_s_AkT)

        SED_BB_ERR = vE ** 2 * np.sqrt(vvfBB_err2) * keVtoErg
        SED_BB_high = SED_BB + SED_BB_ERR
        SED_BB_low = SED_BB - SED_BB_ERR
        #plt.fill_between(vE, SED_BB_high, SED_BB_low, color=color10, alpha=.4)

    ip=ip+Nparameters[ic]        

errs = 0.
for iparam1 in range (sum(Nparameters)):
    param1 = parnames[iparam1]
    for iparam2 in range(sum(Nparameters)):
        param2 = parnames[iparam2]
        sigmaparam = covariance_matrix[iparam1,iparam2]
        #print(param1, param2, sigmaparam)
        errs = errs+(eval('partial'+param1))*(eval('partial'+param2))*sigmaparam

totalerrs = vE**2*np.sqrt(errs)*keVtoErg
totalerrs_high = SED_sum+ totalerrs
totalerrs_low = SED_sum- totalerrs



if OPT_SHADE:
    plt.fill_between(vE[ie_low:ie_high], totalerrs_high[ie_low:ie_high], totalerrs_low[ie_low:ie_high], color=color10, alpha=.25)
    plt.loglog(vE[ie_low:ie_high],totalerrs_high[ie_low:ie_high],'-', linewidth=1.,c=color10,alpha=0.3)
    plt.loglog(vE[ie_low:ie_high],totalerrs_low[ie_low:ie_high],'-', linewidth=1.,c=color10,alpha=0.3)
else:
    plt.loglog(vE[ie_low:ie_high],totalerrs_high[ie_low:ie_high],'-', linewidth=1.,c=color10,alpha=0.7)
    plt.loglog(vE[ie_low:ie_high],totalerrs_low[ie_low:ie_high],'-', linewidth=1.,c=color10,alpha=0.7)
    print('')
plt.loglog(vE[ie_low:ie_high],SED_sum[ie_low:ie_high],'-', linewidth=3.,label=str_label,c=color10)



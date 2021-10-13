# ErrorCalculationTool
An error calculation tool for fits of gamma-ray burst (GRB) spectral energy distributions.

The application of empirical functions to the analysis of gamma-ray burst (GRB) spectra plays
an important role in understanding their fundamental nature. These functions, called spectral
models, are fitted to spectral energy distributions (SEDs). In order to obtain a true fit, two or
more models are often used in conjunction, and are referred to as the base model and additional
components.

The Gamma-ray Spectral Fitting Package (RMfit) calculates best-fit parameters for user-chosen
models and outputs a FITS (Flexible Image Transport System) file which contains information
about the quality of the fit. In order to calculate the propagation of errors in multiple models,
the following information is required from the FITS file:

(i) covariances between all the different parameters involved (i.e. covariance matrix)

(ii) values of the parameters

(iii) errors in the parameter values


Inputs:

The code calls the astropy.io.fits package in order to provide access to FITS files. The FITS
file type is one that is commonly used in astronomy and is capable of storing tables as well
as figures. More importantly, the FITS format allows the code to easily access the covariance
matrix, parameter values, and their errors.
FITS files can be obtained from RMfit, and should be named following the pattern:
BASEMODEL_ADDITIONALCOMPONENT1_ADDITIONALCOMPONENT2_EAC.fit. If the Effective Area
Correction (EAC) is not present in the particular fit, it should be omitted.

Opening a FITS file:

from astropy.io import fits
hdu=fits.open(str_file)
The names of spectral models being used in the analyses should be input in the list variable
models, with the first one being the base model. The portion of the FITS file name containing
the models is stored in the variable model_complete. If the FITS file contains an EAC term, the
eac variable should be set to 'True'. This adds '_EAC' to model_complete, which, along with
the directory path and file extension, form the file name (e.g. /path/BAND_PL_BB_EAC.fit).
The open() function returns an HDUList, which is a collection of HDU (Header Data Unit)
objects, usually consisting of a data array. For more information on HDULists and opening
FITS files, please see the astropy documentation.

Reading the HDU data:

Since we are interested only in obtaining the covariance matrix and mean values of the
parameters, we must only deal with the relevant part of the HDUList. Note that HDULists
can be indexed like arrays, and the covariance matrix lies in the second index. The following
command accesses the covariance HDU:
tmpcovariance_matrix=hdu[2].data.field('COVARMAT')
Similarly, we obtain the fitting parameters' errors and mean values using a for loop which
iterates over the total number of parameters:
for ip in range(Nparameters_model_only):
param_mean.append(hdu[2].data.field('PARAM'+str(ip))[0][0])
param_err.append(hdu[2].data.field('PARAM'+str(ip))[0][1])
where the zeroth element of the array gives the parameter's mean value and the rst element
gives its error.

Dealing with the models and their components:

The base model is separated from the additional components using the split() function. We
define a new function get_model_details which passes the names of all models and returns
them (models), reordered according to the RMfit convention. It also returns the number of
parameters for each model (Nparameters), their respective names (parnames), and whether
they are variable or fixed.

Calculating the SEDs:

We define functions which calculate the photon number 
ux (N) in s􀀀1cm􀀀2keV 􀀀1. The values
of each model's parameters are then passed to their respective functions, which return N.

Calculating errors:

The covariance matrix from the FITS file gives us the covariances between each of the different
parameters. The matrix's diagonal elements are the standard deviations of the parameters.

In order to calculate the errors for the sum of the spectral model components, we make use
of a nested for loop which pairs two parameters together and obtains the covariance between
them from the covariance matrix.
sigmaparam = covariance_matrix[iparam1,iparam2]
The correlated error of two parameters is calculated using the error propagation formula and
successively added to the variable errs.
errs = errs+(eval('partial'+param1))*(eval('partial'+param2))*sigmaparam
The error is then added and subtracted to the sum of the components, in order to get an
upper and lower bound.

Output:

The tool outputs a butterfly SED, displaying individual model components, as well as the best-
fit model. The transparent regions are representative of the errors calculated by the tool.

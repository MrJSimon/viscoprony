##############################################################################
##
## Author:      Jamie E. Simon
##
## Description: The script creates a shiftfactors, mastercurve and pronyseries
##              to desribe the viscoelastic behaviour of polymers and other
##              rate and temperature sensitive material systems.
##
##############################################################################


## Load in packages
import numpy as np
from PythonFunctions.Mastercurves.mastercurve import mastercurve_tanh
from PythonFunctions.PronySeries.prony_series import prony_1
from PythonFunctions.PlottingFunctions.plotting_functions import plot_storageloss_curves
from PythonFunctions.PlottingFunctions.plotting_functions import plot_mastercurve
from PythonFunctions.PlottingFunctions.plotting_functions import plot_pronyseries
from PythonFunctions.PlottingFunctions.plotting_functions import plot_shiftfactorsVStemperatures
from PythonFunctions.OutputFunctions.output_functions import numpysavetxt

## Load in values
values = np.loadtxt('Data/output_file_1.txt',delimiter = ',')

## Set X, Y1, Y2, Z
X  = values[:, 0] # Frequencies [rad/s] (converted from Hz later)
Ys = values[:, 1] # Storage modulus [MPa]
Yl = values[:, 2] # Loss modulus [MPa]
Z  = values[:, 3] # Temperature [°C] 

## Set temperature threshold values
T1, T2 = -30, 20  # Temperatures between -30°C and 20°C

## Set frequency threshold values
f1, f2 = 10**(-24), 10**24  # Frequencies between 10⁻²⁴ and 10²⁴ [1/s]

## Get indicides of values within thresshold
indicides_Z = np.where((Z>=T1) & (Z<=T2))

## Remove data-based on temperature index
X, Ys, Yl, Z = X[indicides_Z], Ys[indicides_Z], Yl[indicides_Z], Z[indicides_Z]

## Convert to radians pr. sec
X = 2.0*np.pi*X

## Plot storage and loss moduli curves
plot_storageloss_curves(np.column_stack((Z,X,Ys,Yl)),logscale = True, extra_stuff= True, picture_name = 'Output/StorageLossCurves.png')

## Set reference temperature
Tref = 20

## Get mastercurve solution
mastercurve_solution = mastercurve_tanh(X,Ys,Z,Tref,'tanh')

## Get X and Y output
Xout, a_T, temperatures = mastercurve_solution.Xout, mastercurve_solution.shift, mastercurve_solution.temperature_shift 

## Plot mastercurve
plot_mastercurve(np.column_stack((Z,X,Ys,Xout)),logscale = True, extra_stuff= True, picture_name = 'Output/MasterCurve.png')

## Plot shiftfactors vs. temperatures
plot_shiftfactorsVStemperatures(temperatures, np.log10(a_T), logscale = False, extra_stuff=True, picture_name = 'Output/ShiftfactorVsTemperatures.png')

## Get indicides
indicides_X = np.where((Xout>=f1) & (Xout<=f2))

## Remove data-based on frequency index
Xr, Ysr, Ylr, Zr = Xout[indicides_X], Ys[indicides_X], Yl[indicides_X], Z[indicides_X]

## Sort data based on incresing frequencies
index_rate = np.argsort(Xr)

## Sorted
Xr, Ysr, Ylr, Zr = Xr[index_rate], Ysr[index_rate], Ylr[index_rate], Zr[index_rate]

## Get prony series
PA = prony_1(Xr,Ysr,Ylr, w1 = 0.8, w2 = 0.2, maxiter=50)

## Get best fitting parameters
optimized_e, tau, optimized_esum, optimized_E0 = PA.result.x[:PA.nprony], PA.tau, PA.result.x[-2], PA.result.x[-1]

## Get optimized e based on the sum - The variable substitution method
optimized_e =PA.e_r(optimized_e,optimized_esum)

## Compute prediction
fitted_storage = PA.prony_storage(Xr, optimized_e, tau, optimized_E0)
fitted_loss = PA.prony_loss(Xr, optimized_e, tau, optimized_E0)

## Plot pronyseries
A, B = np.column_stack((Z,Xout,Ys)), np.column_stack((Xr,fitted_storage))
plot_pronyseries(A,B,logscale = True, extra_stuff= True, picture_name = 'Output/PronyseriesFit.png')

## Print the prony series
print('Prony series optimized values of e with a total sum = ' + str(np.sum(optimized_e)))
print(optimized_e[np.argsort(tau)])

print('Prony series optimized values of tau')
print(tau[np.argsort(tau)])

print('Prony series optimized values of E0')
print(optimized_E0)

## Create output format for numpy saving .txt - PRONY SERIES PARMATERS
A = np.column_stack((np.repeat(Tref,len(tau)),tau,optimized_e,np.repeat(optimized_E0,len(tau))))

## Create header for output
column_names = ['Tref [Celcius]','tau_i [s]', 'g_i [-]', 'E0 [MPa]']

## Save to output
numpysavetxt(filename='Output/PronySeriesCoefficientsDMTA.txt',data=A,headers=column_names,decimals=12,delimiter=',')

## Create output format for numpy saving .txt - temperatures and shiftfactors
A = np.column_stack((temperatures,a_T,np.repeat(Tref,len(a_T))))

## Create header for output
column_names = ['T [Celcius]','a_t [-]', 'Tref [Celcius]']

## Save to output
numpysavetxt(filename='Output/ShiftFactorsDMTA.txt',data=A,headers=column_names,decimals=12,delimiter=',')

## Create output format for numpy saving . txt - data and unshifted data
A = np.column_stack((Z,X,Xout,Ys,np.repeat(Tref,len(Z))))

## Create header for output
column_names = ['T [Celcius]','X [rad/sec]', 'Xs [rad/sec]', 'Y [MPa]', 'Tref [Celcius]']

## Save shifted data
numpysavetxt(filename='Output/ShiftedDataDMTA.txt',data=A,headers=column_names,decimals=12,delimiter=',')



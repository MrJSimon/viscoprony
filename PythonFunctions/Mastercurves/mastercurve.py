##############################################################################
##
## Author:      Theo Laurent, last modified by Jamie E. Simon
##
## Description: The script creates a master-curve using the tanh model
##              Ypred = a * np.tanh(b * np.log10(x) + c) + d
##
##############################################################################


## Load in modulues
import numpy as np
from scipy.optimize import least_squares

class mastercurve_tanh:
    """ This module fits a tanh model onto each temperature segment
        consisting of storage and a frequency data obtained from a
        DMTA machine.
        
    input
    ---------
    X: numpy array N x 1, frequency data measured [Hz]
    
    Y: numpy array N x 1, storage or loss modulus [MPa]
    
    Z: numpy array N x 1, temperature [celcius]
    
    output
    ---------
    class structure: contains shiftfactors and the fitted shifted data.
    
    """
    
    def __init__(self,Xi,Yi,Zi,Tref,dtype = 'linear'):
    
        ## Construct array types to ease computation
        Xarray,Yarray,pointer_t,pointer_f = self.construct_array_type(Xi,Yi,Zi)
        
        ## Build shiftfactor array
        shift = np.ones(len(pointer_t))
        
        ## Build tanh_model coefficient array
        Coeff = [[] for i in range(len(pointer_t))]
        
        if dtype == 'linear':
            ## Set initial guess
            p0 = [1,1]
            ## Set prediction statement
            Ypred = self.linear_model
            ## Compute model coefficients
            self.Coeff = self.get_linear_coeff(Xarray,Yarray,p0,Coeff,method='trf')
        
        if dtype == 'tanh':
            ## Set initial guess
            p0 = [1,1,1,1]
            ## Set prediction statement
            Ypred = self.tanh_model
            ## Compute model coefficients
            self.Coeff = self.get_tanh_coeff(Xarray,Yarray,p0,Coeff,method='trf')
        
        ## Get shiftfactor
        shift = self.get_shift_factors(Xarray,Yarray,Ypred,Coeff,shift)
        
        ## shift to desired reference temperature
        self.shift = self.shift_to_reference_temp(pointer_t,Tref,shift)
        
        ## Return fitted function and shifted x-values
        self.Xout,self.Yout,self.Xout_e,self.Yout_e = self.construct_shifted_vector(Xarray,Yarray,Ypred,Coeff,self.shift)
        
        ## Output temperatures associated with shift-factors
        self.temperature_shift = pointer_t
    
    def linear_model(self,x,p):
        """ 
        input
        ---------
        X: numpy array or float, frequency data measured [rad/s]
        
        p: list, fitting parameters
        
        output
        ---------
        tanh_model output
        """
        return (p[0]*np.log10(x)+p[1])
    
    
    def tanh_model(self,x,p):
        """ 
        input
        ---------
        X: numpy array or float, frequency data measured [rad/s]
        
        p: list, fitting parameters
        
        output
        ---------
        tanh_model output
        """
        return(p[0] * np.tanh(p[1] * np.log10(x) + p[2]) + p[3])
    
    def objfunc_1(self,p,Xi,Yi,Ypred):
        return(Yi- Ypred(Xi,p))
    
    def objfunc_2(self,x,coeff_n,Yr,Ypred):
        return(Yr - Ypred(x,coeff_n))
    
    def construct_array_type(self,Xi,Yi,Zi):
        """ 
        input
        ---------
        X: numpy array N x 1, frequency data measured [Hz]
        
        Y: numpy array N x 1, storage or loss modulus [MPa]
        
        Z: numpy array N x 1, temperature [celcius]
        
        output
        ---------
        Xarray: numpy array K x M, increasing frequency row-wize
        Yarray: numpy array K x M, increasing temperature and frequency along row and column respectively
        
        pointer_t: numpy array 
        
        """
        ## Get unique ID by identifying temperature
        pointer_t = np.unique(Zi)
        pointer_f = np.unique(Xi)
        ## Create new array format 
        Yarray = np.zeros((len(pointer_f),len(pointer_t)))
        Xarray = np.zeros((len(pointer_f),len(pointer_t)))
        ## Run through each element in unique pointer t and find associated values
        for i in range(0,len(pointer_t)):
            ## Set index
            index = np.isin(element = Zi,test_elements = pointer_t[i])
            ## Assign correct index to X- and Yarray
            Xarray[:,i] = Xi[index]
            Yarray[:,i] = Yi[index]
        return Xarray, Yarray, pointer_t, pointer_f

    def get_tanh_coeff(self,Xarray,Yarray,p0,Coeff,method='trf',apply_smoothing=False):
        """
        This function get the coefficients, used in the tanh trend curve.
        
        Input
        ----------
        Xarray: numpy array N x K, frequency data
        
        Yarray: numpy array N x K, storage or loss modulus
        
        p0: list, containing initial guess [a0,b0,c0,d0]

        Output
        -------
        Coeff : list, containing the optimized value for each temperature curve.

        """
        ## Run through all temperature curves
        for i in range(Yarray.shape[1]):
            ## Set prediction function
            Ypred = self.tanh_model
            ## Set current working X and Y values
            Xi, Yi = Xarray[:,i], Yarray[:,i]
            ## Apply smoothing if data points are too few (threshold can be adjusted)
            if apply_smoothing and len(Yi) < 5:
                Xi, Yi = self.smooth_data(Xi, Yi, window_size=3)
            ## Compute least squares results
            result = least_squares(self.objfunc_1,p0,args=(Xi,Yi,Ypred),method=method)
            ## Get coefficients
            a,b,c,d = result.x
            Coeff[i].append(a)
            Coeff[i].append(b)
            Coeff[i].append(c)
            Coeff[i].append(d)       
        return Coeff
        
    def get_linear_coeff(self,Xarray,Yarray,p0,Coeff,method='trf'):
        """
        This function get the coefficients, used in the linear trend curve.
        
        Input
        ----------
        Xarray: numpy array N x K, frequency data
        
        Yarray: numpy array N x K, storage or loss modulus
        
        p0: list, containing initial guess [a0,b0,c0,d0]

        Output
        -------
        Coeff : list, containing the optimized value for each temperature curve.

        """
        ## Run through all temperature curves
        for i in range(Yarray.shape[1]):
            ## Set prediction function
            Ypred = self.linear_model
            ## Set current working X and Y values
            Xi, Yi = Xarray[:,i], Yarray[:,i]
            ## Compute least squares results
            result = least_squares(self.objfunc_1,p0,args=(Xi,Yi,Ypred),method=method)
            ## Get coefficients
            a,b = result.x
            Coeff[i].append(a)
            Coeff[i].append(b)      
        return Coeff
    
    def identify_overlapping_region(self,Xi1,Xi2,Yi1,Yi2,tolerance=0.01):
        """
        This function identify the overlapping region between two curves, if any.
        
        Input
        ----------
        Xi1, Yi1 : Arrays for the first curve (x and y values)
        
        Xi2, Yi2 : Arrays for the second curve (x and y values)
        
        tolerance : float, optional
            Maximum allowable difference between y-values for overlap (default is 0.01)

        Output
        -------
        Xr, Xc : Arrays representing the x-values of the overlapping region

        """
        ## Identify the miminum value between Yi1 and Yi2
        ## Step 1: Identify the highest maximum y-value between both curves
        ymax = min(np.max(Yi1), np.max(Yi2))
        
        ## Step 2: Identify the lowest minimum y-value between both curves
        ymin = max(np.min(Yi1), np.min(Yi2))
        
        ## Step 3: Find common x-values where y-values fall within [ymin, ymax]
        Xr = Xi1[np.logical_and(Yi1 >= ymin, Yi1 <= ymax)]
        Xc = Xi2[np.logical_and(Yi2 >= ymin, Yi2 <= ymax)]
        
        ## Step 4: Handle cases where no overlapping region is found
        if len(Xr) == 0 or len(Xc) == 0:
            print('No overlapping region found')
            # Find the x-values corresponding to ymin and ymax for Xi1 and Xi2
            Xr_ymin = Xi1[np.argmin(np.abs(Yi1 - ymin))]
            Xr_ymax = Xi1[np.argmin(np.abs(Yi1 - ymax))]
            
            Xc_ymin = Xi2[np.argmin(np.abs(Yi2 - ymin))]
            Xc_ymax = Xi2[np.argmin(np.abs(Yi2 - ymax))]
            
            # Return the boundary x-values based on ymin and ymax
            return np.array([Xr_ymin, Xr_ymax]), np.array([Xc_ymin, Xc_ymax])
        print('Overlapping region identified')
        return Xr, Xc
        
    
    def get_shift_factors(self,Xarray,Yarray,Ypred,Coeff,shift):
        """
        This function computes the shiftfactors to generate the mastercurve
        
        Input
        ----------
        Xarray: numpy array N x K, frequency data
        
        Yarray: numpy array N x K, storage or loss modulus
        
        Ypred: function, describing the trend-curve.
        
        Coeff: list, containing the best fit parameters for each temperature curve.

        Output
        -------
        shift : numpy array N x 1, containing the shiftfactors.

        """
        ## Run through all temperatures
        for i in range(0,Yarray.shape[1]-1):
            ## Compute x values between two succesive temperature increments
            Xr = Xarray[:,i]
            Xc = Xarray[:,i+1]
            Xr = np.logspace(np.log10(np.min(Xr)),np.log10(np.max(Xr)),endpoint=True,num=int(len(Xr)*5))
            Xc = np.logspace(np.log10(np.min(Xc)),np.log10(np.max(Xc)),endpoint=True,num=int(len(Xc)*5))
            ## Compute prediction statement between two succesive temperature increments
            Yr = Ypred(Xr,Coeff[i])
            Yc = Ypred(Xc,Coeff[i+1])
            ## Identify overlapping region between curves if any
            Xr, Xc = self.identify_overlapping_region(Xr,Xc,Yr,Yc)
            ## Compute reference value
            Yr = Ypred(Xr,Coeff[i])
            Yc = Ypred(Xc,Coeff[i+1])
            ## Compute least squares results
            result = least_squares(self.objfunc_2,np.mean(Xc),args=(Coeff[i+1],Yr,Ypred),method='trf')
            ## Set resulting/optimum shift values
            X_obj = result.x
            ## append to output shift list
            shift[i+1:] = shift[i+1:] * np.mean(Xr)/X_obj
        return shift
    
    def shift_to_reference_temp(self,pointer_t,Tref,shift):
        """
        This function shift to the desired reference temperature.
        
        Input
        ----------
        pointer_t: numpy array N x 1, containing unique temperatures.
        
        Tref: float, reference temperature.
        
        shift: numpy array N x 1, containing shiftfactors
        
        Output
        -------
        shift: numpy array N x 1 , containing the shiftfactors, in relation to Tref

        """
        return shift/shift[list(pointer_t).index(Tref)]
    
    def construct_shifted_vector(self,Xarray,Yarray,Ypred,Coeff,shift):
        """
        This function computes the shifted shifted frequency values and 
        prediction statement using the trend-curves.
        
        Input
        ----------
        Xarray: numpy array N x K, frequency data
        
        Yarray: numpy array N x K, storage or loss modulus
        
        Ypred: function, describing the trend-curve.
        
        Coeff: list, containing the best fit parameters for each temperature curve.
        
        shift : numpy array N x 1, containing the shiftfactors.

        Output
        -------
        Xout: numpy array J x 1, shifted frequency data
        
        Yout: numpy array J x 1, trendline data
        
        Xenhanced: numpy array L x 1, ....
        
        Yenhanced: numpy array L x 1, ....
            
        """
        
        ## Initiate output list
        Xout,Yout,Yenhanced,Xenhanced = [], [], [], []
        
        ## Create an output based on the fit using original X data
        for i in range(0,Yarray.shape[1]):
            for val in Xarray[:,i]:
                Xout.append(val*shift[i])
                Yout.append(Ypred(val,Coeff[i]))
        
        ## Create an output based on the fit using original X enhanced range Xrr
        for i in range(0,Yarray.shape[1]):
            ## Create enhanced xrange
            Xrr = np.linspace(np.min(Xarray[:,i]),np.max(Xarray[:,i]),endpoint=True,num=50)
            for val in Xrr:
                Xenhanced.append(val*shift[i])
                Yenhanced.append(Ypred(val,Coeff[i]))    
        
        ## Covert to numpy array and return
        return np.array(Xout), np.array(Yout), np.array(Xenhanced),np.array(Yenhanced)
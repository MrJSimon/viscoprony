##############################################################################
##
## Author:      Mahdi Tayabeti and Jamie E. Simon
##
## Description: The script creates a prony series
##
##############################################################################


## Load in modulues
import numpy as np
#from scipy.optimize import minimize
from scipy.optimize import differential_evolution

class prony_1:
    """ This module creates a prony series
        
    input
    ---------
    X: numpy array N x 1, frequency data measured [rad/sec]
    
    Y1: numpy array N x 1, storage modulus [MPa]
    
    Y2: numpy array N x 1, storage loss modulus [MPa]
    
    w1: float, weight value, govering storage has to be between 0,1
    
    w2: float, weight value, govering loss has to be between 0,1
    
    maxiter: int, number of iteration to the differential_evo algorithm
    
    output
    ---------
    class structure: contains prony parameters
    
    """

    def __init__(self,Xi,Yi1,Yi2,w1=1.0,w2=1.0,maxiter=1000):
        
        # Define material properties: short-term modulus (E0)
        E_storage_short, E_loss_short = np.max(Yi1), np.max(Yi2)
       
        ## Compute logarithmic decades
        self.log_decades = self.get_log_decades(Xi)
        
        ## Initiate porny parameters
        initial_params,self.nprony,self.tau = self.init_prony_parameters(Xi,self.log_decades)
        
        ## Set bounds for e_i between 0 and 1, and for E_0 to be positive
        bounds = [(1/np.max(Xi), 1/np.min(Xi))] * 1 * self.nprony + [(0, 0.9998)]  + [(0.5*np.max(Yi1), 1.5*np.max(Yi1))]
        
        ## Utilize differential evolutionary algorithm 
        self.result = differential_evolution(self.objective_function,
        bounds,
        args=(Xi, Yi1, Yi2, E_storage_short, E_loss_short, self.tau, w1, w2),
        strategy='best1bin',
        maxiter=maxiter,
        popsize=15,
        tol=1e-7,
        mutation=(0.5, 1),
        recombination=0.7,
        disp = True)
               
    def get_log_decades(self,Xi):
        """ This module computes the number of logorithmic decades.
            
        input
        ---------
        Xi: numpy array N x 1, frequency data
          
        output
        ---------
        log_decades: float, number of logarithmic decades.
        
        """
        ## Get number of logorithmic decades
        log_decades = np.round(np.log10(np.max(Xi)/np.min(Xi)))
        return log_decades

    def init_prony_parameters(self,Xi,log_decades):  
        """ This module computes the initial parameters used in the
            optimization procedure.
            
        input
        ---------
        Xi: numpy array N x 1, frequency data
        
        log_decades: float, number of logarithmic decades.
          
        output
        ---------
        Inital_params: numpy array, initial guess [e^h_i0,...,e^h_Nprony0,esum,E0]
        
        nprony: int: number of prony series, based on logorithmic decades
        
        tau: numpy array, nprony x 1, relaxation timies.
        
        """
        ## Set number of logorithmic decades        
        nprony = int(1.0*log_decades)
        
        ## Assume tau can be predetermined based on omega input
        tau_min = np.log10(1.0/(np.max(Xi)/2.0*np.pi))
        tau_max = np.log10(1.0/(np.min(Xi)/2.0*np.pi))
        
        ## Compute constant relaxation parameters
        tau = np.logspace(tau_min,tau_max,nprony,endpoint=True)
        
        # Initial guess for e_i values (same length as tau) and E_0
        e, ehat, E0, esum = np.ones(nprony) * 0.5, np.ones(nprony) * 0.5, 1.0, 0.1
 
        ## Combine into a single array (to accomondate minimimzation syntax)
        initial_params = np.concatenate((ehat,np.array([esum,E0])))
        
        return initial_params,nprony,tau

    # Substitution parameter
    def e_r(self,ehat,e_sum):
        """ This module ensures that the sum of e_i is below 1, by introducing
            a numerical substitution parameter.
            
        input
        ---------
        ehat: numpy array, nprony x 1 - optimization parameter
        
        e_sum: float, optimization parameter with restriction e_sum <= 1.0
          
        output
        ---------
        e_i: numpy array, nprony x 1
        
        """
        e_i = ehat*(e_sum/np.sum(ehat))
        return e_i
    
    # Prony series functions for Storage (E') and Loss (E'') moduli
    def prony_storage(self,omega, e, tau, E_0):
        """ This module computes the storage modulus using a prony series function
        
        input
        ---------
        omega: numpy arra, N x 1, radial frequency
        
        e: numpy array nprony x 1, normalized e = E_i/E_0
        
        tau: numpy array nprony x 1, relaxation times
        
        E0: float, Short term modulus
          
        output
        ---------
        numpy array N x 1 prony, prony series storage.
        
        """
        return E_0 * (1 - np.sum(e)) + E_0 * np.sum([e_i * tau_i**2 * omega**2 / (1 + tau_i**2 * omega**2) for e_i, tau_i in zip(e, tau)], axis=0)

    def prony_loss(self,omega, e, tau, E_0):
        """ This module computes the loss modulus using a prony series function
        
        input
        ---------
        omega: numpy arra, N x 1, radial frequency
        
        e: numpy array nprony x 1, normalized e = E_i/E_0
        
        tau: numpy array nprony x 1, relaxation times
        
        E0: float, Short term modulus
          
        output
        ---------
        numpy array N x 1 prony, prony series loss.
        
        """
        return E_0 * np.sum([e_i * tau_i * omega / (1 + tau_i**2 * omega**2) for e_i, tau_i in zip(e, tau)], axis=0)

    # Objective function λ^2 to minimize with regularization
    def objective_function(self,params, omega, E_storage_data, E_loss_data, E_storage_short, E_loss_short, tau, w_s, w_l):
        """ This module computes the objective functionm used to minimize
            the prony-prediction to the data.
            
        input
        ---------
        params: numpy array, [e^h_i,...,e^h_Nprony,esum,E0]
        
        omega: numpy array N x 1, radial frequency
        
        E_storage_data: numpy array N x 1, storage modulus [MPa]
        
        E_loss_data: numpy array N x 1, storage loss modulus [MPa]
        
        E_storage_short: float, maximum value in E_storage_data [MPa]
        
        E_loss_short: float, maximum value in E_loss_data [MPa]
        
        tau: Numpy array nprony x 1, relaxation times [s]
        
        w1: float, weight value, govering storage has to be between 0,1
        
        w2: float, weight value, govering loss has to be between 0,1
        
        output
        ---------
        O_total: float, residual/objective
        
        """
        ## Set number of prony series
        nprony = self.nprony
        
        ## Initiate fitting parameters
        ehat,esum,E_0 = params[:nprony],params[-2],params[-1]
        ehat,esum,E_0 = params[:nprony],params[-2],params[-1]
        
        ## Empose variable substitution
        esubs = self.e_r(ehat,esum)
        
        ## Get the prony fit associated with the storage and loss modulus
        E_s = self.prony_storage(omega, esubs, tau, E_0)
        E_l = self.prony_loss(omega, esubs, tau, E_0)
        
        ## Compute prediction for each term
        O_s = (1.0/E_storage_short**2)*np.sum((E_storage_data - E_s)**2)
        O_l = (1.0/E_loss_short**2)*np.sum((E_loss_data - E_l)**2)
        
        ## Compute total prediction
        O_total = O_s*w_s + O_l*w_l
        
        return O_total
import numpy as np
from syrtis.tools import *

class Material:
    """
    Placeholder
    """

class Solid(Material):
    """
    Object for a homogenous solid Shell

    Args:
        name (str):         Name of the material
        k (float):          Thermal conductivity (W/m/K)
        rho (float):        Density (kg/m3)
        cp (float):         Specific heat capacity (kJ/kg/K)
        absorb (float):     Absorptivity, between 0 and 1
        emit (float):       Emissivity, between 0 and 1. Optional, defaults equal to absorbivity
        transmit (float):   Transmissivity, between 0 and 1, Optional, defaults to zero
    """
    def __init__(self, name, k, rho, cp, absorb, emit=0,  transmit=0):
        assert is_numeric(k, positive=True),   "Material 'k' must be a positive numerical value"
        assert is_numeric(rho, positive=True),   "Material 'rho' must be a positive numerical value"
        assert is_numeric(cp, positive=True), "Material 'cp' must be a positive numerical value"
        assert is_numeric(absorb, unit=True), "Material 'absorb' must be between 0 and 1"
        assert is_numeric(emit, unit=True), "Material 'emit' must be between 0 and 1"
        assert is_numeric(transmit, unit=True), "Material 'alpha' must be between 0 and 1"

        self.name = name
        self.k = k
        self.rho = rho
        self.cp = cp
        self.absorb = absorb

        if emit == 0:
            self.emit = self.absorb
        else:
            self.emit = emit
            
        self.transmit = transmit

class ConstrainedIdealGas(Material):
    """
    Object for a constrained ideal gas with no imposed flow (undergoing natural convection) and at constant pressure
    All temperature-dependent properties (density, beta) are calculated in functions of the same name
    Reference specific heat capacity cp_ref is given at standard temperature and pressure (101.3kPa, 298K)

    Args:
        p (float):      Pressure (Pa)
        M (float):      Molecular mass number
        Pr (float):     Prandtl number
        mu (float):     Absolute viscosity (Pa s)
        cp_ref (float): Reference isobaric specific heat capacity (kJ/kg/K)
        k (float):      Thermal conductivity (W/m/K)
        p (float):      Pressure (Pa) - optional, used as constant pressure if filled
        T (float):      Temperature (K) - optional, used as constant temperature if filled
        beta (float):   Thermal expansion coefficient (1/K) - optional, set to 1/T if left blank
    """
    def __init__(self, name, p, M, Pr, mu, cp, k, T=0, beta=0):
        assert is_numeric(p, positive=True),     "Material 'p' must be a positive numerical value"
        assert is_numeric(M, positive=True),     "Material 'M' must be a positive numerical value"
        assert is_numeric(Pr, positive=True),   "Material 'Pr' must be a positive numerical value"
        assert is_numeric(mu, positive=True),   "Material 'mu' must be a positive numerical value"
        assert is_numeric(cp, positive=True),   "Material 'cp' must be a positive numerical value"
        assert is_numeric(k, positive=True),     "Material 'k' must be a positive numerical value"

        if beta != 0:
            assert is_numeric(beta, positive=True),   "Material 'beta' must be a positive numerical value"

        if T != 0:
            assert is_numeric(T, positive=True),   "Material 'T' must be a positive numerical value"

        self.name = name
        self._p = p
        self._M = M
        self._Pr = Pr
        self._mu = mu
        self._cp = cp
        self._k = k

        if beta == 0:
            self.input_beta = True
            self._beta = 0.0
        else:
            self.input_beta = False
            self._beta = beta
        
        if T == 0:
            self.input_T = True
            self._T = 0.0
        else:
            self.input_T = False
            self._T = T
    
    def rho(self, T=0):
        """
        Gas density. If T is held constant, the input values will be ignored

        Args:
            T (float):  Temperature (K) 
        """
        
        if self.input_T == False:
            T = self._T
        else:
            assert is_numeric(T, positive=True),   "Input 'T' must be a positive numerical value"
        
        rho = self._p / ((8314 / self._M) * T) 

        return(rho)
    
    def beta(self, T):
        """
        Coefficient of thermal expansion, equal to 1/T. If T is held constant, the input value is ignored
        
        Args:
            T (float):  Temperature (K) 
        """

        if self.input_T == False:
            T = self._T
        else:
            assert is_numeric(T, positive=True),   "Input 'T' must be a positive numerical value"

        beta = 1/T

        return(beta)

    def mu(self, T):
        """
        Dynamic viscosity

        Args:
            T (float):  Temperature (K) 
        """

        if self.input_T == False:
            T = self._T
        else:
            assert is_numeric(T, positive=True),   "Input 'T' must be a positive numerical value"

        return(self._mu)

    def cp(self, T):
        """
        Specific heat capacity at constant pressure

        Args:
            T (float):  Temperature (K) 
        """

        if self.input_T == False:
            T = self._T
        else:
            assert is_numeric(T, positive=True),   "Input 'T' must be a positive numerical value"

        return(self._cp)

    def Pr(self, T):
        """
        Prandtl number

        Args:
            T (float):  Temperature (K) 
        """

        if self.input_T == False:
            T = self._T
        else:
            assert is_numeric(T, positive=True),   "Input 'T' must be a positive numerical value"

        return(self._Pr)
    
    def k(self, T):
        """
        Thermal conductivity

        Args:
            T (float):  Temperature (K) 
        """

        if self.input_T == False:
            T = self._T
        else:
            assert is_numeric(T, positive=True),   "Input 'T' must be a positive numerical value"

        return(self._k)
    
    def Re(self, T, L, v):
        """
        Dynamic viscosity

        Args:
            T (float):  Temperature (K) 
        """

        if self.input_T == False:
            T = self._T
        else:
            assert is_numeric(T, positive=True),   "Input 'T' must be a positive numerical value"
        
        assert is_numeric(L, positive=True),   "Input 'L' must be a positive numerical value"
        assert is_numeric(v, positive=True),   "Input 'v' must be a positive numerical value"

        Re = self.rho(T) * v * L / self.mu(T)

        return(Re)

    def Ra(self, T_avg, T_delta, g, length):
        """
        Calculate the Rayleigh number, using constant fluid properties

        
        """
        
        if self.input_T == False:
            T_avg = self._T
        else:
            assert is_numeric(T_avg, positive=True),    "Input 'T_avg' must be a positive numerical value"

        assert is_numeric(T_delta),      "Input 'T_delta' must be a numerial value"

        assert is_numeric(g, positive=True),            "Input 'g' must be a positive numerical value"
        assert is_numeric(length, positive=True),       "Input 'length' must be a positive numerical value"

        Ra = ((g * np.power(length, 3) * self.beta(T_avg) * abs(T_delta) * np.power(self.rho(T_avg), 2)) 
        / np.power(self.mu(T_avg), 2))

        assert Ra > 0, "Input has produced an invalid output"

        return(Ra)




material_classes = [Solid, ConstrainedIdealGas]

ambient_atmosphere = ConstrainedIdealGas("Martian ambient CO2", 580, 44, 0.71, 10.9e-6, 749, 0.0153, T=210)

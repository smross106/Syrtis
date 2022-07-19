

from numbers import Number
import numpy as np

class Material:
    """
    Placeholder
    """

class Solid(Material):
    """
    Object for a homogenous solid Shell

    Args:
        k (float):      Thermal conductivity (W/m/K)
        rho (float):    Density (kg/m3)
        cp (float):     Specific heat capacity (kJ/kg/K)

    """
    def __init__(self, k, rho, cp):
        assert isinstance(k, Number) and k > 0,   "Material 'k' must be a positive numerical value"
        assert isinstance(rho, Number) and rho > 0,   "Material 'rho' must be a positive numerical value"
        assert isinstance(cp, Number) and cp > 0, "Material 'cp' must be a positive numerical value"

        self.k = k
        self.rho = rho
        self.cp = cp

class ConstrainedIdealGas(Material):
    """
    Object for a constrained ideal gas with no imposed flow (undergoing natural convection)
    All temperature-dependent properties (density, beta) are calculated in functions of the same name
    Reference specific heat capacity cp_ref is given at standard temperature and pressure (101.3kPa, 298K)

    Args:
        M (float):      Molecular mass number
        Pr (float):     Prandtl number
        mu (float):     Absolute viscosity (Pa s)
        cp_ref (float): Reference isobaric specific heat capacity (kJ/kg/K)
        k (float):      Thermal conductivity (W/m/K)
        p (float):      Pressure (Pa) - optional, used as constant pressure if filled
        T (float):      Temperature (K) - optional, used as constant temperature if filled
        beta (float):   Thermal expansion coefficient (1/K) - optional, set to 1/T if left blank
    """
    def __init__(self, M, Pr, mu, cp, k, T=0, p=0, beta=0):
        assert isinstance(M, Number) and M > 0,     "Material 'M' must be a positive numerical value"
        assert isinstance(Pr, Number) and Pr > 0,   "Material 'Pr' must be a positive numerical value"
        assert isinstance(mu, Number) and mu > 0,   "Material 'mu' must be a positive numerical value"
        assert isinstance(cp, Number) and cp > 0,   "Material 'cp' must be a positive numerical value"
        assert isinstance(k, Number) and k > 0,     "Material 'k' must be a positive numerical value"

        if beta != 0:
            assert isinstance(beta, Number) and beta > 0,   "Material 'beta' must be a positive numerical value"

        if T != 0:
            assert isinstance(T, Number) and T > 0,   "Material 'T' must be a positive numerical value"

        if p != 0:
            assert isinstance(p, Number) and p > 0,   "Material 'p' must be a positive numerical value"

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
        
        if p == 0:
            self.input_p = True
            self._p = 0.0
        else:
            self.input_p = False
            self._p = p
    
    def rho(self, T=0, p=0):
        """
        Gas density. If T and p are held constant, the input values will be ignored

        Args:
            T (float):  Temperature (K) 
            p (float):  Pressure (Pa)
        """
        
        if self.input_p == False:
            p = self._p
        else:
            assert isinstance(p, Number) and p > 0,   "Input 'p' must be a positive numerical value"
        
        if self.input_T == False:
            T = self._T
        else:
            assert isinstance(T, Number) and T > 0,   "Input 'T' must be a positive numerical value"
        
        rho = p / ((8314 / self._M) * T) 

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
            assert isinstance(T, Number) and T > 0,   "Input 'T' must be a positive numerical value"

        beta = 1/T

        return(beta)

    def Ra(self, T_avg, T_delta, p, g, length):
        """
        Calculate the Rayleigh number, using constant fluid properties

        
        """
        
        if self.input_T == False:
            T_avg = self._T
        else:
            assert isinstance(T_avg, Number) and T_avg > 0,   "Input 'T_avg' must be a positive numerical value"

        assert isinstance(T_delta, Number), "Input 'T_delta' must be a numerial value"

        if self.input_p == False:
            p = self._p
        else:
            assert isinstance(p, Number) and p > 0,   "Input 'p' must be a positive numerical value"

        assert isinstance(g, Number) and g > 0,   "Input 'g' must be a positive numerical value"
        assert isinstance(length, Number) and length > 0,   "Material 'T' must be a positive numerical value"

        Ra = ((g * np.power(length, 3) * self.beta(T_avg) * abs(T_delta) * np.power(self.rho(T_avg, p), 2)) 
        / np.power(self._mu, 2))

        assert Ra > 0, "Input has produced an invalid output"

        return(Ra)




material_classes = [Solid, ConstrainedIdealGas]

ambient_atmosphere = ConstrainedIdealGas(44, 0.71, 10.9e-6, 749, 0.0153, T=210, p=580)

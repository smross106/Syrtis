

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

class ConstrainedGas(Material):
    """
    Object for a constrained gas with no imposed postion (undergoing natural convection)

    Args:
        T (float):      Temperature (K)
        p (float):      Pressure (Pa)
        Pr (float):     Prandtl number
        mu (float):     Absolute viscosity (Pa s)
        cp (float):     Isobaric specific heat capacity (kJ/kg/K)
        k (float):      Thermal conductivity (W/m/K)
        rho (float):    Density (kg/m3) 
        beta (float):   Thermal expansion coefficient (1/K) - optional, set to 1/T if left blank

    """
    def __init__(self, T, p, Pr, mu, cp, k, rho, beta=0.0):
        assert isinstance(T, Number) and T > 0,   "Material 'T' must be a positive numerical value"
        assert isinstance(p, Number) and p > 0,   "Material 'p' must be a positive numerical value"
        assert isinstance(Pr, Number) and Pr > 0, "Material 'Pr' must be a positive numerical value"
        assert isinstance(mu, Number) and mu > 0, "Material 'mu' must be a positive numerical value"
        assert isinstance(cp, Number) and cp > 0, "Material 'cp' must be a positive numerical value"
        assert isinstance(k, Number) and k > 0,   "Material 'k' must be a positive numerical value"
        assert isinstance(rho, Number) and rho > 0,   "Material 'rho' must be a positive numerical value"

        if beta != 0.0:
            assert isinstance(beta, Number) and beta > 0,   "Material 'beta' must be a positive numerical value"

        self.T = T
        self.p = p
        self.Pr = Pr
        self.mu = mu
        self.cp = cp
        self.k = k
        self.rho = rho

        if beta == 0.0:
            self.beta = 1 / self.T
        else:
            self.beta = beta


    def Ra(self, T_delta, g, length):
        """
        Calculate the Rayleigh number, using constant fluid properties
        """
        
        assert isinstance(T_delta, Number), "Input 'T_delta' must be a numerial value"
        assert isinstance(g, Number) and g > 0,   "Input 'g' must be a positive numerical value"
        assert isinstance(length, Number) and length > 0,   "Material 'T' must be a positive numerical value"

        Ra = ((g * np.power(length, 3) * self.beta * abs(T_delta) * np.power(self.rho, 2)) 
        / np.power(self.mu, 2))

        assert Ra > 0, "Input has produced an invalid output"

        return(Ra)




material_classes = [Solid, ConstrainedGas]

ambient_atmosphere = ConstrainedGas(210, 580, 0.71, 10.9e-6, 749, 8.74e-3, 0.0143)
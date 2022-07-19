"""
Classes and inherited subclasses for Shells, the general term for layers of material included in the simulation.
In most cases, a Shell will be added as a subclass with more specialised properties

References:

 - [1] - https://en.wikipedia.org/wiki/Thermal_resistance
 - [2] - Y Cengel, Heat Transfer

"""


from material import *



class Shell:
    """
    Object to store a single Shell. Parent to other subtypes of Shell, mostly empty
    """
    def __init__(self):
        pass

class StaticShell(Shell):
    """
    A Shell which forms a complete cylinderical layer around the centre of the object, without external flow.
    This can be a solid material (eg a layer of insulation or structural hull) or a layer of gas without forced flow.

    Args:
        orientation (str):        either 'horizontal' or 'vertical' to specify orientation of axis
        material (syrtis.Material): the material of the Shell
        radius_inner (float):       the internal radius of the Shell (m)
        thickness (float):          the thickness of the Shell (m)
        length (float):             the length of the Shell (m)
        external (bool) :           if True, the outer surface of the Shell is taken as the outer surface of the whole system. Default False
        thermal_resistance (float): an override for simple calculation of the Shell's thermal resistance. Ignored when not a positive number
    """

    def __init__(self, orientation, material, radius_inner, thickness, length, external=False, thermal_resistance=0):

        assert orientation == "horizontal" or orientation == "vertical", "'orientation' must be either 'horizontal' or 'vertical'"

        assert type(material) in material_classes, "'material' input is not valid Material-inherited object"

        assert isinstance(radius_inner, Number), "Shell 'internal_radius' must be a numerical value"
        assert radius_inner >= 0, "Shell 'internal_radius' must be a positive value"

        assert isinstance(thickness, Number), "Shell 'thickness' must be a numerical value"
        assert thickness > 0, "Shell 'thickness' must be a positive value"

        assert isinstance(length, Number), "Shell 'length' must be a numerical value"

        self.orientation = orientation
        self.material = material

        self.radius_inner = radius_inner
        self.thickness = thickness
        self.radius_outer = self.radius_inner + self.thickness

        self.length = length

        self.external = external

        if thermal_resistance == 0:
            self.calculate_thermal_resistance = True
            self.thermal_resistance = -1
        else:
            self.calculate_thermal_resistance = False
            self.thermal_resistance = thermal_resistance
        
    def __le__(self, other):
        if type(other) != StaticShell: raise TypeError

        return (self.radius_inner <= other.radius_inner)
    
    def __lt__(self, other):
        if type(other) != StaticShell: raise TypeError

        return (self.radius_inner < other.radius_inner)
    
    def __ge__(self, other):
        if type(other) != StaticShell: raise TypeError

        return (self.radius_inner >= other.radius_inner)
    
    def __gt__(self, other):
        if type(other) != StaticShell: raise TypeError

        return (self.radius_inner > other.radius_inner)
    
    def __eq__(self, other):
        if type(other) != StaticShell: raise TypeError

        return (self.radius_inner == other.radius_inner and self.radius_outer == other.radius_outer)

    def __ne__(self, other):
        if type(other) != StaticShell: raise TypeError

        return (self.radius_inner != other.radius_inner and self.radius_outer != other.radius_outer)

    def __repr__(self):
        return("A shell of {}, with an inner radius of {:.2f} and thickness of {:.3f} \n".format(
            self.material.name, self.radius_inner, self.thickness))

    def thermal_resistance(self, T_avg, T_delta, p, g):
        """
        Calculates the thermal resistance of the Shell
        Determines which equation to use and uses it
        """
        if self.calculate_thermal_resistance == False:
            pass
            
        elif self.radius_inner == 0:
            self.thermal_resistance = 0

        elif type(self.material) == Solid:
            self.thermal_resistance = self.thermal_resistance_solid()
        
        elif type(self.material) == ConstrainedIdealGas:
            self.thermal_resistance = self.thermal_resistance_annulus(T_avg, T_delta, p, g)
            

    def thermal_resistance_solid(self):
        # Using equation from Reference [1]
        R_th = np.log(self.radius_outer / self.radius_inner) / (2 * np.pi * self.material.k * self.length)

        return(R_th)

    def thermal_resistance_annulus(self, T_avg, T_delta, p, g):
        """
        Find the thermal resistance across a gas-filled annulus using natural convection correlations
        
        Args:
            T_avg (float):      Average temperature in the annulus (K)
            T_delta (float):    Temperature difference across the annulus (K)
            p (float):          Pressure in the annulus (Pa)
            g (float)           Gravitational acceleration (m/s2)
        """
        Ra = self.material.Ra(T_avg, T_delta, p, g, self.thickness)

        if self.orientation == "horizontal":
            # Using equation 9-56 from Reference [2] (via Raithby and Hollands, 1975)

            F_cyl = np.power(np.log(self.radius_outer/self.radius_inner), 4) / (
                np.power(self.thickness, 3) * np.power(
                    np.power(2 * self.radius_inner, -0.6) + np.power(2 * self.radius_outer, -0.6), 5))
            # Shape factor of the cylinders

            if F_cyl * Ra > 100:
                k_eq = (0.386 * np.power(self.material.Pr / (0.861 + self.material.Pr), 0.25) 
                * np.power(F_cyl * Ra, 0.25))
                k_eq = max(k_eq, 1)
            else:
                k_eq = 1
            
        elif self.orientation == "vertical":
            # Using equation 9-52 to 9-53 from Reference [2] from Berkovsky and Polevikov, 1977
            ratio = self.length / self.thickness

            if 1 <= ratio and 2 < ratio:
                k_eq = 0.18 * np.power((self.material.Pr * Ra) / (0.2 + self.material.Pr), 0.29)

                # Ensure model is working in the correct regime
                if (Ra * self.material.Pr) / (0.2 * self.material.Pr) < 1e3 and k_eq > 1:
                    print("Warning: heat transfer in vertical enclosure correlation is out of validation range. Proceed with caution")
                 
            elif 2 <= ratio < 10:
                k_eq = 0.22 * np.power((self.material.Pr * Ra) / (0.2 + self.material.Pr), 0.28) * np.power(ratio, -0.25)

                # Ensure model is working in the correct regime
                if Ra > 1e10 and k_eq > 1:
                    print("Warning: heat transfer in vertical enclosure correlation is out of validation range. Proceed with caution")
                
            elif 10 <= ratio:
                k_eq = 0.22 * np.power((self.material.Pr * Ra) / (0.2 + self.material.Pr), 0.28) * np.power(ratio, -0.25)

                if k_eq > 1:
                    print("Warning: heat transfer in vertical enclosure correlation is out of validation range. Proceed with caution")
            
            if k_eq < 1:
                k_eq = 1

        R_th = np.log(self.radius_outer / self.radius_inner) / (2 * np.pi * self.material.k * k_eq * self.length)

        return(R_th)

if __name__ == "__main__":
    co2 = ConstrainedIdealGas(44, 0.71, 10.9e-6, 749, 0.0153, p=580)
    steel = Solid(150, 8700, 500)

    s = StaticShell("horizontal", co2, 1, 0.2, 1)
    s.thermal_resistance(10, 9.81)

    print(s.thermal_resistance)
"""
Classes and inherited subclasses for Shells, the general term for layers of material included in the simulation.
In most cases, a Shell will be added as a subclass with more specialised properties

References:

 - [1] - https://en.wikipedia.org/wiki/Thermal_resistance
 - [2] - Y Cengel, Heat Transfer

"""

if __name__ == "__main__":
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
        configuration (str):        either 'horizontal' or 'vertical' to specify orientation of axis
        material (syrtis.Material): the material of the Shell
        radius_inner (float):       the internal radius of the Shell (m)
        thickness (float):          the thickness of the Shell (m)
        length (float):             the length of the Shell (m)
        external (bool) :           if True, the outer surface of the Shell is taken as the outer surface of the whole system. Default False
        thermal_resistance (float): an override for simple calculation of the Shell's thermal resistance. Ignored when not a positive number
    """

    def __init__(self, configuration, material, radius_inner, thickness, length, external=False, thermal_resistance=0):

        assert configuration == "horizontal" or configuration == "vertical", "'configuration' must be either 'horizontal' or 'vertical'"

        assert type(material) in material_classes, "'material' input is not valid Material-inherited object"

        assert isinstance(radius_inner, Number), "Shell 'internal_radius' must be a numerical value"
        assert radius_inner > 0, "Shell 'internal_radius' must be a positive value"

        assert isinstance(thickness, Number), "Shell 'thickness' must be a numerical value"
        assert thickness > 0, "Shell 'thickness' must be a positive value"

        assert isinstance(length, Number), "Shell 'length' must be a numerical value"

        self.configuration = configuration
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
        


    def thermal_resistance(self, T_avg, T_delta, g):
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
            self.thermal_resistance = self.thermal_resistance_wide_annulus(T_avg, T_delta, g)
            

    def thermal_resistance_solid(self):
        # Using equation from Reference [1]
        R_th = np.log(self.radius_outer / self.radius_inner) / (2 * np.pi * self.material.k * self.length)

        return(R_th)

    def thermal_resistance_wide_annulus(self, T_avg, T_delta, g):
        Ra = self.material.Ra(T_avg, T_delta, g, self.thickness)

        if self.configuration == "horizontal":
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
            
        elif self.configuration == "vertical":
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
    co2 = ConstrainedGas(210, 580, 0.71, 10.9e-6, 749, 8.74e-3, 0.0143)
    steel = Solid(150, 8700, 500)

    s = StaticShell("horizontal", co2, 1, 0.2, 1)
    s.thermal_resistance(10, 9.81)

    print(s.thermal_resistance)
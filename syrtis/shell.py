"""
Classes and inherited subclasses for Shells, the general term for layers of material included in the simulation.
In most cases, a Shell will be added as a subclass with more specialised properties

References:

 - [1] - https://en.wikipedia.org/wiki/Thermal_resistance
 - [2] - Y Cengel, Heat Transfer
 - [3] - Radiation View Factors, http://webserver.dmt.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf
 - [4] - Rapid Method for Determining Concentric Cylinder Radiation View Factors, Rea 1975

"""

from syrtis.material import *
import numpy as np

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

    def __init__(self, orientation, material, radius_inner, thickness, length, 
    external=False, thermal_resistance=0, parallel_thermal_resistance=0):

        assert orientation == "horizontal" or orientation == "vertical", "'orientation' must be either 'horizontal' or 'vertical'"

        assert type(material) in material_classes, "'material' input is not valid Material-inherited object"

        assert is_numeric(radius_inner, not_negative=True), "Shell 'internal_radius' must be a positive numerical value"

        assert is_numeric(thickness, positive=True), "Shell 'thickness' must be a positive numerical value"

        assert is_numeric(length, not_negative=True), "Shell 'length' must be a positive numerical value"

        self.orientation = orientation
        self.material = material

        self.radius_inner = radius_inner
        self.thickness = thickness
        self.radius_outer = self.radius_inner + self.thickness

        self.length = length

        self.external = external

        if thermal_resistance == 0:
            self.thermal_resistance_constant = 0
        else:
            self.thermal_resistance_constant = thermal_resistance
        
        self.parallel_thermal_resistance = parallel_thermal_resistance
        
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

    def calculate_thermal_resistance(self, T_avg, T_delta, g):
        """
        Calculates the thermal resistance of the Shell
        Determines which equation to use and uses it
        """
        if self.thermal_resistance_constant != 0:
            thermal_resistance = self.thermal_resistance_constant
            
        elif self.radius_inner == 0:
            thermal_resistance = 0

        elif type(self.material) == Solid:
            thermal_resistance = self.thermal_resistance_solid()
        
        elif type(self.material) == ConstrainedIdealGas:
            thermal_resistance = self.thermal_resistance_annulus(T_avg, T_delta, g)
        
        if self.parallel_thermal_resistance == 0:
            self.thermal_resistance = thermal_resistance
        else:
            if thermal_resistance == 0:
                self.thermal_resistance = thermal_resistance
            else:
                self.thermal_resistance = 1 / ((1 / thermal_resistance) + (1 / self.parallel_thermal_resistance))
            
    def thermal_resistance_solid(self):
        # Using equation from Reference [1]
        R_th = np.log(self.radius_outer / self.radius_inner) / (2 * np.pi * self.material.k * self.length)

        return(R_th)

    def thermal_resistance_annulus(self, T_avg, T_delta, g):
        """
        Find the thermal resistance across a gas-filled annulus using natural convection correlations
        
        Args:
            T_avg (float):      Average temperature in the annulus (K)
            T_delta (float):    Temperature difference across the annulus (K)
            p (float):          Pressure in the annulus (Pa)
            g (float)           Gravitational acceleration (m/s2)
        """
        Ra = self.material.Ra(T_avg, T_delta, g, self.thickness)
        Pr = self.material.Pr(T_avg)
        k = self.material.k(T_avg)

        if self.orientation == "horizontal":
            # Using equation 9-56 from Reference [2] (via Raithby and Hollands, 1975)

            F_cyl = np.power(np.log(self.radius_outer/self.radius_inner), 4) / (
                np.power(self.thickness, 3) * np.power(
                    np.power(2 * self.radius_inner, -0.6) + np.power(2 * self.radius_outer, -0.6), 5))
            # Shape factor of the cylinders

            if F_cyl * Ra > 100:
                k_eq = (0.386 * np.power(Pr / (0.861 + Pr), 0.25) 
                * np.power(F_cyl * Ra, 0.25))
                k_eq = max(k_eq, 1)
            else:
                k_eq = 1
            
        elif self.orientation == "vertical":
            # Using equation 9-52 to 9-53 from Reference [2] from Berkovsky and Polevikov, 1977
            ratio = self.length / self.thickness

            if 1 <= ratio and 2 < ratio:
                k_eq = 0.18 * np.power((Pr * Ra) / (0.2 + Pr), 0.29)

                # Ensure model is working in the correct regime
                if (Ra * Pr) / (0.2 * Pr) < 1e3 and k_eq > 1:
                    print("Warning: heat transfer in vertical enclosure correlation is out of validation range. Proceed with caution")
                 
            elif 2 <= ratio < 10:
                k_eq = 0.22 * np.power((Pr * Ra) / (0.2 + Pr), 0.28) * np.power(ratio, -0.25)

                # Ensure model is working in the correct regime
                if Ra > 1e10 and k_eq > 1:
                    print("Warning: heat transfer in vertical enclosure correlation is out of validation range. Proceed with caution")
                
            elif 10 <= ratio:
                k_eq = 0.22 * np.power((Pr * Ra) / (0.2 + Pr), 0.28) * np.power(ratio, -0.25)

                if k_eq > 1:
                    print("Warning: heat transfer in vertical enclosure correlation is out of validation range. Proceed with caution")
            
            else:
                k_eq = 1
            
            if k_eq < 1:
                k_eq = 1

        R_th = np.log(self.radius_outer / self.radius_inner) / (2 * np.pi * k * k_eq * self.length)

        return(R_th)

    def thermal_energy(self, endcap_type, T_avg, T_ref=0):
        """
        Calculates the thermal energy (J) stored in the Shell relative to T_ref
        
        Args:
            endcap_type (str):      either "hemisphere" or "flat"
            T_avg (float):          average temperature of the Shell (K)
            T_ref (float):          reference temperature for the energy state (K)
        """

        volume = np.pi * (np.power(self.radius_outer, 2) - np.power(self.radius_inner, 2)) * self.length

        if endcap_type == "hemisphere":
            volume += (4/3) * np.pi * (np.power(self.radius_outer, 3) - np.power(self.radius_inner, 3))
        elif endcap_type == "flat":
            volume += 2 * np.pi * np.power(self.radius_outer, 2) * self.thickness
        
        if type(self.material) == Solid:
            thermal_capacity = self.material.rho * self.material.cp * volume
        elif type(self.material) == ConstrainedIdealGas:
            thermal_capacity = self.material.rho(T_avg) * self.material.cp(T_avg) * volume
        
        thermal_energy = thermal_capacity * (T_avg - T_ref)

        return(thermal_energy)

class GroundLevel(Shell):
    """
    An object to store the ground, whenever the Habitat is in contact with the ground

    Args:
        habitat_axis_height (float):        height of the habitat central axis above ground level, for horizontal orientation
                                            height of the lower end of the cylinder above ground level, for vertical orientation
                                            Defaults to a large value, indicating no thermal conduction contact (m)
        thermal_resistance (float):         thermal resistance between Habitat outer wall, default to 0 = no thermal contact (K/W)
    
    """
    def __init__(self, habitat_axis_height=1e3, thermal_resistance=0):
        
        assert is_numeric(habitat_axis_height), "GroundLevel 'height_below_habitat' must be a numerical value"

        assert is_numeric(thermal_resistance, not_negative=True), "GroundLevel 'thermal_resistance' must be a positive numerical value"

        self.habitat_axis_height = habitat_axis_height
        self.thermal_resistance = thermal_resistance

class Earthworks(Shell):
    """
    An object to store cylinderical hole in the ground in which the habitat sits. 
    Can include any large open cavity in the ground, including trenches or holes
    Also includes tunnels or lava tubes, natural or otherwise

    TODO: find a better name. Earthworks? Excavations? Cavities?

    Args:
        radius_inner (float):       radius of the inside of the of the cavity (m)
        depth_of_axis (float):      depth below ground level that the central axis of the cavity sits (m) 
                                    if less than radius_inner, the cavity will not be entirely below ground
    
    """
    def __init__(self, radius_inner, depth_of_axis):
        assert is_numeric(radius_inner, positive=True), "Cavity 'radius_inner' must be a positive number"
        assert is_numeric(depth_of_axis), "Cavity 'depth_of_axis' must be a number"
        assert depth_of_axis > -radius_inner, "Cavity 'depth_of_axis' is too low, cavity does not intersect the ground"

        self.radius_inner = radius_inner
        self.depth_of_axis = depth_of_axis

        if self.depth_of_axis > self.radius_inner:
            self.exit_strip_width = 0
        else:
            self.exit_strip_width = 2 * self.radius_inner * np.sin(
                np.arccos(self.depth_of_axis / self.radius_inner))
    
    def view_factor_ground_cylinder(self, orientation, radius_outer, length_outer, axis_height_from_ground):
        """
        Calculates the view factor from a cylinder inside the Earthworks to the ground 
            (including the inner surface of the earthworks)
        
        Args:
            orientation (str):                  one of "horizontal" or "vertical"
            radius_outer (float):               outer radius of the cylinder (m)
            length_outer (float):               length of the cylinder
            axis_height_from_ground (float):    height of the cylinder axis (horizontal) or bottom of cylinder (vertical)
                                                above the ground
        """

        view_factor_ground = 0

        if self.depth_of_axis > self.radius_inner:
            # Earthworks has no view to the sky
            view_factor_ground = 1
        
        elif orientation == "horizontal":
            # Using Equation from [3] for view factor from planar strip to cylinder
            v = self.exit_strip_width / (2 * radius_outer)
            h = (self.radius_inner + self.depth_of_axis - axis_height_from_ground) / radius_outer

            if h > 0:
                # The top of the habitat is below the top rim of the hole
                # An virtual surface is drawn across the rim
                view_factor_ground = 1 - (np.arctan(v / h) / np.pi)
            
            else:
                # The top of the habitat extends above the top rim of the hole
                # Virtual surfaces are drawn above and to the sides of the cylinder

                v_top = 2 * radius_outer
                h_top = 1
                view_factor_sky_top = np.arctan(v_top / h_top) / np.pi

                # Side boxes
                w1 = self.radius_inner + self.depth_of_axis - axis_height_from_ground
                w2 = radius_outer
                h_side = 1

                v1_side = w1 / radius_outer
                v2_side = w2 / radius_outer

                view_factor_sky_side = (np.arctan(v2_side / h_side) - np.arctan(v1_side / h_side)) / np.pi

                view_factor_ground = 1 - (view_factor_sky_top + view_factor_sky_side)

        elif orientation == "vertical":

            exit_strip_radius = self.exit_strip_width / 2
            D = ((self.radius_inner + self.depth_of_axis) - (axis_height_from_ground + length_outer)) / exit_strip_radius
            Y = 1e10 / exit_strip_radius
            X = ((self.radius_inner + self.depth_of_axis) - axis_height_from_ground) / exit_strip_radius
            L = length_outer / exit_strip_radius
            R = radius_outer / exit_strip_radius

            def A(x):
                return(np.power(x, 2) + np.power(R, 2) - 1) 
            
            def B(x):
                return(np.power(x, 2) - np.power(R, 2) + 1)
            
            def F(x):
                f = np.sqrt(np.power((A(x) + 2) / R, 2 ) - 4) * np.arccos(A(x) * x / B(x)) - (
                    A(x) / (2 * x * R)) * np.arcsin(R)
                f *= (-0.5 / x) 
                f += B(x) / (8 * R * x)

                return(f)

            if self.radius_inner + self.depth_of_axis > (axis_height_from_ground + length_outer):
                # The top of the habitat is below the top rim of the hole
                # A virtual cylinder is drawn at the top of the rim
                # Equation 8 from [4] is used, with the length of the outer cylinder tending to infinity
                view_factor_sky = ((L + D) / D) * F(L + D)
                view_factor_sky += ((Y + D) / L) * F(Y + D)
                view_factor_sky -= (D / L) * F(D)
                view_factor_sky -= ((Y + D + L) / L) * F(Y + D + L)
            
            else:
                # The top of the habitat is above the top rim of the hole
                # A virtual cylinder is drawn at the top of the rim
                # Equation 6 from [4] is used, with the length of the outer cylinder tending to infinity

                view_factor_sky = (X / L) * F(X)
                view_factor_sky += ((L - X) / L) * (1 - F(L - X))
                view_factor_sky += ((Y + X - L) / L) * F(Y + X - L)
                view_factor_sky -= ((X + Y) / L) * F(X + Y)
            
            view_factor_ground = 1 - view_factor_sky

            
        return(view_factor_ground)

    def view_factor_ground_hemisphere(self, orientation, radius_outer, length_outer, axis_height_from_ground):
        """
        Calculates the view factor from a hemispherical endcap inside the Earthworks to the ground 
            (including the inner surface of the earthworks)
        
        Args:
            orientation (str):                  one of "horizontal" or "vertical"
            radius_outer (float):               outer radius of the cylinder (m)
            length_outer (float):               length of the cylinder
            axis_height_from_ground (float):    height of the cylinder axis (horizontal) or bottom of cylinder (vertical)
                                                above the ground
        """
        view_factor_ground = 0

        if self.depth_of_axis > self.radius_inner:
            # Earthworks has no view to the sky
            view_factor_ground = 1
        
        return(view_factor_ground)
    
    def view_factor_ground_disc(self, orientation, radius_outer, length_outer, axis_height_from_ground):
        """
        Calculates the view factor from a flat endcap inside the Earthworks to the ground 
            (including the inner surface of the earthworks)
        
        Args:
            orientation (str):                  one of "horizontal" or "vertical"
            radius_outer (float):               outer radius of the cylinder (m)
            length_outer (float):               length of the cylinder
            axis_height_from_ground (float):    height of the cylinder axis (horizontal) or bottom of cylinder (vertical)
                                                above the ground
        """
        view_factor_ground = 0

        if self.depth_of_axis > self.radius_inner:
            # Earthworks has no view to the sky
            view_factor_ground = 1
        
        return(view_factor_ground)
    
    def direct_solar_intensity_scaling(self, orientation, radius_outer, length_outer, axis_height_from_ground, solar_altitude, solar_azimuth):
        """
        Scales the direct solar intensity based on how much of the habitat can see the Sun at a given value of altitude and azimuth
        Passes value back to solar heat gain functions in Habitat 

        Args:
            orientation (str):                  one of "horizontal" or "vertical"
            radius_outer (float):               outer radius of the cylinder (m)
            length_outer (float):               length of the cylinder
            axis_height_from_ground (float):    height of the cylinder axis (horizontal) or bottom of cylinder (vertical)
                                                above the ground
            solar_altitude (float):             vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):              horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
                                                assumed that a horizontal habitat lies in the same axis as the Earthworks tunnel
        
        TODO: work out the geometry for this function. The percentage of the habitat which is hit with direct sunlight.
            Probably requires working out the depth into the Earthworks which the sunlight reaches, and what fraction of the
            Habitat is above a given angle.
            Can't currently wrap head around how the solar azimuth affects the depth reached. When azimuth=0, sunlight will reach
            the bottom. When azimuth=90, it will be the direct projection from the rim. What about in between?
        """
        solar_multiplier = 1

        if self.depth_of_axis > self.radius_inner:
            # Earthworks has no view to the sky
            solar_multiplier = 1
            return(solar_multiplier)

        angle_over_rim = 0
        
        if orientation == "horizontal":
            pass
        
        elif orientation == "vertical":
            pass
        
        return(solar_multiplier)
    
    def indirect_solar_intensity_scaling(self, orientation, radius_outer, length_outer, axis_height_from_ground, solar_altitude, solar_azimuth):
        """
        Scales the indirect solar intensity based on how much of the habitat can see the Sun at a given value of altitude and azimuth
        Passes value back to solar heat gain functions in Habitat 

        Args:
            orientation (str):                  one of "horizontal" or "vertical"
            radius_outer (float):               outer radius of the cylinder (m)
            length_outer (float):               length of the cylinder
            axis_height_from_ground (float):    height of the cylinder axis (horizontal) or bottom of cylinder (vertical)
                                                above the ground
            solar_altitude (float):             vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):              horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
                                                assumed that a horizontal habitat lies in the same axis as the Earthworks tunnel

        TODO: work out the geometry required for this function. The percentage of the surroundings lit with direct sunlight
        """
        solar_multiplier = 1

        if self.depth_of_axis > self.radius_inner:
            # Earthworks has no view to the sky
            solar_multiplier = 1
            return(solar_multiplier)
        
        if orientation == "horizontal":
            pass
        
        elif orientation == "vertical":
            pass
        
        return(solar_multiplier)

"""
Top-level object and tools for the whole Habitat, composed of individual Shells

References:

 - [1] - Y Cengel, Heat Transfer
 - [2] - Sky temperature modelisation and applications in building simulation, Adelard et al 1998
 - [3] - Long-wave  radiation from clear skies, Swinbank 1963
 - [4] - Radiation View Factors, http://webserver.dmt.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf
 - [5] - Thermal control of MSL Rover "Curiosity" using an Active Fluid Loop, Birur 2013
 - [6] - The thermal control system of NASA’s Curiosity rover: a case study, Quattrocchi et al, 2022 
"""

from gettext import translation
import numpy as np
from numbers import Number


if __name__ == "__main__":
    from shell import *
    #from material import *
else:
    #from shell import *
    from syrtis.shell import *
    from syrtis.material import *

endcap_types = ["hemisphere", "flat"]

class Habitat:
    """
    Object to store the whole habitat geometry

    Args:
        orientation (str):        either 'horizontal' or 'vertical' to specify orientation of axis
        length (float):             the length of the central axis of the Habitat (m)    
    """

    def __init__(self, orientation, length, endcap_type):
        assert orientation == "horizontal" or orientation == "vertical", "'orientation' must be either 'horizontal' or 'vertical'"

        assert isinstance(length, Number), "Habitat 'length' must be a numerical value"
        assert length > 0, "Habitat 'length' must be a positive value"

        assert endcap_type in endcap_types, "'endcap_type' must be a valid keyword"

        self.orientation = orientation
        self.length = length

        self.endcap_type = endcap_type

        self._shells = []

        self.radius_outer = 0

        self.verified = False
    
    def create_static_shell(self, material, thickness, thermal_resistance=0):
        """
        Create a StaticShell that conforms around the outside of the outermost one
        """
        
        new_shell = StaticShell(self.orientation, material, self.radius_outer, 
        thickness, self.length, False, thermal_resistance)

        self.radius_outer += thickness
        self.verified = False
        self._shells.append(new_shell)
    
    def append_shell(self, shell):
        """
        Append a Shell of any type to the Habitat, checking for no overlaps
        """

        self.verified = False
        
        self._shells.append(shell)

        self.radius_outer = max(self.radius_outer, shell.radius_outer)
    
    def verify_geometry(self):
        """
        Checks the geometry of the habitat for errors
        """
        
        self._shells.sort()

        assert len(self._shells) >= 2, "Habitat requires at least two Shells, one enclosing the other"
        
        # Do any of the layers overlap, or do they have any gaps?
        
        overlap_error = False
        gap_error = False
        length_error = False
        for shell_count, shell in enumerate(self._shells):

            if shell_count != 0:
                # Check the shell inside this one to ensure a snug fit

                if round(self._shells[shell_count - 1].radius_outer, 8) > round(shell.radius_inner, 8):
                    # The Shell overlaps the inner edge of the current one
                    overlap_error = True
                    break

                if round(self._shells[shell_count - 1].radius_outer, 8) < round(shell.radius_inner, 8):
                    # The Shell overlaps the inner edge of the current one
                    gap_error = True
                    break
                    
                if round(self._shells[shell_count - 1].length, 8) > round(shell.length, 8):
                    length_error = True
                    break
            
            else:
                assert round(shell.radius_inner, 8) == 0, "There is no Shell at the centre of the Habitat. It must be present for a full calculation. Use a ConstrainedIdealGas Shell to simulate the pressurised space."
        
        assert not overlap_error, "Some Shells overlap each other"
        assert not gap_error, "Some Shells have gaps betweem them"
        assert not length_error, "Some Shells become shorter, not longer, moving away from the centre of the Habitat"
        

        # Is the outermost layer a solid and set to external?

        assert type(self._shells[-1].material) == Solid, "Outermost shell must be a Solid"

        for shell in self._shells:
            shell.external = False

        self._shells[-1].external = True

        self.exposed_convective_area()

        self.verified = True
    
    def exposed_convective_area(self):
        """
        Find the area exposed to convective losses, endcap and cylinder.
        Called during geometry verification
        """
        self.exposed_area_cylinder = 2 * np.pi * self._shells[-1].radius_outer * self._shells[-1].length
        
        if self.orientation == "horizontal":
            
            if self.endcap_type == "hemisphere":
                self.exposed_area_endcap = 4 * np.pi * np.power(self._shells[-1].radius_outer, 2)
            
            elif self.endcap_type == "flat":
                self.exposed_area_endcap = 2 * np.pi * np.power(self._shells[-1].radius_outer, 2)
        
        elif self.orientation == "vertical":

            if self.endcap_type == "hemisphere":
                self.exposed_area_endcap = 2 * np.pi * np.power(self._shells[-1].radius_outer, 2)
            
            elif self.endcap_type == "flat":
                self.exposed_area_endcap = 1 * np.pi * np.power(self._shells[-1].radius_outer, 2)

    def exposed_radiative_area(self, solar_altitude, solar_azimuth):
        """
        Find the area of the habitat in direct and indirect sunlight
        Called by solar_gain_direct and solar_gain_indirect to find areas
        
        Args:
            solar_altitude (float): vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):  horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
        """
        solar_altitude_rad = np.deg2rad(solar_altitude)
        solar_azimuth_rad = np.deg2rad(solar_azimuth)

        direct_solar_area = 0
        indirect_solar_area = 0

        """
        Direct solar area - area that can see the Sun
        """

        if self.orientation == "horizontal":
            # Area of cylinder projected to the plane perpendicular to the Sun
            # Two components are axial and radial-ish directions
            direct_solar_area += 2 * self._shells[-1].length * self._shells[-1].radius_outer * (
                np.cos(solar_altitude_rad) * np.cos(solar_azimuth_rad) 
                + np.sin(solar_azimuth_rad))
            
            if self.endcap_type == "flat":
                direct_solar_area += np.pi * np.power(self._shells[-1].radius_outer, 2) * (
                    np.cos(solar_altitude_rad) * abs(np.cos(solar_azimuth_rad)))
            
            elif self.endcap_type == "hemisphere":
                # Area of hemisphere projected onto the plane perpendicular to the Sun
                # Two components are axial (plan view) and radial-ish (side view)
                direct_solar_area += np.pi * np.power(self._shells[-1].radius_outer, 2) * (
                    np.cos(solar_altitude_rad) + 0.5 * np.sin(solar_altitude_rad)) * abs(np.cos(solar_azimuth_rad))
                         
        elif self.orientation == "vertical":
            # Area of cylinder projected to the plane perpendicular to the Sun
            direct_solar_area += 2 * self._shells[-1].length * self._shells[-1].radius_outer * (
            np.sin(solar_azimuth_rad))

            if self.endcap_type == "flat":
                direct_solar_area += self.exposed_area_endcap * np.cos(solar_altitude_rad)
            
            elif self.endcap_type == "hemisphere":
                # Area of hemisphere projected onto the plane perpendicular to the Sun
                # Two components are axial (plan view) and radial-ish (side view)
                direct_solar_area += np.pi * np.power(self._shells[-1].radius_outer, 2) * (
                    np.cos(solar_altitude_rad) + 0.5 * np.sin(solar_altitude_rad))
        
        """
        Indirect solar area - area that can see the ground
        """
        if self.orientation == "horizontal":
            indirect_solar_area += 2 * self._shells[-1].length * self._shells[-1].radius_outer

            if self.endcap_type == "flat":
                indirect_solar_area += 2 * np.pi * np.power(self._shells[-1].radius_outer, 2)
            
            elif self.endcap_type == "hemisphere":
                indirect_solar_area += 4 * np.pi * np.power(self._shells[-1].radius_outer, 2)
        
        elif self.orientation == "vertical":
            indirect_solar_area += 2 * self._shells[-1].length * self._shells[-1].radius_outer

            if self.endcap_type == "flat":
                pass

            elif self.endcap_type == "hemisphere":
                indirect_solar_area += 2 * np.pi * np.power(self._shells[-1].radius_outer, 2)
        
        return(direct_solar_area, indirect_solar_area)

    def build_thermal_resistances(self, shell_temperatures, g):
        """
        Outputs the thermal resistances of each shell layer
        
        Args:
            shell_temperatures (list of floats):    list of len(self._shells)+1, referring to temperatures at the boundary
                                                    of each material. 0th is the inner wall temperature of innermost shell,
                                                    1st is outer wall of innermost/inner wall of second.
            g (float):                              gravity
        """
        assert len(shell_temperatures) == len(self._shells)+1, "length of 'shell_temperatures' must be equal to number of Shells plus one"

        shell_thermal_resistances = np.zeros((len(self._shells)))

        for shell_count, shell in enumerate(self._shells):
            T_avg = (shell_temperatures[shell_count] + shell_temperatures[shell_count+1]) / 2
            T_delta = shell_temperatures[shell_count+1] - shell_temperatures[shell_count]

            shell.calculate_thermal_resistance(T_avg, T_delta, g)

            shell_thermal_resistances[shell_count] = shell.thermal_resistance
        
        return(shell_thermal_resistances)

    def placeholder_convective_loss(self, T_wall, T_air):
        """
        Placeholder constant-h model for testing purposes
        """
        h = 1
        A = 2 * np.pi * self.length * self.radius_outer

        R_th = 1 / (h * A)

        Q = (T_wall - T_air) / R_th
        
        return(Q)
    
    def nusselt_sphere(self, air, v_air, T_air, T_wall, D):
        """
        Find the Nusselt number for a sphere in uniform flow
        Whitaker correlation, all properties evaulated at freestream temperature
        From Reference [2], Equation (7-36)

        Args:
            air (ConstrainedIdealGas):  object for the external flow
            v_air (float):              velocity of airflow (m/s)
            T_air (float):              temperature of the freestream (K)
            T_wall (float):             temperature of the wall (K)
            D (float):                  diameter (m)
        """
        Re = air.Re(T_air, D, v_air)
        Pr = air.Pr(T_air)
        mu_air = air.mu(T_air)
        mu_wall = air.mu(T_wall)

        if 0.7 >= Pr or Pr >= 380:
            print("Warning: heat transfer on hemispherical endcap correlation is out of validation range. Proceed with caution")
        
        if 3.5 >= Re or Re >= 80000:
            print("Warning: heat transfer on hemispherical endcap correlation is out of validation range. Proceed with caution")

        Nu_D = 2 + (
            (0.4 * np.power(Re, 0.4) + 0.06 * np.power(Re, 2/3)) * np.power(Pr, 0.4) * 
            np.power(mu_air / mu_wall, 0.25))
        
        return(Nu_D)
    
    def nusselt_plate_crossflow(self, air, v_air, T_air, T_wall, D):
        """
        Find the Nusselt number for crossflow over a plate
        Does the correct checking for laminar, turbulent etc
        
        Args:
            air (ConstrainedIdealGas):  object for the external flow
            v_air (float):              velocity of airflow (m/s)
            T_air (float):              temperature of the freestream (K)
            T_wall (float):             temperature of the wall (K)
            D (float):                  length scale (m)
        """
        transition_Re = 5e5

        T_film = (T_wall + T_air) / 2

        Re = air.Re(T_film, D, v_air)
        Pr = air.Pr(T_film)

        if Re < transition_Re:
            # Laminar correlation, Reference [2] Equation (7-21)
            Nu_D = 0.664 * np.power(Re, 0.5) * np.power(Pr, 1/3)
        
        elif Re > transition_Re and Re < 1e7:
            # Use partially-laminar correlation, Reference [2] Equation (7-24)
            Nu_D = (0.037 * np.power(Re, 0.8) - 871) * np.power(Pr, 1/3)
        
        else:
            print("Warning: heat transfer on flat plate correlation is out of validation range. Proceed with caution")
            Nu_D = (0.037 * np.power(Re, 0.8) - 871) * np.power(Pr, 1/3)


        return(Nu_D)

    def convective_loss_endcap_cross(self, air, v_air, T_air, T_wall):
        """
        Convective heat loss from the endcaps with uniform crossflow
        
        Args:
            air (ConstrainedIdealGas):  object for the external flow
            v_air (float):              velocity of airflow (m/s)
            T_air (float):              temperature of the freestream (K)
            T_wall (float):             temperature of the wall (K)
        """
        D = self._shells[-1].radius_outer * 2

        if self.endcap_type == "hemisphere":
            Nu_D = self.nusselt_sphere(air, v_air, T_air, T_wall, D)
            
        elif self.endcap_type == "flat":
            Nu_D = self.nusselt_plate_crossflow(air, v_air, T_air, T_wall, D)
            
        h = Nu_D * air.k((T_air + T_wall) / 2) / D

        Q_conv = h * self.exposed_area_endcap * (T_wall - T_air)

        return(Q_conv)
    
    def convective_loss_endcap_axial(self, air, v_air, T_air, T_wall):
        """
        Convective heat loss from the endcaps with uniform crossflow
        
        Args:
            air (ConstrainedIdealGas):  object for the external flow
            v_air (float):              velocity of airflow (m/s)
            T_air (float):              temperature of the freestream (K)
            T_wall (float):             temperature of the wall (K)
        """
        D = self._shells[-1].radius_outer * 2

        if self.endcap_type == "hemisphere":
            Nu_D = self.nusselt_sphere(air, v_air, T_air, T_wall, D)
        
        elif self.endcap_type == "flat":
            if self.orientation == "horizontal":
                T_film = (T_wall + T_air) / 2

                Re = air.Re(T_film, D, v_air)
                Pr = air.Pr(T_film)
                if 4000 >= Re or Re >= 15000:
                    print("Warning: heat transfer on flat endcap correlation is out of validation range. Proceed with caution")
                
                Nu_D = 0.228 * np.power(Re, 0.731) * np.power(Re, 1/3)
            
            elif self.orientation == "vertical":
                Nu_D = self.nusselt_plate_crossflow(v_air, T_air, T_wall, D)

        h = Nu_D * air.k((T_air + T_wall) / 2) / D

        Q_conv = h * self.exposed_area_endcap * (T_wall - T_air)

        return(Q_conv)

    def convective_loss_cylinder_cross(self, air, v_air, T_air, T_wall):
        """
        Convective heat loss from a cylinder with uniform crossflow.

        Args:
            air (ConstrainedIdealGas):  object for the external flow
            v_air (float):              velocity of airflow (m/s)
            T_air (float):              temperature of the freestream (K)
            T_wall (float):             temperature of the wall (K)
        """
        T_film = (T_wall + T_air) / 2

        D = self._shells[-1].radius_outer * 2

        Re = air.Re(T_film, D, v_air)
        Pr = air.Pr(T_film)

        # Churchill-Bernstein correlation, from Reference [1] Equation (7-35)

        Nu_D = 0.3 + (((0.62 * np.power(Re, 0.5) * np.power(Pr, 1/3.)) / (np.power(1 + np.power(0.4 / Pr, 2/3), 0.25))) *
        np.power(1 + np.power(Re / 282000, 0.625), 0.8))

        h = Nu_D * air.k(T_film) / D

        Q = h * self.exposed_area_cylinder * (T_wall - T_air)

        return(Q)

    def convective_loss_cylinder_axial(self, air, v_air, T_air, T_wall):
        """
        Convective heat loss from a cylinder with uniform axial flow.

        Args:
            air (ConstrainedIdealGas):  object for the external flow
            v_air (float):              velocity of airflow (m/s)
            T_air (float):              temperature of the freestream (K)
            T_wall (float):             temperature of the wall (K)
        """
        T_film = (T_wall + T_air) / 2

        L = self._shells[-1].length

        Nu_D = self.nusselt_plate_crossflow(air, v_air, T_air, T_wall, L)

        h = Nu_D * air.k(T_film) / L

        Q_conv = h * self.exposed_area_cylinder * (T_wall - T_air)

        return(Q_conv)
    
    def sky_temperature(self, T_air):
        """
        Calculate the sky temperature for radiative heat transfer
        Uses equation from References [2] and [3], validated against data from Reference [5]

        Args: 
            T_air:      air temperature (K)
        """

        T_sky = 0.0552 * np.power(T_air, 1.5)

        return(T_sky)

    def view_factor_ground(self):
        """
        Calculate the view factor to the ground from the cylinder
        Uses equations from Reference [4]
        """

        if self.orientation == "vertical":
            vf_cylinder = 0.5

            vf_endcap = 0.5

            vf = ((vf_cylinder * self.exposed_area_cylinder) +
            (vf_endcap * self.exposed_area_endcap)) / (self.exposed_area_cylinder + self.exposed_area_endcap)

        elif self.orientation == "horizontal":

            vf_cylinder_top = 0.5 - (1 / np.pi)
            vf_cylinder_bottom = 0.5 + (1 / np.pi)

            vf_endcap = 0.5

            vf = ((vf_cylinder_top * 0.5 * self.exposed_area_cylinder) + 
            (vf_cylinder_bottom * 0.5 * self.exposed_area_cylinder) +
            (vf_endcap * self.exposed_area_endcap)) / (self.exposed_area_cylinder + self.exposed_area_endcap)
        
        return(vf)

    def radiative_loss_sky(self, T_wall):
        """
        Calculate radiative loss to the sky.

        Args:
            T_wall (float):             temperature of the wall (K)
        """
        
        vf_sky = 1 - self.view_factor_ground()

        Q_sky = 5.67e-8 * self._shells[-1].material.emit * np.power(T_wall, 4) * (
            self.exposed_area_cylinder + self.exposed_area_endcap) * vf_sky
        
        return(Q_sky)
    
    def radiative_gain_sky(self, T_air):
        """
        Calculate radiative loss to the sky.

        Args:
            T_air (float):              temperature of the ambient surrounding air (K)
        """

        T_sky = self.sky_temperature(T_air)
        
        vf_sky = 1 - self.view_factor_ground()
        
        Q_sky = 5.67e-8 * self._shells[-1].material.absorb * np.power(T_sky, 4) * (
            self.exposed_area_cylinder + self.exposed_area_endcap) * vf_sky

        # A negative sign is used for consistency with convention that +ve Q = heat loss
        return(Q_sky)
    
    def radiative_loss_ground(self, T_wall):
        """
        Calculate radiative loss to the sky.

        Args:
            T_wall (float):             temperature of the wall (K)
        """

        vf_ground = 1 - self.view_factor_ground()

        Q_ground = 5.67e-8 * self._shells[-1].material.emit * np.power(T_wall, 4) * (
            self.exposed_area_cylinder + self.exposed_area_endcap) * vf_ground

        return(Q_ground)
    
    def radiative_gain_ground(self, T_ground):
        """
        Calculate radiative loss to the sky.

        Args:
            T_ground (float):              temperature of the ground(K)
        """

        vf_ground = 1 - self.view_factor_ground()

        Q_ground = 5.67e-8 * self._shells[-1].material.absorb * np.power(T_ground, 4) * (
            self.exposed_area_cylinder + self.exposed_area_endcap) * vf_ground

        # A negative sign is used for consistency with convention that +ve Q = heat loss
        return(-Q_ground)

    def solar_gain_direct(self, solar_altitude, solar_azimuth, solar_intensity):
        """
        Calculate heat gain to the outermost layer of the habitat, or into the centre if all outer layers have sufficient transparency
        
        Args:
            solar_altitude (float): vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):  horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
            solar_intensity (float):power delivered by solar radiation after dust absorption, W/m2
        """

        direct_lit_area, indirect_lit_area = self.exposed_radiative_area(solar_altitude, solar_azimuth)
        
        if self._shells[-1].material.transmit == 0:
            # Outer layer is opaque - all solar energy is delivered to outermost shell
            Q_solar_direct = solar_intensity * direct_lit_area * self._shells[-1].material.absorb
        
        # A negative sign is used for consistency with convention that +ve Q = heat loss
        return(-Q_solar_direct)
    
    def solar_gain_indirect(self, solar_altitude, solar_azimuth, solar_intensity, albedo_ground):
        """
        Calculate heat gain to the outermost layer of the habitat, or into the centre if all outer layers have sufficient transparency
        
        Args:
            solar_altitude (float): vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):  horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
            solar_intensity (float):power delivered by solar radiation after dust absorption, W/m2
            albedo_ground (float):  albedo of the surface around the habitat
        """

        direct_lit_area, indirect_lit_area = self.exposed_radiative_area(solar_altitude, solar_azimuth)
        
        if self._shells[-1].material.transmit == 0:
            # Outer layer is opaque - all solar energy is delivered to outermost shell
            Q_solar_indirect = solar_intensity * albedo_ground * indirect_lit_area * self._shells[-1].material.absorb
        
        # A negative sign is used for consistency with convention that +ve Q = heat loss
        return(Q_solar_indirect)


    




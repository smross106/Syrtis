"""
Top-level object and tools for the whole Habitat, composed of individual Shells

References:

 - [1] - Y Cengel, Heat Transfer

"""

from gettext import translation
import numpy as np
from numbers import Number


if __name__ == "__main__":
    from shell import *
    #from material import *
else:
    from shell import *
    #from syrtis.shell import *
    #from syrtis.material import *

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

        self.exposed_area_convection_cylinder = 2 * np.pi * self._shells[-1].radius_outer * self._shells[-1].length
        
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

        self.verified = True
    
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

        Q = h * self.exposed_area_endcap * (T_wall - T_air)

        return(Q)
    
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

        Q = h * self.exposed_area_endcap * (T_wall - T_air)

        return(Q)

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

        Q = h * self.exposed_area_convection_cylinder * (T_wall - T_air)

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

        Q = h * self.exposed_area_convection_cylinder * (T_wall - T_air)

        return(Q)







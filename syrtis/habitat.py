"""
Top-level object and tools for the whole Habitat, composed of individual Shells

References:

 - [1] - Y Cengel, Heat Transfer

"""

import numpy as np
from numbers import Number


if __name__ == "__main__":
    from shell import *
    #from material import *
else:
    from shell import *
    #from syrtis.shell import *
    #from syrtis.material import *

class Habitat:
    """
    Object to store the whole habitat geometry

    Args:
        orientation (str):        either 'horizontal' or 'vertical' to specify orientation of axis
        length (float):             the length of the central axis of the Habitat (m)    
    """

    def __init__(self, orientation, length):
        assert orientation == "horizontal" or orientation == "vertical", "'orientation' must be either 'horizontal' or 'vertical'"

        assert isinstance(length, Number), "Habitat 'length' must be a numerical value"
        assert length > 0, "Habitat 'length' must be a positive value"

        self.orientation = orientation
        self.length = length

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
            
            else:
                assert round(shell.radius_inner, 8) == 0, "There is no Shell at the centre of the Habitat. It must be present for a full calculation. Use a ConstrainedIdealGas Shell to simulate the pressurised space."
        
        assert not overlap_error, "Some Shells overlap each other"
        assert not gap_error, "Some Shells have gaps betweem them"
        

        # Is the outermost layer a solid and set to external?

        assert type(self._shells[-1].material) == Solid, "Outermost shell must be a Solid"

        for shell in self._shells:
            shell.external = False

        self._shells[-1].external = True

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
    
    def convective_loss_cylinder_cross(self, air, v_air, T_air, T_wall):
        """
        Convective heat loss from a cylinder with uniform crossflow.

        Args:
            air (ConstrainedIdealGas):  object for the external flow
        """

        T_film = (T_wall + T_air) / 2

        D = self._shells[-1].radius_outer * 2

        Re = air.Re(T_film, D, v_air)
        Pr = air.Pr(T_film)

        # Churchill-Bernstein correlation, from Reference [1] Equation (7-35)

        Nu_D = 0.3 + (((0.62 * np.power(Re, 0.5) * np.power(Pr, 1/3.)) / (np.power(1 + np.power(0.4 / Pr, 2/3), 0.25))) *
        np.power(1 + np.power(Re / 282000, 0.625), 0.8))

        h = Nu_D * air.k(T_film) / D

        Q = h * (2 * np.pi * D * self._shells[-1].length) * (T_wall - T_air)

        return(Q)









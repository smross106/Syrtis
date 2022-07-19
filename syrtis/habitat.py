"""
Top-level object and tools for the whole Habitat, composed of individual Shells
"""

import numpy as np


if __name__ == "__main__":
    from shell import *
    #from material import *
else:
    from syrtis.shell import *
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
        
        new_shell = StaticShell(self.orientation, material, self.radius_outer, thickness, self.length, False, thermal_resistance)

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
                    
                
            """for i in self._shells[shell_count+1:]:
                
                # All shells after the current one. This reduces compute time
                if round(i.radius_outer, 8) > inner and round(i.radius_inner, 8) < inner:
                    # The Shell overlaps the inner edge of the current one
                    overlap_error = True
                    break
                
                elif round(i.radius_inner, 8) < outer and round(i.radius_inner, 8) < outer:
                    # The Shell overlaps the outer edge of the current one
                    overlap_error = True
                    break"""
        
        assert not overlap_error, "Some Shells overlap each other"
        assert not gap_error, "Some Shells have gaps betweem them"
        

        # Is the outermost layer a solid and set to external?

        assert type(self._shells[-1].material) == Solid, "Outermost shell must be a Solid"

        for shell in self._shells:
            shell.external = False

        self._shells[-1].external = True

        self.verified = True
    
    def build_thermal_resistances(self, shell_temperatures, shell_pressures, g):
        """
        Outputs the thermal resistances of each shell layer
        
        Args:
            shell_temperatures (list of floats):    list of len(self._shells)+1, referring to temperatures at the boundary
                                                    of each material. 0th is the inner wall temperature of innermost shell,
                                                    1st is outer wall of innermost/inner wall of second.
            shell_pressures (list of floats):       list of len(self._shells), referring to pressures in each shell
            g (float):                              gravity
        """
        assert len(shell_temperatures) == len(self._shells)+1, "length of 'shell_temperatures' must be equal to number of Shells plus one"
        assert len(shell_pressures) == len(self._shells), "length of 'shell_pressures' must be equal to number of Shells"

        shell_thermal_resistances = np.zeros((len(self._shells)))

        for shell_count, shell in enumerate(self._shells):
            T_avg = (shell_temperatures[shell_count] + shell_temperatures[shell_count+1]) / 2
            T_delta = shell_temperatures[shell_count+1] - shell_temperatures[shell_count]
            pressure = shell_pressures[shell_count]

            shell_thermal_resistances[shell_count] = shell.thermal_resistance(T_avg, T_delta, p, g)
        
        return(shell_temperatures)








if __name__ == "__main__":
    steel = Solid("Steel", 150, 8700, 500)
    co2 = ConstrainedIdealGas("STP CO2", 44, 0.71, 10.9e-6, 749, 0.0153, p=101325)

    steel_a = StaticShell("vertical", steel, 1, 0.1, 5)
    steel_b = StaticShell("vertical", steel, 1.1, 0.1, 5)
    steel_c = StaticShell("vertical", steel, 1.2, 0.1, 5)
    steel_d = StaticShell("vertical", steel, 1.3, 0.1, 5)


    starship = Habitat("vertical", 70)
    starship.create_static_shell(co2, 1)

    #starship.create_static_shell(steel, 0.001)
    starship.append_shell(steel_a)
    starship.append_shell(steel_c)
    starship.append_shell(steel_b)

    print(starship._shells)
    starship.verify_geometry()
    print(starship._shells)


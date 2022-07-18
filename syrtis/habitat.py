"""
Top-level object and tools for the whole Habitat, composed of individual Shells
"""
if __name__ == "__main__":
    from shell import *
    from material import *
else:
    from syrtis.shell import *
    from syrtis.material import *

class Habitat:
    """
    Object to store the whole habitat geometry

    Args:
        configuration (str):        either 'horizontal' or 'vertical' to specify orientation of axis
        length (float):             the length of the central axis of the Habitat (m)    
    """

    def __init__(self, configuration, length):
        assert configuration == "horizontal" or configuration == "vertical", "'configuration' must be either 'horizontal' or 'vertical'"

        assert isinstance(length, Number), "Habitat 'length' must be a numerical value"
        assert length > 0, "Habitat 'length' must be a positive value"

        self.configuration = configuration
        self.length = length

        self._shells = []

        self.radius_outer = 0
    
    def create_static_shell(self, material, thickness, thermal_resistance=0):
        """
        Create a StaticShell that conforms around the outside of the outermost one
        """
        self._shells[-1].external = False
        self.radius_outer += thickness

        new_shell = StaticShell(self.configuration, material, self.radius_outer, thickness, self.length, True, thermal_resistance)

        self._shells.append(new_shell)
    
    def append_shell(self, shell):
        """
        Append a Shell of any type to the Habitat, checking for no overlaps
        """

        assert shell.radius_inner >= self.radius_outer, "Shell intersects with existing Shells"

        self._shells[-1].external = False

        if shell.radius_inner > self.radius_outer:
            print("Gap exists between inserted Shell and existing Shells. Filling gap with ambient Martian air")

            ambient_gap = StaticShell(self.configuration, ambient_atmosphere, self.radius_outer, 
            shell.radius_inner-self.radius_outer, self.length)
            self._shells.append(ambient_gap)
        
        self._shells.append(shell)

        self.radius_outer = shell.radius_outer
    
    def build_thermal_resistance(self, inner_temp):
        pass


steel = Solid(150, 8700, 500)
co2 = ConstrainedGas(210, 580, 0.71, 10.9e-6, 749, 8.74e-3, 0.0143)

starship = Habitat("vertical", 70)
starship.create_static_shell(co2, 4.5)
starship.create_static_shell(steel, 0.001)


"""
Classes and inherited subclasses for Shells, the general term for layers of material included in the simulation.
In most cases, a Shell will be added as a subclass with more specialised properties

"""



if __name__ == "__main__":
    from material import *
else:
    from syrtis.material import *


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
        material (syrtis.Material): the material of the Shell
        radius_inner (float):       the internal radius of the Shell (m)
        thickness (float):          the thickness of the Shell (m)
        external (bool) :           if True, the outer surface of the Shell is taken as the outer surface of the whole system. Default False
        thermal_resistance (float): an override for simple calculation of the Shell's thermal resistance. Ignored when not a positive number
    """

    def __init__(self, configuration, material, radius_inner, thickness, length, external=False, thermal_resistance=0.0):

        assert configuration == "horizontal" or configuration == "vertical", "'configuration' must be either 'horizontal' or 'vertical'"

        assert type(material) in material_classes, "'material' input is not valid Material-inherited object"

        assert isinstance(radius_inner, Number), "Shell 'internal_radius' must be a numerical value"
        assert radius_inner > 0, "Shell 'internal_radius' must be a positive value"

        assert isinstance(thickness, Number), "Shell 'thickness' must be a numerical value"
        assert thickness > 0, "Shell 'thickness' must be a positive value"

        self.material = material
        self.radius_inner = radius_inner
        self.thickness = thickness
        self.radius_outer = self.radius_inner + self.thickness


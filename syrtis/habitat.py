"""
Top-level object and tools for the whole Habitat, composed of individual Shells

References:

 - [1] - Y Cengel, Heat Transfer 2nd edition
 - [2] - Sky temperature modelisation and applications in building simulation, Adelard et al 1998
 - [3] - Long-wave  radiation from clear skies, Swinbank 1963
 - [4] - Radiation View Factors, http://webserver.dmt.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf
 - [5] - Thermal control of MSL Rover "Curiosity" using an Active Fluid Loop, Birur 2013
 - [6] - The thermal control system of NASA’s Curiosity rover: a case study, Quattrocchi et al, 2022 
 - [7] - Conduction heat transfer solutions, VanSant 1980
 - [8] - Thermal Spreading Resistance of Arbitrary-Shape Heat Sources on a Half-Space: A Unified Approach, Sadeghi et al 2010
 - [9] - Heat transfer from partially buried pipes, Morud & Simonsen 2007
 - [10]- Proposed OHTC Formula for Subsea Pipelines Considering Thermal Conductivities of Multi-Layered Soils, Park et al 2018
 - [11]- Heat Transfer, J P Holman 10th Edition
"""

import numpy as np
from random import choice
import matplotlib.patches as patches


if __name__ == "__main__":
    from shell import *
    from material import *
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

        assert is_numeric(length, not_negative=True), "Habitat 'length' must be a positive numerical value"

        assert endcap_type in endcap_types, "'endcap_type' must be a valid keyword"

        self.orientation = orientation
        self.length = length

        self.endcap_type = endcap_type

        self._shells = []
        self.groundlevel = None
        self.earthworks = None

        self.radius_outer = 0
        self.length_outer = 0
        self.ground_contact_angle = 0

        self.verified = False
    
    def draw(self,size=(20,20)):
        """
        Sets up a matplotlib object to draw the Habitat
        """
        self.verify_geometry()

        fig = plt.figure(figsize=size)

        if self.orientation == "horizontal":
            max_x = self.radius_outer * 1.25
            max_y = self.radius_outer * 1.25
            if self.endcap_type == "hemisphere":
                max_z = (self.length / 2) + (self.radius_outer * 1.25)
            else:
                max_z = (self.length / 2) * 1.25
        elif self.orientation == "vertical":
            max_x = self.radius_outer * 1.25
            if self.endcap_type == "hemisphere":
                max_y = (self.length / 2) + (self.radius_outer * 1.25)
            else:
                max_y = (self.length / 2) * 1.25
            max_z = self.radius_outer

        spec = fig.add_gridspec(2, 2)

        side_view = fig.add_subplot(spec[0, 0], aspect="equal")
        top_view = fig.add_subplot(spec[0, 1], aspect="equal")
        section = fig.add_subplot(spec[1,:],aspect="equal")
        #section = fig.add_subplot(spec[1,:], projection='polar')

        
        top_y = max_y
        if self.groundlevel != None and self.groundlevel.habitat_axis_height < -max_y:
            top_y = -self.groundlevel.habitat_axis_height+0.5
        if self.earthworks != None:
            top_y = self.earthworks.radius_inner - self.earthworks.axis_height_from_ground + self.earthworks.depth_of_axis + 2

        side_view_max_x = max(top_y, max_x)

        side_view.set_xlim(-side_view_max_x, side_view_max_x)
        side_view.set_ylim(-max_y, top_y)
        side_view.title.set_text("Elevation section view")

        top_view.set_xlim(-max_x, max_x)
        top_view.set_ylim(-max_z, max_z)
        top_view.title.set_text("Plan section view")

        section.set_xlim(-0.2 * max_x, 0.2 * max_x)
        section.set_ylim(0.95*self._shells[0].radius_outer, max_x*0.82)
        section.title.set_text("Wall cross-section detail")

        # Martian atmosphere
        side_view.fill_between(x=[-side_view_max_x, side_view_max_x], 
            y1=[top_y, top_y],
            y2=[-max_y, -max_y],
            color="#FFBEA0")

        # Ground level
        if self.groundlevel != None:
            if self.groundlevel.habitat_axis_height == 1e3:
                side_view.fill_between(x=[-side_view_max_x, side_view_max_x], 
                y1=[-0.9*max_y, -0.9*max_y],
                y2=[-max_y, -max_y],
                color="#BE4628")
            elif self.earthworks != None:
                side_view.fill_between(x=[-side_view_max_x, side_view_max_x], 
                y1=[top_y-0.5, top_y-0.5],
                y2=[-max_y, -max_y],
                color="#BE4628")

            else:
                if self.orientation == "horizontal":
                    side_view.fill_between(x=[-side_view_max_x, side_view_max_x], 
                    y1=[-self.groundlevel.habitat_axis_height, -self.groundlevel.habitat_axis_height],
                    y2=[-max_y, -max_y],
                    color="#BE4628")
                elif self.orientation == "vertical":
                    side_view.fill_between(x=[-side_view_max_x, side_view_max_x], 
                    y1=[-self.length/2 - self.groundlevel.habitat_axis_height, -self.length/2 - self.groundlevel.habitat_axis_height],
                    y2=[-max_y, -max_y],
                    color="#BE4628")
        elif self.earthworks != None:
                side_view.fill_between(x=[-side_view_max_x, side_view_max_x], 
                y1=[top_y-2, top_y-2],
                y2=[-max_y, -max_y],
                color="#BE4628")

        top_view.set_facecolor("#BE4628")

        # Earthworks if present
        if self.earthworks != None:
            earthworks_axis = self.earthworks.radius_inner - self.earthworks.axis_height_from_ground
            earthworks_circ = plt.Circle([0, earthworks_axis], self.earthworks.radius_inner, color="#FFBEA0")
            side_view.add_patch(earthworks_circ)
        

        for shell_count, shell in enumerate(reversed(self._shells)):
            col = choice(["r", "b", "g", "c", "y"])
            if "steel" in shell.material.name:
                col = "#777777"
            elif "aluminium" in shell.material.name:
                col = "#999999"
            elif "STP Air" in shell.material.name:
                col = "#87CEEB"
            elif "plastic" in shell.material.name:
                col = "b"
            elif "foam" in shell.material.name:
                col = "#333333"
            elif "CO2" in shell.material.name:
                col = "#FFA073"
            
            # Small section
            circ = plt.Circle([0, 0], shell.radius_outer, color=col)

            # Side section
            rect = patches.Rectangle([-shell.radius_outer, -shell.length/2], 2*shell.radius_outer, shell.length, color=col)

            # Endcaps
            if self.endcap_type == "hemisphere":
                endcap_1 = plt.Circle([0, -shell.length/2], shell.radius_outer, color=col)
                endcap_2 = plt.Circle([0, shell.length/2], shell.radius_outer, color=col)

            elif self.endcap_type == "flat" and shell.radius_inner != 0:
                endcap_thickness = sum(shell.thickness for shell in self._shells[1:len(self._shells) - shell_count])
                endcap_1 = patches.Polygon([
                    [-shell.radius_outer, -shell.length/2], [-0.2*shell.radius_outer, -shell.length/2 - endcap_thickness], 
                    [0.2*shell.radius_outer, -shell.length/2 - endcap_thickness], [shell.radius_outer, -shell.length/2]],
                    color=col)

                endcap_2 = patches.Polygon([
                    [-shell.radius_outer, shell.length/2], [-0.2*shell.radius_outer, shell.length/2 + endcap_thickness], 
                    [0.2*shell.radius_outer, shell.length/2 + endcap_thickness], [shell.radius_outer, shell.length/2]],
                    color=col)

            if self.orientation == "horizontal":
                side_view.add_patch(circ)
                
                top_view.add_patch(endcap_1)
                top_view.add_patch(endcap_2)
                top_view.add_patch(rect)
            
            elif self.orientation == "vertical":
                top_view.add_patch(circ)

                side_view.add_patch(endcap_1)
                side_view.add_patch(endcap_2)
                side_view.add_patch(rect)

            # Wedge
            wedge = plt.Circle([0, 0], shell.radius_outer, color=col)
            section.add_patch(wedge)
        
    def create_static_shell(self, material, thickness, thermal_resistance=0, parallel_thermal_resistance=0):
        """
        Create a StaticShell that conforms around the outside of the outermost one
        """
        
        new_shell = StaticShell(self.orientation, material, self.radius_outer, 
        thickness, self.length, False, thermal_resistance, parallel_thermal_resistance)

        self.radius_outer += thickness
        self.verified = False
        self._shells.append(new_shell)
    
    def create_ground_level(self, habitat_axis_height=1e3, thermal_resistance=0):
        """
        Create a GroundLevel object

        Args:
            habitat_axis_height (float):        height of the habitat central axis above ground level, for horizontal orientation
                                                height of the lower end of the cylinder above ground level, for vertical orientation
                                                Defaults to a large value, indicating no thermal conduction contact (m)
            thermal_resistance (float):         thermal resistance between Habitat outer wall, default to 0 = no thermal contact (K/W)
        """
        assert self.groundlevel==None, "A GroundLevel has already been set"

        groundlevel = GroundLevel(habitat_axis_height, thermal_resistance)

        self.groundlevel = groundlevel

    def create_earthworks(self, radius_inner, depth_of_axis, axis_height_from_ground=None):
        """
        Create an Earthworks object for the habitat, to represent a tunnel or similar

        Args:
            radius_inner (float):               internal radius of the earthworks cylinder (m)
            depth_of_axis (float):              depth of the central axis of the earthworks cylinder below ground level (m)
            axis_height_from_ground (float):    height of the habitat axis from the ground (m). 
                                                Overwritten by a GroundLevel object, if it exists
        """
        if self.groundlevel != None:
            axis_height_from_ground = self.groundlevel.habitat_axis_height
        
        habitat_outer_radius = max([shell.radius_outer for shell in self._shells])

        if axis_height_from_ground + habitat_outer_radius > (2 * radius_inner):
            return(ValueError("Habitat collides with inside of Earthworks. Make the Earthworks bigger or lower the habitat"))
        
        if depth_of_axis < radius_inner:
            return(ValueError("Syrtis does not currently support Earthworks which breach the surface"))
        
        self.earthworks = Earthworks(radius_inner, depth_of_axis, axis_height_from_ground)

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

        # Set whole-Habitat radius_outer and length_outer
        self.radius_outer = self._shells[-1].radius_outer
        self.length_outer = self._shells[-1].length


        # If ground level is set, use it to calculate angle of contact
        if self.groundlevel != None and self.orientation == "horizontal":
            if self.groundlevel.habitat_axis_height > self.radius_outer:
                # Habitat is entirely above the ground
                self.ground_contact_angle = 0

            elif self.groundlevel.habitat_axis_height < (-self.radius_outer):
                # Habitat is entirely below ground
                self.ground_contact_angle = 180
            
            else:
                self.ground_contact_angle = np.rad2deg(np.arccos(
                    (self.radius_outer - self.groundlevel.habitat_axis_height) / self.radius_outer))

        
        self.exposed_convective_area()

        self.verified = True
    
    def exposed_convective_area(self):
        """
        Find the area exposed to convective losses, endcap and cylinder.
        Called during geometry verification
        """

        if self.groundlevel == None:
            self.exposed_area_cylinder = 2 * np.pi * self.radius_outer * self.length_outer

            if self.endcap_type == "hemisphere":
                self.exposed_area_endcap = 4 * np.pi * np.power(self.radius_outer, 2)
            
            elif self.endcap_type == "flat":
                self.exposed_area_endcap = 2 * np.pi * np.power(self.radius_outer, 2)
        
        elif self.groundlevel.habitat_axis_height > self.radius_outer:
            self.exposed_area_cylinder = 2 * np.pi * self.radius_outer * self.length_outer

            if self.endcap_type == "hemisphere":
                self.exposed_area_endcap = 4 * np.pi * np.power(self.radius_outer, 2)
            
            elif self.endcap_type == "flat":
                self.exposed_area_endcap = 2 * np.pi * np.power(self.radius_outer, 2)
        
        elif (self.groundlevel.habitat_axis_height < self.radius_outer and self.groundlevel.habitat_axis_height > -self.radius_outer and 
        self.orientation == "horizontal"):
            # Habitat is partially below ground level, so exposed area of cylinder and endcaps are reduced

            self.exposed_area_cylinder = 2 * np.pi * self.radius_outer * self.length_outer * 0.5 * (
                1 + np.cos(np.deg2rad(self.ground_contact_angle)))
            
            if self.endcap_type == "hemisphere":
                # Spherical cap area formula
                self.exposed_area_endcap = 2 * np.pi * self.radius_outer * (
                    self.radius_outer + self.groundlevel.habitat_axis_height) 
            
            elif self.endcap_type == "flat":
                # Segment formula area formula
                self.exposed_area_endcap = 2 * np.power(self.radius_outer, 2) * (np.pi - 
                    np.deg2rad(self.ground_contact_angle) + 0.5 * np.cos(np.deg2rad(2 * self.ground_contact_angle)))
        
        elif self.groundlevel.habitat_axis_height < -self.radius_outer and self.orientation == "horizontal":
            # The entire habitat is below ground level
            self.exposed_area_cylinder = 0
            self.exposed_area_endcap = 0
            
        elif (self.groundlevel.habitat_axis_height < self.radius_outer and self.groundlevel.habitat_axis_height > 0 and 
        self.orientation == "vertical"):
            # The cylinder of the habitat is above ground level but a hemsipherical endcap would be truncated
                self.exposed_area_cylinder = 2 * np.pi * self.radius_outer * self.length_outer

                if self.endcap_type == "hemisphere":
                    # One truncated cap, one full cap
                    # Spherical cap area formulae
                    self.exposed_area_endcap = 2 * np.pi * self.radius_outer * self.groundlevel.habitat_axis_height
                    self.exposed_area_endcap += 2 * np.pi * np.power(self.radius_outer, 2)
                
                elif self.endcap_type == "flat":
                    self.exposed_area_endcap = 2 * np.pi * np.power(self.radius_outer, 2)
        
        elif (self.groundlevel.habitat_axis_height < 0 and self.groundlevel.habitat_axis_height > -self.length_outer and
            self.orientation == "vertical"):
            # A section of the cylinder is below ground
            self.exposed_area_cylinder = 2 * np.pi * self.radius_outer * (
                    self.length_outer + self.groundlevel.habitat_axis_height)
                
            if self.endcap_type == "hemisphere":
                self.exposed_area_endcap = 2 * np.pi * np.power(self.radius_outer, 2)

            elif self.endcap_type == "flat":
                self.exposed_area_endcap = np.pi * np.power(self.radius_outer, 2)
        
        elif (self.groundlevel.habitat_axis_height < -self.length_outer and self.groundlevel.habitat_axis_height < -(self.length_outer + self.radius_outer) and
            self.orientation == "vertical"):
            # The whole cylinder is below ground, but a hemispherical endcap may protrude
            self.exposed_area_cylinder = 0

            if self.endcap_type == "hemisphere":
                self.exposed_area_endcap = 2 * np.pi * self.radius_outer * (
                    self.length_outer + self.radius_outer + self.groundlevel.habitat_axis_height)
            
            elif self.endcap_type == "flat":
                self.exposed_area_endcap = 0
        
        elif self.groundlevel.habitat_axis_height < -(self.length_outer + self.radius_outer) and self.orientation == "vertical":
            # Habitat is entirely underground
            self.exposed_area_cylinder = 0
            self.exposed_area_endcap = 0

    def direct_solar_area(self, solar_altitude, solar_azimuth):
        """
        Find the area of the habitat in direct sunlight
        Called by solar_gain_direct to find areas
        
        Args:
            solar_altitude (float): vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):  horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
        """
        solar_altitude_rad = np.deg2rad(solar_altitude)
        solar_azimuth_rad = np.deg2rad(solar_azimuth)

        direct_solar_area = 0

        if self.groundlevel == None:
            if self.orientation == "horizontal":
                direct_solar_area += 2 * self.length_outer * self.radius_outer * (
                np.cos(solar_altitude_rad) * np.cos(solar_azimuth_rad) + np.sin(solar_azimuth_rad))
                
                if self.endcap_type == "flat":
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                        np.cos(solar_altitude_rad) * abs(np.cos(solar_azimuth_rad)))
                
                elif self.endcap_type == "hemisphere":
                    # Area of hemisphere projected onto the plane perpendicular to the Sun
                    # Two components are axial (plan view) and radial-ish (side view)
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                        np.cos(solar_altitude_rad) + 0.5 * np.sin(solar_altitude_rad)) * abs(np.cos(solar_azimuth_rad))

            elif self.orientation == "vertical":
                direct_solar_area += 2 * self.length_outer * self.radius_outer * (
                np.cos(solar_altitude_rad)) 

                if self.endcap_type == "flat":
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * np.sin(solar_altitude_rad) 
                
                elif self.endcap_type == "hemisphere":
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                        np.sin(solar_altitude_rad) + 0.5 * np.cos(solar_altitude_rad))
        
        elif self.groundlevel.habitat_axis_height > self.radius_outer:
            if self.orientation == "horizontal":
                direct_solar_area += 2 * self.length_outer * self.radius_outer * (
                np.cos(solar_altitude_rad) * np.cos(solar_azimuth_rad) + np.sin(solar_azimuth_rad))
                
                if self.endcap_type == "flat":
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                        np.cos(solar_altitude_rad) * abs(np.cos(solar_azimuth_rad)))
                
                elif self.endcap_type == "hemisphere":
                    # Area of hemisphere projected onto the plane perpendicular to the Sun
                    # Two components are axial (plan view) and radial-ish (side view)
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                        np.cos(solar_altitude_rad) + 0.5 * np.sin(solar_altitude_rad)) * abs(np.cos(solar_azimuth_rad))

            elif self.orientation == "vertical":
                direct_solar_area += 2 * self.length_outer * self.radius_outer * (
                np.cos(solar_altitude_rad)) 

                if self.endcap_type == "flat":
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * np.sin(solar_altitude_rad) 
                
                elif self.endcap_type == "hemisphere":
                    direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                        np.sin(solar_altitude_rad) + 0.5 * np.cos(solar_altitude_rad))
        
        elif (self.groundlevel.habitat_axis_height < self.radius_outer and self.groundlevel.habitat_axis_height > -self.radius_outer and 
        self.orientation == "horizontal"):
            # Habitat is partially below ground level, so exposed area of cylinder and endcaps are reduced
            direct_solar_area += 2 * self.length_outer * self.radius_outer * (
                np.cos(solar_altitude_rad) * np.cos(solar_azimuth_rad) * 0.5 * (1 + np.cos(np.deg2rad(self.ground_contact_angle)))
                + np.sin(solar_azimuth_rad)) 
            
            if self.endcap_type == "flat":
                direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                    np.cos(solar_altitude_rad) * abs(np.cos(solar_azimuth_rad))) * 0.5 * (
                1 + np.cos(np.deg2rad(self.ground_contact_angle)))
            
            elif self.endcap_type == "hemisphere":
                # Area of hemisphere projected onto the plane perpendicular to the Sun
                # Two components are axial (plan view) and radial-ish (side view)
                direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                    np.cos(solar_altitude_rad) + 0.5 * np.sin(solar_altitude_rad)) * abs(np.cos(solar_azimuth_rad)) * 0.5 * (
                1 + np.cos(np.deg2rad(self.ground_contact_angle)))
        
        elif self.groundlevel.habitat_axis_height < -self.radius_outer and self.orientation == "horizontal":
            # Habitat is entirely below ground level
            pass

        elif (self.groundlevel.habitat_axis_height < 0 and self.groundlevel.habitat_axis_height > -self.length_outer and
            self.orientation == "vertical"):
            direct_solar_area += 2 * (self.length_outer + self.groundlevel.habitat_axis_height) * self.radius_outer * (
                np.cos(solar_altitude_rad)) 

            if self.endcap_type == "flat":
                direct_solar_area += np.pi * np.power(self.radius_outer, 2) * np.sin(solar_altitude_rad) 
            
            elif self.endcap_type == "hemisphere":
                direct_solar_area += np.pi * np.power(self.radius_outer, 2) * (
                    np.sin(solar_altitude_rad) + 0.5 * np.cos(solar_altitude_rad))
        
        elif (self.groundlevel.habitat_axis_height < -self.length_outer and self.groundlevel.habitat_axis_height < -(self.length_outer + self.radius_outer) and
            self.orientation == "vertical"):
            # The whole cylinder is below ground, but a hemispherical endcap may protrude
            if self.endcap_type == "hemisphere":
                direct_solar_area += np.pi * np.power(self.radius_outer * 
                np.cos(np.arcsin((self.length_outer + self.groundlevel.habitat_axis_height)/self.radius_outer)), 2) * (
                    np.sin(solar_altitude_rad) + 0.5 * np.cos(solar_altitude_rad))
            
            elif self.endcap_type == "flat":
                pass
        
        elif self.groundlevel.habitat_axis_height < -(self.length_outer + self.radius_outer) and self.orientation == "vertical":
            pass
    
        return(direct_solar_area)
        
    def indirect_solar_area(self):
        """
        Find the area of the habitat indirect sunlight
        Called by solar_gain_direct and solar_gain_indirect to find areas

        Very similar to convective_area, except a vertically oriented Habitat with flat endcaps will not have the flat endcap counted
        """
        indirect_solar_area = self.exposed_area_cylinder + self.exposed_area_endcap

            
        if self.orientation == "horizontal":
            pass
        elif self.groundlevel == None:
            pass
        elif (self.groundlevel.habitat_axis_height < -self.length_outer and self.groundlevel.habitat_axis_height < -(self.length_outer + self.radius_outer) and
            self.orientation == "vertical"):
            # One flat endcap is above ground level but isn't visible to the surface of the ground: exclude it
            indirect_solar_area -= np.pi * np.power(self.radius_outer, 2)
        
        return(indirect_solar_area)

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
                
                Nu_D = 0.228 * np.power(Re, 0.731) * np.power(Pr, 1/3)
            
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
        Calculate the view factor to the habitat to flat ground
        Uses equations from Reference [4]
        """

        if self.orientation == "vertical" and self.earthworks == None:
            vf_cylinder = 0.5

            vf_endcap = 0.5

        elif self.orientation == "horizontal" and self.earthworks == None:
            vf_cylinder = 0.5 

            vf_endcap = 0.5
        
        elif self.earthworks != None:

            vf_cylinder = self.earthworks.view_factor_ground_cylinder(
                self.orientation, self.radius_outer, self.length_outer)
            
            if self.endcap_type == "hemisphere":
                vf_endcap = self.earthworks.view_factor_ground_hemisphere(
                    self.orientation, self.radius_outer, self.length_outer)
            
            elif self.endcap_type == "flat":
                vf_endcap = self.earthworks.view_factor_ground_disc(
                    self.orientation, self.radius_outer, self.length_outer)
            
        vf = ((vf_cylinder * self.exposed_area_cylinder) + 
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
        return(-Q_sky)
    
    def radiative_loss_ground(self, T_wall):
        """
        Calculate radiative loss to the sky.

        Args:
            T_wall (float):             temperature of the wall (K)
        """

        vf_ground = self.view_factor_ground()

        Q_ground = 5.67e-8 * self._shells[-1].material.emit * np.power(T_wall, 4) * (
            self.exposed_area_cylinder + self.exposed_area_endcap) * vf_ground

        return(Q_ground)
    
    def radiative_gain_ground(self, T_ground):
        """
        Calculate radiative loss to the sky.

        Args:
            T_ground (float):              temperature of the ground(K)
        """

        vf_ground = self.view_factor_ground()

        Q_ground_gain = 5.67e-8 * self._shells[-1].material.absorb * np.power(T_ground, 4) * (
            self.exposed_area_cylinder + self.exposed_area_endcap) * vf_ground

        # A negative sign is used for consistency with convention that +ve Q = heat loss
        return(-Q_ground_gain)

    def solar_gain_direct(self, solar_altitude, solar_azimuth, solar_intensity):
        """
        Calculate heat gain to the outermost layer of the habitat, or into the centre if all outer layers have sufficient transparency
        
        Args:
            solar_altitude (float): vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):  horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
            solar_intensity (float):power delivered by solar radiation after dust absorption, W/m2
        """

        if self.earthworks == None:
            solar_intensity_multipler = 1
        else:
            if self.orientation == "vertical":
                axis_height_from_ground = 0
            elif self.orientation == "horizontal":
                axis_height_from_ground = self.radius_outer

            if self.groundlevel != None:
                axis_height_from_ground += self.groundlevel.habitat_axis_height
            
            solar_intensity_multipler = self.earthworks.direct_solar_intensity_scaling(
                self.orientation, self.radius_outer, self.length_outer, axis_height_from_ground, 
                solar_altitude, solar_azimuth)

        direct_lit_area = self.direct_solar_area(solar_altitude, solar_azimuth)
        
        if self._shells[-1].material.transmit == 0:
            # Outer layer is opaque - all solar energy is delivered to outermost shell
            Q_solar_direct = solar_intensity * solar_intensity_multipler * direct_lit_area * self._shells[-1].material.absorb
        
        # A negative sign is used for consistency with convention that +ve Q = heat loss
        return(-abs(Q_solar_direct))
    
    def solar_gain_indirect(self, solar_altitude, solar_azimuth, solar_intensity, albedo_ground):
        """
        Calculate heat gain to the outermost layer of the habitat, or into the centre if all outer layers have sufficient transparency
        
        Args:
            solar_altitude (float): vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):  horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
            solar_intensity (float):power delivered by solar radiation after dust absorption, W/m2
            albedo_ground (float):  albedo of the surface around the habitat
        """

        if self.earthworks == None:
            solar_intensity_multipler = 1
        else:
            if self.orientation == "vertical":
                axis_height_from_ground = 0
            elif self.orientation == "horizontal":
                axis_height_from_ground = self.radius_outer

            if self.groundlevel != None:
                axis_height_from_ground += self.groundlevel.habitat_axis_height
            
            solar_intensity_multipler = self.earthworks.indirect_solar_intensity_scaling(
                self.orientation, self.radius_outer, self.length_outer, axis_height_from_ground, 
                solar_altitude, solar_azimuth)

        indirect_lit_area = self.indirect_solar_area()
        
        if self._shells[-1].material.transmit == 0:
            # Outer layer is opaque - all solar energy is delivered to outermost shell
            Q_solar_indirect = solar_intensity * solar_intensity_multipler * albedo_ground * indirect_lit_area * self._shells[-1].material.absorb
        
        # A negative sign is used for consistency with convention that +ve Q = heat loss
        return(-Q_solar_indirect)

    def conductive_loss_fixed_resistance(self, T_wall, T_ground, R_ground):
        """
        Calculate thermal resistance losses for a fixed thermal resistance value

        Args:
            T_wall (float):     temperature of the outer wall (K)
            T_ground (float):   temperature of the far-field ground (K)
            R_ground (float):   thermal resistance from the outer wall to the ground (K/W)
        """

        Q_ground = (T_wall - T_ground) / R_ground
        
        return(Q_ground)

    def conductive_loss_horizontal_cylinder_steady(self, T_wall, T_ground, k_ground, R_wall):
        """
        Calculate the conductive loss from a horizontal cylinder in contact with the ground, steady-state solution
        
        Uses formulae from [1] for fully buried.
        Uses formulae from [9] and [10] for partially buried

        Args:
            T_wall (float):         temperature of the outer wall (K)
            T_ground (float):       temperature of the far-field ground (K)
            k_ground (float):       thermal conductivity of the far-field ground (W/m/K)
            R_wall (float):         thermal resistance of the entire wall, required for Biot number (K/W)
        """

        if self.groundlevel == None:
            return(0)
        if self.groundlevel.habitat_axis_height > self.radius_outer:
            # Habitat is entirely above the ground
            return (0)
        
        elif self.groundlevel.habitat_axis_height < (-self.radius_outer):
            # Habitat is entirely below the ground - use formula from [1] (and [9] if shallow buried)
            buried_depth = -self.groundlevel.habitat_axis_height
            depth_radius = buried_depth / self.radius_outer
            
            if depth_radius > 3:
                # Formula (1) from Table 3.5 in [1]
                S = 2 * np.pi * self.length_outer / np.log(4 * buried_depth / self.radius_outer)

                Q_ground = k_ground * S * (T_wall - T_ground)

            else:
                # Formulae 2-5 from [9]
                U_wall = 1 / (R_wall * 2 * np.pi * self.length_outer * self.radius_outer)
                Bi = U_wall * self.radius_outer / k_ground
                alpha_0 = np.log(depth_radius + np.sqrt(np.power(depth_radius, 2) - 1))

                U = Bi * (k_ground / self.radius_outer) / np.sqrt(
                    1 + np.power(Bi * alpha_0, 2) + (2 * Bi * alpha_0 / np.tanh(alpha_0)))
                
                Q_ground = U * (2 * np.pi * self.length_outer * self.radius_outer) * (T_wall - T_ground)

        
        else:
            # Habitat is partially buried
            # Use formulae from [9] and [10]

            buried_depth = -self.groundlevel.habitat_axis_height
            depth_radius = buried_depth / self.radius_outer

            U_wall = 1 / (R_wall * 2 * np.pi * self.length_outer * self.radius_outer)

            Bi = U_wall * self.radius_outer / k_ground

            theta_b = np.arccos(depth_radius)

            C1 = np.sqrt(1 - np.power(depth_radius, 2))
            C2 = depth_radius + (C1 / (Bi * theta_b))

            if C2 > 1:
                Nu = 2 / (theta_b * (np.pi - theta_b)) * C1 / np.sqrt(np.power(C2, 2) - 1) * (
                np.pi/2 - np.arctan(np.sqrt((C2 + 1) / (C2 - 1)) * np.tan(theta_b / 2)))
            else:# C2 <= 1:
                Nu = 1 / (theta_b * (np.pi - theta_b)) * C1 / np.sqrt(1 - np.power(C2, 2)) * np.log(
                (np.tan(theta_b / 2) + np.sqrt((1 - C2) / (1 + C2))) / 
                (np.tan(theta_b / 2) - np.sqrt((1 - C2) / (1 + C2))))
            
            U = Nu * k_ground / self.radius_outer

            Q_ground = U * (2 * np.pi * self.length_outer * self.radius_outer * (
                1 - (theta_b / np.pi))) * (T_wall - T_ground)
        
        return(Q_ground)

    def conductive_loss_vertical_cylinder_steady(self, T_wall, T_ground, k_ground):
        """
        Calculate the conductive loss from a vertical cylinder in contact with the ground, steady-state solution
        
        Uses formulae from [1] for partially buried, and for fully buried with a correction term

        Args:
            T_wall (float):         temperature of the outer wall (K)
            T_ground (float):       temperature of the far-field ground (K)
            k_ground (float):       thermal conductivity of the far-field ground (W/m/K)
        """

        if self.groundlevel == None:
            return(0)
        
        elif self.groundlevel.habitat_axis_height >= 0:
            return (0)

        elif self.groundlevel.habitat_axis_height < (-self.length_outer):
            # Habitat is partially below 
            buried_depth = -self.groundlevel.habitat_axis_height

            S = 2 * np.pi * buried_depth / np.log(2 * buried_depth / self.radius_outer)

            Q_ground = k_ground * S * (T_wall - T_ground)

        elif self.groundlevel.habitat_axis_height > (-self.length_outer):
            # Habitat is entirely below ground
            buried_depth = -self.groundlevel.habitat_axis_height

            S_spread = 2 * np.pi * buried_depth / np.log(2 * buried_depth / self.radius_outer)
            R_spread = 1 / (k_ground * S_spread)

            R_overburden = abs(buried_depth - self.length_outer) / (k_ground * self.radius_outer * 2)

            R = R_spread + R_overburden

            Q_ground = (T_wall - T_ground) / R

        
        return(Q_ground)

    def conductive_loss_hemisphere_steady(self, T_wall, T_ground, k_ground):
        """
        Calculate the conductive loss from a sphere, steady-state solution. Used for spherical endcaps
        
        Uses equation 3.1.7 from [7]

        Args:
            T_wall (float):         temperature of the outer wall (K)
            T_ground (float):       temperature of the far-field ground (K)
            k_ground (float):       thermal conductivity of the far-field ground (W/m/K)
        """
        if self.groundlevel == None:
            return(0)
        
        if self.orientation == "horizontal":
            if self.groundlevel.habitat_axis_height > self.radius_outer:
                return (0)
            else:
                buried_depth = -self.groundlevel.habitat_axis_height

                S = 4 * np.pi * self.radius_outer / (1 + 1 / (buried_depth/self.radius_outer + 1))

                Q_ground = k_ground * S * (T_wall - T_ground)
        
        elif self.orientation == "vertical":
            if self.groundlevel.habitat_axis_height > self.radius_outer:
                return(0)
            elif self.groundlevel.habitat_axis_height < -self.length_outer:
                # Bottom hemisphere only - assume the heat flux is half that of a sphere at equivalent depth
                buried_depth = -self.groundlevel.habitat_axis_height

                S = 2 * np.pi * self.radius_outer / (1 + 1 / (buried_depth/self.radius_outer + 1))

                Q_ground = k_ground * S * (T_wall - T_ground)

            elif self.groundlevel.habitat_axis_height > -self.length_outer:
                # Both hemisphere are partially buried
                # As a very janky approximation, assume this is equal to a sphere with centre at the mid-length
                equiv_buried_depth = -self.groundlevel.habitat_axis_height - (self.length_outer / 2)

                S = 4 * np.pi * self.radius_outer / (1 + 1 / (equiv_buried_depth/self.radius_outer + 1))

                Q_ground = k_ground * S * (T_wall - T_ground)


        return(Q_ground)

    def conductive_loss_disc_steady(self, T_wall, T_ground, k_ground):
        """
        Calculate the conductive loss from a disc in various orientations, steady-state solution. Used for flat endcaps
        
        Uses formula from [1]

        Args:
            T_wall (float):         temperature of the outer wall (K)
            T_ground (float):       temperature of the far-field ground (K)
            k_ground (float):       thermal conductivity of the far-field ground (W/m/K)
        """

        if self.groundlevel == None:
            return(0)
        
        elif self.orientation == "horizontal":
            if self.groundlevel.habitat_axis_height > self.radius_outer:
                # Habitat is entirely above the ground
                return (0)
            else:
                #print("Conduction loss from disc: flux patterns assumed to be the same as a sphere, with reduced area")
                buried_depth = -self.groundlevel.habitat_axis_height

                S = 2 * np.pi * self.radius_outer / (1 + 1 / (buried_depth/self.radius_outer + 1))

                Q_ground = k_ground * S * (T_wall - T_ground)
        
        elif self.orientation == "vertical":
            if self.groundlevel.habitat_axis_height > 0:
                # Habitat is entirely above the ground
                return (0)
            
            elif self.groundlevel.habitat_axis_height == 0:
                # Disc is resting on the surface
                S = 4 * self.radius_outer

                Q_ground = k_ground * S * (T_wall - T_ground)
            
            elif self.groundlevel.habitat_axis_height > -self.length:
                # Just the bottom disc 
                # Use Equation from Table 3-1 in [11], halved

                buried_depth = -self.groundlevel.habitat_axis_height

                S = 2 * np.pi * self.radius_outer / (
                    0.5 * np.pi - np.arctan(self.radius_outer / (2 * buried_depth)))
                
                Q_ground = k_ground * S * (T_wall - T_ground)
            
            else:
                # Both discs buried
                buried_depth_lower = -self.groundlevel.habitat_axis_height
                buried_depth_upper = -self.groundlevel.habitat_axis_height - self.length

                S_lower = 2 * np.pi * self.radius_outer / (
                    0.5 * np.pi - np.arctan(self.radius_outer / (2 * buried_depth_lower)))
                
                S_upper = 2 * np.pi * self.radius_outer / (
                    0.5 * np.pi - np.arctan(self.radius_outer / (2 * buried_depth_upper)))

                Q_ground = k_ground * (S_lower + S_upper) * (T_wall - T_ground)

        return(Q_ground) 

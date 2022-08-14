"""
Stores the Configuration for a given simulation: habitat heat loading and external environment
"""

from syrtis.material import *


solution_types = ["constant power", "constant temperature"]
air_directions = ["axial", "cross"]

class Configuration:
    def __init__(self, name, solution_type,
    T_ground, k_ground, albedo_ground, T_air, p_air, v_air, air_direction, solar_altitude, solar_azimuth, solar_intensity,
    Q_habitat=0, T_habitat=0):
        """
        An object to store the Configuration for a given simulation

        Args:
            name (str):             an identifier for the configuration
            solution_type(str):     either "constant power" or "constant temperature" for now, more may be added
            
            T_ground (float):       ground/soil temperature (K)
            k_ground (float):       soil thermal conductivity (W/m/K)
            albedo_ground (float)   albedo of the ground, between 0 and 1
            
            T_air (float):          air temperature (K)
            p_air (float):          air pressure (Pa)
            v_air (float):          speed of the air (m/s)
            air_direction (float):  direction of the air, either "axial" for along a horizontal cyclinder or "cross" for crossflow

            solar_altitude (float): vertical angle of the sun above the horizon (degrees)
            solar_azimuth (float):  horizontal angle of the sun RELATIVE TO HABITAT AXIS (degrees) - 0=directly along axis
            solar_intensity (float):power delivered by solar radiation after dust absorption, W/m2

            Q_habitat (float):      constant power generation in the Habitat (W) - optional, defaults to zero
            T_habitat (float):      constant temeprature in the Habitat (K) - optional unless solution_type is "constant temperature"  
        """

        assert solution_type in solution_types, "Configuration 'solution_type' must be a valid keyword"

        #assert solution_type == "constant temperature" and T_habitat != 0, "If 'constant temperature' has been selected, T_habitat must be set"

        assert air_direction in air_directions, "Configuration 'air_direction' must be a valid keyword"

        assert is_numeric(albedo_ground, unit=True), "Configuration 'albedo_ground' must be between 0 and 1"

        self.name = name
        self.solution_type = solution_type

        self.T_ground = T_ground
        self.k_ground = k_ground
        self.albedo_ground = albedo_ground

        self.T_air = T_air
        self.p_air = p_air
        self.v_air = v_air
        self.air_direction = air_direction

        self.solar_altitude = solar_altitude
        self.solar_azimuth = solar_azimuth
        self.solar_intensity = solar_intensity

        self.Q_habitat = Q_habitat
        self.T_habitat = T_habitat

        self.GRAVITY = 3.71
        # Assume we're only using Mars for now

        self.air = ConstrainedIdealGas("Surface air", self.p_air, 44, 0.71, 10.9e-6, 749, 0.0153)
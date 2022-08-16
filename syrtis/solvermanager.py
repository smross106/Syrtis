"""
Object for sweeping through different parameters in the Solver

References:
 - [1] - A solar azimuth formula that renders circumstantial treatment unnecessary without compromising mathematical rigor: 
        Mathematical setup, application and extension of a formula based on the subsolar point and atan2 function, Zhang et al 2021
 - [2] - Über die Extinktion des Lichtes in der Erdatmosphäre, Schoenberg 1929
 - [3] - Mars Fact Sheet, https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
 - [4] - Thermal control of MSL Rover "Curiosity" using an Active Fluid Loop, Birur 2013 
"""
import numpy as np
import copy

from syrtis.configuration import Configuration
from syrtis.solver import Solver
from itertools import product

SOL_HRS = 24.65
MARS_RADIUS = 3396.2e3
MARS_ATMOSPHERE_OPTICAL_HEIGHT = 8624
MARS_ECC = 0.093377


class SolverManager:
    """
    Parent class for all other solver managers
    """


    def __init__(self):
        pass

    def run_all_configurations(self, verbose=False):
        """
        Run all the configurations
        """

        heat_losses = []
        reports = []
        completed_solvers = []
        
        for configuration_index, configuration in enumerate(self.configurations):
            self.habitat.verify_geometry()

            solver = Solver(str(self.each_configuration_inputs_dicts[configuration_index]["name"]),
                self.habitat, configuration)
            
            solver.solve()
            heat_losses.append(solver.heat)
            reports.append(solver.report)
            completed_solvers.append(solver)
            
        self.heat_losses = heat_losses
        self.reports = reports
        self.completed_solvers = completed_solvers

        if verbose:
            return(self.each_configuration_inputs_dicts, self.heat_losses, self.reports)
        else:
            return(self.each_configuration_inputs_dicts, self.heat_losses)

class ConfigurationManager(SolverManager):
    """
    Sweep through different parameters in the Configuration
    """
    configuration_variables = ["name", "solution_type",
        "T_ground", "k_ground", "albedo_ground", "T_air", "p_air", "v_air", "air_direction", "solar_altitude", "solar_azimuth", "solar_intensity",
        "Q_habitat", "T_habitat"]

    def __init__(self, habitat, configuration, variable_parameters, all_configurations=[]):
        
        assert all(parameter in self.configuration_variables for parameter in list(variable_parameters.keys())), "Not all configuration parameters set"

        self.habitat = habitat
        self.configuration_baseline = configuration

        self.variable_parameters = variable_parameters

        if all_configurations == []:
            self.configurations = self.create_configurations()
        else:
            self.configurations = all_configurations

    def create_configurations(self):
          
        iterations_number = []
        inputs_keys = []
        iterations_values = []

        for key in self.configuration_variables:
            if not key in self.variable_parameters:
                pass
            elif type(self.variable_parameters[key]) == list or type(self.variable_parameters[key]) == tuple:
                iterations_number.append(len(self.variable_parameters[key]))
                inputs_keys.append(key)
                iterations_values.append(self.variable_parameters[key])
        
        self.all_configuration_inputs = np.zeros(iterations_number)

        each_configuration_inputs = list(product(*iterations_values))


        self.each_configuration_inputs_dicts = []
        
        for config_input in each_configuration_inputs:
            input_dict = copy.deepcopy(self.configuration_baseline.__dict__)

            for input_index, input in enumerate(config_input):
                input_dict[inputs_keys[input_index]] = input

            self.each_configuration_inputs_dicts.append(input_dict)

        all_configurations = []

        for input_dict in self.each_configuration_inputs_dicts:
            config = Configuration(
                name = input_dict["name"], 
                solution_type = input_dict["solution_type"],
                T_ground = input_dict["T_ground"], 
                k_ground = input_dict["k_ground"], 
                albedo_ground = input_dict["albedo_ground"], 
                T_air = input_dict["T_air"], 
                p_air = input_dict["p_air"], 
                v_air = input_dict["v_air"], 
                air_direction = input_dict["air_direction"], 
                solar_altitude = input_dict["solar_altitude"], 
                solar_azimuth = input_dict["solar_azimuth"],
                solar_intensity = input_dict["solar_intensity"],
                Q_habitat = input_dict["Q_habitat"], 
                T_habitat = input_dict["T_habitat"]
            )
            all_configurations.append(config)
        
        return(all_configurations)

class DayManager(SolverManager):
    """
    Calculate heat flux at each point in a day
    
    Args:
        habitat (syrtis.Habitat):               habitat under consideration
        configuration (syrtis.Configuration):   baseline configuration - temperatures and solar properties will be overwritten

    """
    def __init__(self, habitat, configuration, num_points,
        atmosphere_tau, latitude, areocentric_longitude,
        T_air_max, T_air_min, T_ground_max, T_ground_min, time_air_peak=12, time_ground_peak=12):

        self.habitat = habitat
        self.configuration = configuration
        self.num_points = num_points

        self.atmosphere_tau = atmosphere_tau
        self.latitude = latitude
        self.areocentric_longitude = areocentric_longitude

        self.T_air_max = T_air_max
        self.T_air_min = T_air_min
        self.T_ground_max = T_ground_max
        self.T_ground_min = T_ground_min

        self.time_air_peak = time_air_peak
        self.time_ground_peak = time_ground_peak

        self.create_configurations()
              
    def create_configurations(self):
        """
        Create all the configurations for the day
        """
        self.times = np.linspace(0, SOL_HRS, self.num_points)
        self.generate_solar_data()
        self.T_airs = self.generate_temperatures(self.T_air_max, self.T_air_min, self.time_air_peak)
        self.T_grounds = self.generate_temperatures(self.T_ground_max, self.T_ground_min, self.time_ground_peak)

        baseline_configuration = copy.deepcopy(self.configuration)

        configurations = []
        each_configuration_inputs_dicts = []
        
        for point, time in enumerate(self.times):

            time_configuration = copy.deepcopy(baseline_configuration)
            time_configuration.name = time_configuration.name + " {:.2f}".format(time)

            time_configuration.T_air = self.T_airs[point]
            time_configuration.T_ground = self.T_grounds[point]

            time_configuration.solar_intensity = self.solar_intensities[point]
            time_configuration.solar_altitude = self.solar_altitudes[point]
            time_configuration.solar_azimuth = self.solar_azimuths[point]

            configuration_dict = copy.deepcopy(time_configuration.__dict__)

            configurations.append(time_configuration)
            each_configuration_inputs_dicts.append(configuration_dict)
        
        self.configurations = configurations
        self.each_configuration_inputs_dicts = each_configuration_inputs_dicts
    
    def generate_solar_data(self):
        """
        Generate lists of solar intensities, altitudes and azimuths for a single day.

        Solar position implements method from Reference [1] with necessary modifications for Martian orbital parameters.

        Martian atmospheric optimal height, for use in calculation from Source [2]. Calculation is based on constant-density
        hydrostatic argument, with data taken from [3]
        """
        
        axis_obliquity = np.deg2rad(24.936)
        subsolar_latitude = np.arcsin(np.sin(axis_obliquity) * np.sin(np.deg2rad(self.areocentric_longitude)))

        solar_intensity_space = 590 * np.power(
        (1 + MARS_ECC * np.cos(np.deg2rad(self.areocentric_longitude - 248))) / (1 - MARS_ECC), 2)

        solar_intensities = np.zeros((self.num_points))
        solar_altitudes = np.zeros((self.num_points))
        solar_azimuths = np.zeros((self.num_points))


        for point in range(self.num_points):
            subsolar_longitude = np.deg2rad(-15 * (self.times[point] - (SOL_HRS/2)))

            Sx = np.cos(subsolar_latitude) * np.sin(subsolar_longitude)
            Sy = (np.cos(self.latitude) * np.sin(subsolar_latitude)) - (
                np.sin(self.latitude) * np.cos(subsolar_latitude) * np.cos(subsolar_longitude))
            Sz = (np.sin(np.deg2rad(self.latitude)) * np.sin(subsolar_latitude)) + (
                np.cos(np.deg2rad(self.latitude)) * np.cos(subsolar_latitude) * np.cos(subsolar_longitude))

            solar_altitude = np.rad2deg(np.arcsin(Sz))
            solar_azimuth_south_clockwise = np.rad2deg(np.arctan2(-Sx, Sy))
            solar_azimuth = solar_azimuth_south_clockwise + self.configuration.solar_azimuth


            # Calculate relative air mass for scaling the optical depth of the atmosphere
            relative_air_mass = np.sqrt(np.power(MARS_RADIUS + MARS_ATMOSPHERE_OPTICAL_HEIGHT, 2) - 
                np.power(MARS_RADIUS * np.cos(np.deg2rad(solar_altitude)), 2)) - MARS_RADIUS * np.sin(np.deg2rad(solar_altitude))
            relative_air_mass /= MARS_ATMOSPHERE_OPTICAL_HEIGHT

            optical_depth = np.exp(- self.atmosphere_tau * relative_air_mass)

            solar_intensity = solar_intensity_space * optical_depth

            if Sz <= 0:
                solar_intensity = 0
                solar_altitude = 0

            solar_intensities[point] = solar_intensity
            solar_altitudes[point] = solar_altitude
            solar_azimuths[point] = solar_azimuth
        
        self.solar_intensities = solar_intensities
        self.solar_altitudes = solar_altitudes
        self.solar_azimuths = solar_azimuths

    def generate_temperatures(self, T_peak, T_min, time_peak):
        """
        Generate the temperatures at each point in the day - can be called with either ground or ambient

        Model used is a sinusoidal fit during the day, rising from a minimum at sunset to a maximum
            either at local noon, or time_peak if set
        Between sunset and sunrise, a linear fit is used to bridge the temperatures.
        This approximately matches the data from [4]
        """
        
        axis_obliquity = np.deg2rad(24.936)
        solar_declination = np.arcsin(np.sin(axis_obliquity) * np.sin(np.deg2rad(self.areocentric_longitude)))

        time_sunrise = (SOL_HRS / (2 * np.pi)) * np.arccos(-np.tan(np.deg2rad(self.latitude)) * np.tan(solar_declination))
        time_sunset = SOL_HRS - time_sunrise

        # Setting up amplitude for the cos wave
        T_sine_peak = T_peak 
        T_sine_trough = (T_min - T_sine_peak * np.cos((time_sunrise - time_peak) * 2 * np.pi / SOL_HRS)) / (
            1 - np.cos((time_sunrise - time_peak) * 2 * np.pi / SOL_HRS))
        
        def T_sinusoid(time):
            T = (T_sine_peak - T_sine_trough) * (np.cos((time - time_peak) * 2 * np.pi / SOL_HRS)) + T_sine_trough
            return(T)

        # Set up parameters for the linear fit
        T_sunset = T_sinusoid(time_sunset)
        T_sunrise = T_min
        T_midnight = (T_sunset + T_sunrise) / 2

        def T_linear(time):
            if time > time_sunset:
                T = (T_sunset - T_midnight) * ((SOL_HRS - time) / (SOL_HRS - time_sunset)) + T_midnight
            else:
                T = (T_sunrise - T_midnight) * (time / time_sunrise) + T_midnight
            return(T)
        
        temps = np.zeros((self.num_points))


        for point, time in enumerate(self.times):
            if time > time_sunrise and time < time_sunset:
                # Time is during the day
                T = T_sinusoid(time)
            else:
                T = T_linear(time)
            
            temps[point] = T
        
        return(temps)

        






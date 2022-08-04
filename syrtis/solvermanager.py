"""
Object for sweeping through different parameters in the Solver

"""
from itertools import product
import numpy as np
import copy

from syrtis.configuration import Configuration
from syrtis.solver import Solver

class SolverManager:
    pass

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

    def run_all_configurations(self, verbose=False):
        """
        Run all the configurations
        """

        heat_losses = []
        reports = []
        
        for configuration_index, configuration in enumerate(self.configurations):
            self.habitat.verify_geometry()

            solver = Solver(str(self.each_configuration_inputs_dicts[configuration_index]["name"]),
                self.habitat, configuration)
            
            
            
            if verbose:
                heat_loss, report = solver.solve(verbose=True)

                heat_losses.append(heat_loss)

                reports.append(report)
            
            else:
                heat_loss = solver.solve()
                
                heat_losses.append(heat_loss)
        
        if verbose:
            return(self.each_configuration_inputs_dicts, heat_losses, reports)
        else:
            return(self.each_configuration_inputs_dicts, heat_losses)



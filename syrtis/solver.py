"""
Contains the Object which solves for a given Habitat and Configuration

"""


from habitat import *
from configuration import *

class Solver:
    def __init__(self, name, habitat, configuration):

        assert type(habitat) == Habitat, "'habitat' must be a Habitat object"
        assert type(configuration) == Configuration, "'configuration' must be a Configuration object"

        self.name = name
        self.habitat = habitat
        self.configuration = configuration
    
    def iterate_constant_temperature(self):

        T_internal_start = self.configuration.T_habitat
        Q_initial_start = 1000
        
        Q = Q_initial_start
        shell_temperatures = self.generate_initial_state(T_internal_start)

        current_error = Q

        cutoff_ratio = 1e-5
        target_iterations = 100

        iterations = 0

        while (abs(current_error/Q) > cutoff_ratio) and iterations < target_iterations:

            new_shell_temperatures = self.solve_habitat_conduction(shell_temperatures, Q)

            Q_wall = self.solve_wall_loss(new_shell_temperatures[-1])

            new_error = abs(Q - Q_wall)
            previous_error = current_error
            current_error = new_error
            iterations += 1

            Q = (Q_wall + Q) / 2
            shell_temperatures = new_shell_temperatures
        
        if (abs(current_error/Q) > cutoff_ratio):
            print("Did not converge")
            return(np.nan)
        else:
            return(Q)

    def generate_initial_state(self, T_internal_start):
        # Assume there are (N+1) Shells, set the outermost - not attached to any physical Shell has temperature equal to outside air
        # The rest are linearly spaced up to the internal temperature
        T_external_start = self.configuration.T_air
        N_temperatures = len(self.habitat._shells) + 1

        T_external_wall = ((T_internal_start - T_external_start) / (N_temperatures + 1)) + T_external_start
        
        initial_shell_temperatures = np.linspace(T_internal_start, T_external_wall, N_temperatures)

        return(initial_shell_temperatures)

    def solve_habitat_conduction(self, shell_temperatures, Q):
        
        wall_resistances = self.habitat.build_thermal_resistances(shell_temperatures, self.configuration.GRAVITY)

        updated_shell_temperatures = np.zeros((len(shell_temperatures)))

        # For both calculation types, the interior wall temperature is considered constant in this calculation
        # The temperatures throughout the rest of the Shells are found with the thermal resistances

        internal_temperature = shell_temperatures[0]
        updated_shell_temperatures[0] = internal_temperature

        for shell_count in range(1,len(wall_resistances)+1):

            shell_temperature = updated_shell_temperatures[shell_count-1] - wall_resistances[shell_count-1] * Q

            if shell_temperature < 0:
                shell_temperature = shell_temperatures[shell_count]
                break

            updated_shell_temperatures[shell_count] = shell_temperature
        
        return(updated_shell_temperatures)
    
    def solve_wall_loss(self, wall_temperature):

        Q_wall = self.habitat.placeholder_convective_loss(wall_temperature, self.configuration.T_air)

        return(Q_wall)
    
    def generate_error(self, state1, state2):
        error = np.sum((np.ndarray(state1) - np.ndarray(state2))**2)

        return(error)

            


"""
Contains the Object which solves for a given Habitat and Configuration

"""
from scipy.optimize import minimize
import copy

from syrtis.habitat import *
from syrtis.configuration import *


class Solver:
    def __init__(self, name, habitat, configuration):

        assert type(habitat) == Habitat, "'habitat' must be a Habitat object"
        assert type(configuration) == Configuration, "'configuration' must be a Configuration object"

        assert habitat.verified, "'habitat' must be verified with habitat.verify_geometry()"

        self.name = name
        self.habitat = habitat
        self.configuration = configuration

        self.heat = 0
        self.shell_temperatures = []
        self.report = {
                "Name":       self.name,
                "Total heat flux out": 0,
                "Outer wall temperature": 0,
                "Inner wall temperature": 0,
                "Convective loss from cylinder": 0,
                "Convective loss from endcap": 0,
                "Radiative loss to sky": 0,
                "Radiative loss to ground": 0,
                "Radiative gain from sky": 0,
                "Radiative gain from ground": 0,
                "Direct solar gain": 0,
                "Reflected solar gain": 0,
                "Conduction loss to ground": 0
            }
        self.thermal_energy = 0


    def solve(self):
        if self.configuration.solution_type == "constant temperature":
            self.iterate_constant_temperature()
        elif self.configuration.solution_type == "constant power":
            self.iterate_constant_power()

    def iterate_constant_temperature(self):
        """
        Iterate to find the correct properties for a constant inner wall temperature
        Requires a back-and-forth between conduction power loss and external power loss, which must be balanced
        Typically slower
        """
        T_internal_start = self.configuration.T_habitat
        Q_internal_flux = 1000
        
        Q_internal_flux = Q_internal_flux
        shell_temperatures, R_wall = self.generate_initial_state(T_internal_start)

        current_error = Q_internal_flux

        cutoff_ratio = 1e-5
        target_iterations = 1750

        iterations = 0

        while (abs(current_error) > cutoff_ratio) and iterations < target_iterations:

            new_shell_temperatures, R_wall = self.conduction_temperatures(shell_temperatures, Q_internal_flux)

            Q_external_flux = self.external_losses(new_shell_temperatures[-1], R_wall)

            new_error = abs(Q_internal_flux - Q_external_flux) / Q_internal_flux
            current_error = new_error
            iterations += 1

            Q_internal_flux = (Q_external_flux + (98 * Q_internal_flux)) / 99
            shell_temperatures = new_shell_temperatures
        
        if (abs(current_error) > cutoff_ratio):
            print("Did not converge " + str(abs(current_error)))
            return(np.nan)
        else:
            self.external_losses(shell_temperatures[-1], R_wall, verbose=True)

            self.shell_temperatures = shell_temperatures
            self.report["Inner wall temperature"] = shell_temperatures[0]
            self.heat = Q_external_flux
        
        """inputs = [*[Q_internal_flux], *shell_temperatures]

        minimised_outputs = minimize(self.minimisation, inputs, options={"ftol":cutoff_ratio, "xtol":cutoff_ratio})

        print(minimised_outputs)

        Q_internal_flux = minimised_outputs.x[0]
        shell_temperatures = minimised_outputs.x[1:]

        checked_shell_temperatures = self.conduction_temperatures(shell_temperatures, breakdown["Power loss"])
        print(checked_shell_temperatures)

        return(Q_internal_flux)"""
    
    def iterate_constant_power(self):
        """
        Iterate to find the shell temperatures when power loss from the habitat is constant
        """
        Q_target = self.configuration.Q_habitat
        T_internal_start = 298
        T_internal = T_internal_start
        shell_temperatures, R_wall = self.generate_initial_state(T_internal_start)
        Q_loss = 0
        external_temperature = min(self.configuration.T_ground, self.configuration.T_air)

        current_error = 1

        cutoff_ratio = 1e-10
        target_iterations = 1750

        iterations = 0

        while (abs(current_error) > cutoff_ratio) and iterations < target_iterations:
            Q_loss = self.external_losses(shell_temperatures[-1], R_wall)

            current_error = abs((Q_target - Q_loss) / Q_target)

            if Q_loss < Q_target:
                # Actual and target power have the sign, but power loss is too low
                # External shell temperature needs to be increased
                shell_temperatures[-1] *= (1 + 0.02*current_error)
            elif Q_loss > Q_target:
                # Actual and target power have the sign, but power loss is too high
                # External shell temperature needs to be decreased
                shell_temperatures[-1] *= (1 - 0.02*current_error)
            
            new_shell_temperatures, R_wall = self.conduction_temperatures(shell_temperatures, Q_target)

            shell_temperatures = new_shell_temperatures

            #if iterations%10==0:print(shell_temperatures[0], shell_temperatures[-1], Q_loss, Q_target, current_error, R_wall)

            iterations += 1

        if (abs(current_error) > cutoff_ratio):
            print("Did not converge " + str(abs(current_error)))
            self.heat = np.nan
            return(np.nan)
        else:
            self.external_losses(shell_temperatures[-1], R_wall, verbose=True)
            
            self.report["Inner wall temperature"] = shell_temperatures[0]

            self.shell_temperatures = shell_temperatures
            self.heat = Q_target

    def verify_temperatures(self, shell_temperatures):
        """
        Verify that a given set of temperatures are either monotonically increasing or decreasing 
        This is a reasonable metric for a good solution, as it implies that flux is consistent through the Habitat
        May be violated in future if heating/cooling jackets are added as a Shell subclass

        Also verify that the external wall temperature is between the internal and environment temperatures
        """
        non_increasing = all(x>=y for x, y in zip(shell_temperatures, shell_temperatures[1:]))
        non_decreasing = all(x<=y for x, y in zip(shell_temperatures, shell_temperatures[1:]))

        temperature_correct = True
        
        if self.configuration.solution_type == "constant temperature":
            if self.configuration.T_habitat > self.configuration.T_air:
                # Habitat should be losing heat to the environment
                # External shell temperature should not be substantially hotter than the habitat
                if shell_temperatures[-1] > (self.configuration.T_habitat + 50):
                    temperature_correct = False
            
            elif self.configuration.T_habitat < self.configuration.T_air:
                # Habitat should be strictly gaining heat from the environment
                # External shell temperature should be greater than the habitat
                if shell_temperatures[-1] < self.configuration.T_habitat:
                    temperature_correct = False

        return((non_increasing or non_decreasing) and temperature_correct)

    def minimisation(self, inputs):
        """
        Minimise-able function, called into scipy.optimize solvers

        Args:
            inputs (list of floats):    list of all inputs to allow solving
                inputs[0]:              power loss from habitat as calculated by solve_wall_loss
                inputs[1:n]             shell temperatures for shells [0:n-1]
        """
        Q_internal_flux = inputs[0]
        shell_temperatures = inputs[1:]

        temperature_verification = self.verify_temperatures(shell_temperatures)
        if temperature_verification == False:
            return(1e99)

        else:
            new_shell_temperatures = self.conduction_temperatures(shell_temperatures, Q_internal_flux)

            Q_external_flux = self.external_losses(new_shell_temperatures[-1])

            error = abs(Q_internal_flux - Q_external_flux)**2

            return(error)

    def generate_initial_state(self, T_internal_start, T_external_start=0):
        # Assume there are (N+1) Shells, set the outermost - not attached to any physical Shell has temperature equal to outside air
        # The rest are linearly spaced up to the internal temperature
        if T_external_start == 0:
            T_external_start = self.configuration.T_air

        N_temperatures = len(self.habitat._shells) + 1

        T_external_wall = ((T_internal_start - T_external_start) / (N_temperatures + 1)) + T_external_start
        
        initial_shell_temperatures = np.linspace(T_internal_start, T_external_wall, N_temperatures)

        return(initial_shell_temperatures, 1)

    def conduction_temperatures(self, shell_temperatures, Q):
        """
        Find the updated temperatures across the habitat wall, based on a set of starting temperatures and the heat flux
        """

        wall_resistances = self.habitat.build_thermal_resistances(shell_temperatures, self.configuration.GRAVITY)

        total_resistance = sum(wall_resistances)

        updated_shell_temperatures = shell_temperatures

        # For both calculation types, the interior wall temperature is considered constant in this calculation
        # The temperatures throughout the rest of the Shells are found with the thermal resistances

        if self.configuration.solution_type == "constant temperature":
            # Iterate outwards from the inside to the outside - appropriate for a constant-temperature type solution
            updated_shell_temperatures[0] = self.configuration.T_habitat

            for shell_count in range(1,len(wall_resistances)+1):

                shell_temperature = updated_shell_temperatures[shell_count-1] - wall_resistances[shell_count-1] * Q
    
                if shell_temperature <= 0:
                    shell_temperature = shell_temperatures[shell_count]

                updated_shell_temperatures[shell_count] = shell_temperature
    
        elif self.configuration.solution_type == "constant power":
            # Assume the outermost shell temperature is correct and iterate inwards
            updated_shell_temperatures[-1] = shell_temperatures[-1]

            for shell_count in reversed(range(0, len(wall_resistances))):
                shell_temperature = updated_shell_temperatures[shell_count+1] + wall_resistances[shell_count] * Q

                if shell_temperature <= 0:
                    # This means power loss is too high for the wall geometry used
                    updated_shell_temperatures = [shell_temperatures[0] for i in shell_temperatures]
                    break
                
                updated_shell_temperatures[shell_count] = shell_temperature

        return(updated_shell_temperatures, total_resistance)
    
    def external_losses(self, wall_temperature, thermal_resistance_wall, verbose=False):

        #Q_wall = self.habitat.placeholder_convective_loss(wall_temperature, self.configuration.T_air)
        Q_wall = 0
        Q_endcap = 0
        Q_rad_sky_out = 0
        Q_rad_sky_in = 0
        Q_rad_ground_out = 0
        Q_rad_ground_in = 0
        Q_solar_direct = 0
        Q_solar_indirect = 0
        Q_conduction = 0

        if (self.configuration.air_direction == "cross" and self.habitat.orientation == "horizontal") or (
            self.habitat.orientation == "vertical"):
            Q_wall = self.habitat.convective_loss_cylinder_cross(self.configuration.air, 
                                                            self.configuration.v_air,
                                                            self.configuration.T_air,
                                                            wall_temperature)
            
            Q_endcap = self.habitat.convective_loss_endcap_cross(self.configuration.air, 
                                                            self.configuration.v_air,
                                                            self.configuration.T_air,
                                                            wall_temperature)
        
        elif (self.configuration.air_direction == "axial" and self.habitat.orientation == "horizontal"):
            Q_wall = self.habitat.convective_loss_cylinder_axial(self.configuration.air, 
                                                            self.configuration.v_air,
                                                            self.configuration.T_air,
                                                            wall_temperature)
            
            Q_endcap = self.habitat.convective_loss_endcap_axial(self.configuration.air, 
                                                            self.configuration.v_air,
                                                            self.configuration.T_air,
                                                            wall_temperature)
        
        
        # Sign convention is preserved by functions in Habitat
        Q_rad_sky_out = self.habitat.radiative_loss_sky(wall_temperature)
        Q_rad_sky_in = self.habitat.radiative_gain_sky(self.configuration.T_air)

        Q_rad_ground_out = self.habitat.radiative_loss_ground(wall_temperature)
        Q_rad_ground_in = self.habitat.radiative_gain_ground(self.configuration.T_ground)
        
        if self.configuration.solar_altitude > 0:
            Q_solar_direct = self.habitat.solar_gain_direct(self.configuration.solar_altitude,
                                                            self.configuration.solar_azimuth,
                                                            self.configuration.solar_intensity)
            
            Q_solar_indirect = self.habitat.solar_gain_indirect(self.configuration.solar_altitude,
                                                            self.configuration.solar_azimuth,
                                                            self.configuration.solar_intensity,
                                                            self.configuration.albedo_ground)

        Q_convective = Q_wall + Q_endcap
        Q_rad_environment_out = Q_rad_sky_out + Q_rad_ground_out
        Q_rad_environment_in = Q_rad_sky_in + Q_rad_ground_in

        if self.habitat.groundlevel != None:
            if self.habitat.groundlevel.thermal_resistance != 0:
                Q_conduction = self.habitat.conductive_loss_fixed_resistance(wall_temperature,
                                                                            self.configuration.T_ground,
                                                                            self.habitat.groundlevel.thermal_resistance)
            else:
                if self.habitat.orientation == "horizontal":
                    Q_conduction += self.habitat.conductive_loss_horizontal_cylinder_steady(wall_temperature,
                                                                                            self.configuration.T_ground,
                                                                                            self.configuration.k_ground,
                                                                                            thermal_resistance_wall)
                elif self.habitat.orientation == "vertical":
                    Q_conduction += self.habitat.conductive_loss_vertical_cylinder_steady(wall_temperature,
                                                                                        self.configuration.T_ground,
                                                                                        self.configuration.k_ground)
                
                if self.habitat.endcap_type == "hemisphere":
                    Q_conduction += self.habitat.conductive_loss_hemisphere_steady(wall_temperature,
                                                                                self.configuration.T_ground,
                                                                                self.configuration.k_ground)
                elif self.habitat.endcap_type == "flat":
                    Q_conduction += self.habitat.conductive_loss_disc_steady(wall_temperature,
                                                                            self.configuration.T_ground,
                                                                            self.configuration.k_ground)


        Q_total = (Q_convective + 
        Q_rad_environment_out + Q_rad_environment_in + 
        Q_solar_direct + Q_solar_indirect + 
        Q_conduction)
        
        if verbose:
            self.report["Total heat flux out"] = Q_total
            self.report["Outer wall temperature"] = wall_temperature
            self.report["Convective loss from cylinder"] = Q_wall
            self.report["Convective loss from endcap"] = Q_endcap
            self.report["Radiative loss to sky"] = Q_rad_sky_out
            self.report["Radiative loss to ground"] = Q_rad_ground_out
            self.report["Radiative gain from sky"] = Q_rad_sky_in
            self.report["Radiative gain from ground"] = Q_rad_ground_in
            self.report["Direct solar gain"] = Q_solar_direct
            self.report["Reflected solar gain"] = Q_solar_indirect
            self.report["Conduction loss to ground"] = Q_conduction
        else:
            return(Q_total)
    
    def calculate_thermal_energy(self, T_ref=0):
        """
        Calculates the total thermal energy stored in the habitat, compared to a reference temperature

        Args:
            T_ref (float):          reference temperature for the energy state (K)
        """
        assert len(self.shell_temperatures) == len(self.habitat._shells) + 1, "Run solver.solve before solver.thermal_energy"

        self.thermal_energy = 0

        for shell_count, shell in enumerate(self.habitat._shells):
            T_avg = sum(self.shell_temperatures[shell_count:shell_count+1])
            self.thermal_energy += shell.thermal_energy(self.habitat.endcap_type, T_avg, T_ref)

    def generate_error(self, state1, state2):
        error = np.sum((np.ndarray(state1) - np.ndarray(state2))**2)

        return(error)

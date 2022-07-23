"""
Contains the Object which solves for a given Habitat and Configuration

"""
from scipy.optimize import minimize



class Solver:
    def __init__(self, name, habitat, configuration):

        assert type(habitat) == Habitat, "'habitat' must be a Habitat object"
        assert type(configuration) == Configuration, "'configuration' must be a Configuration object"

        assert habitat.verified, "'habitat' must be verified with habitat.verify_geometry()"

        self.name = name
        self.habitat = habitat
        self.configuration = configuration
    
    def iterate_constant_temperature(self):

        T_internal_start = self.configuration.T_habitat
        Q_initial_start = 1000
        
        Q_loss = Q_initial_start
        shell_temperatures = self.generate_initial_state(T_internal_start)

        current_error = Q_loss

        cutoff_ratio = 1e-4
        target_iterations = 1000

        iterations = 0

        while (abs(current_error/Q_loss) > cutoff_ratio) and iterations < target_iterations:

            new_shell_temperatures = self.solve_habitat_conduction(shell_temperatures, Q_loss)

            Q_wall = self.solve_wall_loss(new_shell_temperatures[-1])

            #print(current_error, Q_loss, Q_wall, new_shell_temperatures[-1])

            new_error = abs(Q_loss - Q_wall)
            current_error = new_error
            iterations += 1

            Q_loss = (Q_wall + (99 * Q_loss)) / 100
            shell_temperatures = new_shell_temperatures
        
        print("a")
        
        if (abs(current_error/Q_loss) > cutoff_ratio):
            print("Did not converge " + str(abs(current_error/Q_loss)))
            return(np.nan)
        else:
            return(Q_loss)
        
        """inputs = [*[Q_loss], *shell_temperatures]

        minimised_outputs = minimize(self.minimisation, inputs, method="Powell", options={"ftol":cutoff_ratio, "xtol":cutoff_ratio})

        print(minimised_outputs)

        Q_loss = minimised_outputs.x[0]

        return(Q_loss)"""

    def minimisation(self, inputs):
        """
        Minimise-able function, called into scipy.optimize solvers

        Args:
            inputs (list of floats):    list of all inputs to allow solving
                inputs[0]:              power loss from habitat as calculated by solve_wall_loss
                inputs[1:n]             shell temperatures for shells [0:n-1]
        """
        Q_loss = inputs[0]
        shell_temperatures = inputs[1:]

        new_shell_temperatures = self.solve_habitat_conduction(shell_temperatures, Q_loss)

        Q_wall = self.solve_wall_loss(new_shell_temperatures[-1])

        error = abs(Q_loss - Q_wall)

        return(error)

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

        updated_shell_temperatures = shell_temperatures

        # For both calculation types, the interior wall temperature is considered constant in this calculation
        # The temperatures throughout the rest of the Shells are found with the thermal resistances

        internal_temperature = shell_temperatures[0]
        updated_shell_temperatures[0] = internal_temperature

        for shell_count in range(1,len(wall_resistances)+1):

            shell_temperature = updated_shell_temperatures[shell_count-1] - wall_resistances[shell_count-1] * Q

            if shell_temperature < 0:
                #print("Shell temperature error ", str(shell_temperature))
                shell_temperature = shell_temperatures[shell_count] / 2
                
                #break

            updated_shell_temperatures[shell_count] = shell_temperature
        
        return(updated_shell_temperatures)
    
    def solve_wall_loss(self, wall_temperature, report_full=False):

        #Q_wall = self.habitat.placeholder_convective_loss(wall_temperature, self.configuration.T_air)
        Q_wall = 0
        Q_endcap = 0
        Q_rad_sky_out = 0
        Q_rad_sky_in = 0
        Q_rad_ground_out = 0
        Q_rad_ground_in = 0
        Q_solar_direct = 0
        Q_solar_indirect = 0

        
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
        Q_rad_ground_in = self.habitat.radiative_loss_ground(self.configuration.T_ground)
        
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

        return(Q_convective + Q_rad_environment_out + Q_rad_environment_in + Q_solar_direct + Q_solar_indirect)
    
    def generate_error(self, state1, state2):
        error = np.sum((np.ndarray(state1) - np.ndarray(state2))**2)

        return(error)

            
"""
if __name__ == "__main__":
    steel = Solid("Steel", 150, 8700, 500, 0.55)
    co2 = ConstrainedIdealGas("STP CO2", 101325, 44, 0.71, 10.9e-6, 749, 0.0153)

    bocachica = Configuration("bocachica", "constant temperature",
    300, 1, 0.29, 300, 101325, 1, "cross", 90, 90, 1000, T_habitat=80)
    
    heat_gain = []
    thicknesses = np.logspace(-3, 0, 10)
    #thicknesses = np.linspace(0.001, 1.001, 500)

    for thickness in thicknesses:
        tankfarm = Habitat("vertical", 50, "flat")
        tankfarm.create_static_shell(co2, 4.5)
        tankfarm.create_static_shell(steel, 0.001)
        tankfarm.create_static_shell(co2, thickness)
        tankfarm.create_static_shell(steel, 0.001)

        tankfarm.verify_geometry()

        s = Solver("test", tankfarm, bocachica)
        q = s.iterate_constant_temperature()

        heat_gain.append(-q)
    
    import matplotlib.pyplot as plt

    plt.scatter(thicknesses, heat_gain)
    plt.xscale("log")
    plt.xlim(9e-4, 1.5)
    plt.xlabel("Gap between inner and outer wall (m)")
    plt.ylabel("Heat gain into tank (W)")
    plt.title("Syrtis evaluation case \n Heat gain into GSE tanks at Boca Chica tank farm")
    plt.show()"""


"""if __name__ == "__main__":
    steel = Solid("Steel", 150, 8700, 500, 0.55)
    painted_steel = Solid("Painted Steel", 150, 8700, 500, 0.1)
    internal_air = ConstrainedIdealGas("STP CO2", 101325, 29, 0.71, 10.9e-6, 749, 0.0153)
    co2_ambient = ConstrainedIdealGas("STP CO2", 580, 44, 0.71, 10.9e-6, 749, 0.0153)

    columbus = Habitat("horizontal", 7, "flat")

    columbus.create_static_shell(internal_air, 2.2)
    columbus.create_static_shell(steel, 4e-3)

    columbus.verify_geometry()

    columbus_p = Habitat("horizontal", 7, "flat")

    columbus_p.create_static_shell(internal_air, 2.2)
    columbus_p.create_static_shell(painted_steel, 4e-3)

    columbus_p.verify_geometry()

    temps = np.linspace(273, 313, 40)
    temps_c = np.linspace(0, 40, 40)
    qs = []
    qs_painted = []

    for temp in temps:
        equator = Configuration("equator", "constant temperature",
            210, 0.1, 0.29, 210, 580, 1, "cross", 90, 90, 605, T_habitat=temp)
    
        s = Solver("columbus equator Mars", columbus, equator)

        q = s.iterate_constant_temperature()
    
        qs.append(q)

        s_p = Solver("columbus equator Mars", columbus_p, equator)

        q_p = s_p.iterate_constant_temperature()
    
        qs_painted.append(q_p)
    
    print(qs[0])

    import matplotlib.pyplot as plt

    plt.scatter(qs, temps_c, label="Unpainted steel")
    plt.scatter(qs_painted, temps_c, label="Painted steel")
    plt.ylabel("Internal temperature (C)")
    plt.xlabel("Heat loss")
    plt.title("Syrtis evaluation case \n Heat loss from ISS Columbus on Martian surface")
    plt.legend()
    plt.show()
"""
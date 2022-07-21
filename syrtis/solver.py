"""
Contains the Object which solves for a given Habitat and Configuration

"""


from habitat import *
from configuration import *

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
        
        Q = Q_initial_start
        shell_temperatures = self.generate_initial_state(T_internal_start)

        current_error = Q

        cutoff_ratio = 1e-4
        target_iterations = 200

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

        updated_shell_temperatures = shell_temperatures

        # For both calculation types, the interior wall temperature is considered constant in this calculation
        # The temperatures throughout the rest of the Shells are found with the thermal resistances

        internal_temperature = shell_temperatures[0]
        updated_shell_temperatures[0] = internal_temperature

        for shell_count in range(1,len(wall_resistances)+1):

            shell_temperature = updated_shell_temperatures[shell_count-1] - wall_resistances[shell_count-1] * Q

            if shell_temperature < 0:
                shell_temperature = shell_temperatures[shell_count]
                print("Shell temperature error")
                break

            updated_shell_temperatures[shell_count] = shell_temperature
        
        return(updated_shell_temperatures)
    
    def solve_wall_loss(self, wall_temperature):

        #Q_wall = self.habitat.placeholder_convective_loss(wall_temperature, self.configuration.T_air)
        Q_wall = 0
        Q_endcap = 0
        Q_rad_sky = 0
        Q_rad_ground = 0
        
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
        
        Q_rad_sky = self.habitat.radiative_loss_sky(self.configuration.T_air, wall_temperature)
        Q_rad_ground = self.habitat.radiative_loss_ground(self.configuration.T_ground, wall_temperature)

        return(Q_wall + Q_endcap + Q_rad_sky + Q_rad_ground)
    
    def generate_error(self, state1, state2):
        error = np.sum((np.ndarray(state1) - np.ndarray(state2))**2)

        return(error)

            

"""if __name__ == "__main__":
    steel = Solid("Steel", 150, 8700, 500, 0.55)
    co2 = ConstrainedIdealGas("STP CO2", 101325, 44, 0.71, 10.9e-6, 749, 0.0153)

    equator = Configuration("equator", "constant temperature",
        210, 0.1, 210, 580, 1, "cross", 90, 90, 580, T_habitat=190)

    bocachica = Configuration("bocachica", "constant temperature",
    300, 1, 300, 101325, 1, "cross", 90, 90, 604, T_habitat=80)
    
    heat_gain = []
    thicknesses = np.logspace(-3, 0, 100)
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


if __name__ == "__main__":
    steel = Solid("Steel", 150, 8700, 500, 0.55)
    internal_air = ConstrainedIdealGas("STP CO2", 101325, 29, 0.71, 10.9e-6, 749, 0.0153)
    co2_ambient = ConstrainedIdealGas("STP CO2", 580, 44, 0.71, 10.9e-6, 749, 0.0153)

    columbus = Habitat("horizontal", 7, "flat")

    columbus.create_static_shell(internal_air, 2.2)
    columbus.create_static_shell(steel, 4e-3)
    columbus.create_static_shell(co2_ambient, 0.04)
    columbus.create_static_shell(steel, 4e-3)

    columbus.verify_geometry()


    temps = np.linspace(273, 313, 80)
    temps_c = np.linspace(0, 40, 80)
    qs = []

    for temp in temps:
        equator = Configuration("equator", "constant temperature",
            210, 0.1, 210, 580, 1, "cross", 90, 90, 580, T_habitat=temp)
    
        s = Solver("columbus equator Mars", columbus, equator)

        q = s.iterate_constant_temperature()
    
        qs.append(q)

    import matplotlib.pyplot as plt

    plt.scatter(qs, temps_c)
    plt.ylabel("Internal temperature (C)")
    plt.xlabel("Heat loss")
    plt.title("Syrtis evaluation case \n Heat loss from ISS Columbus on Martian surface")
    plt.show()

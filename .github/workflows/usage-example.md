```python
from syrtis import *

aluminium = Solid("Aluminium", 247, 2700, 900, absorb=0.15, emit=0.04)
plastic = Solid("Generic plastic", 10, 1300, 1420, absorb=0.89, emit=0.84)

internal_air = ConstrainedIdealGas("STP Air", 101325, 29, 0.71, 17.9e-6, 1010, 0.0252)
martian_air = ConstrainedIdealGas("Martian ambient pressure CO2", 580, 44, 0.71, 10.9e-6, 749, 0.0153)

equator = Configuration("Martian equator at noon",
 "constant temperature", 210, 0.2, 0.29, 210, 580, 5, "cross", 90, 90, 605, T_habitat=290)
```

```python
HAB_vertical = Habitat(orientation="vertical", length=8, endcap_type="flat")

HAB_vertical.create_static_shell(internal_air, 4.400)
HAB_vertical.create_static_shell(aluminium, 4.8e-3)

HAB_vertical.create_static_shell(martian_air, 50e-3, parallel_thermal_resistance=8.5e-5)
# Parallel resistance corresponds to 1% of the cross-section being aluminium

HAB_vertical.create_static_shell(plastic, 12e-3)
HAB_vertical.create_static_shell(martian_air, 50e-3, parallel_thermal_resistance=8.4e-5)

HAB_vertical.create_static_shell(aluminium, 2e-3)

HAB_vertical.create_ground_level(thermal_resistance=1)
# Thermal resistance roughly corresponds to six aluminium landing legs, each 5m long and with 50cm2 area
```

```python
cm_vertical = ConfigurationManager(HAB_vertical, equator, {"T_habitat":list(range(273, 313, 1))})
configs_vertical, heats_vertical, reports_vertical = cm_vertical.run_all_configurations(verbose=True)
```
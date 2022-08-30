---
layout: page
title: "usage-example"
permalink: /usage/
---

Import `syrtis` and other relevant modules

```python
from syrtis import *
import numpy as np
import matplotlib.pyplot as plt
```

Create the `Material` objects needed for an example habitat

```python
aluminium = Solid("Aluminium", k=247, rho=2700, cp=900, absorb=0.15, emit=0.04)
plastic = Solid("Generic plastic", k=10, rho=1300, cp=1420, absorb=0.89, emit=0.84)

internal_air = ConstrainedIdealGas("STP Air", 
    p=101325, M=29, Pr=0.71, mu=17.9e-6, cp_ref=1010, k=0.0252)
martian_air = ConstrainedIdealGas("Martian ambient pressure CO2", 
    p=580, M=44, Pr=0.71, mu=10.9e-6, cp_ref=749, k=0.0153)
```

Create the `Configuration` object to store the boundary conditions of the simulation

```python
equator = Configuration("Martian equator at noon", solution_type="constant temperature", 
    T_ground=210, k_ground=0.2, albedo_ground=0.29, 
    T_air=210, p_air=580, v_air=5, air_direction="cross", 
    solar_altitude=90, solar_azimuth=90, solar_intensity=605, 
    T_habitat=290)
```

Create the `Habitat` object and add all the geometry features we want to analyse

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

Create a `ConfigurationManager` to sweep through a range of habitat internal temperatures

```python
cm_vertical = ConfigurationManager(HAB_vertical, equator, {"T_habitat":list(range(273, 313, 1))})
configs_vertical, heats_vertical, reports_vertical = cm_vertical.run_all_configurations(verbose=True)
```

Plot the results of the temperature sweep

```python
plt.scatter(heats_vertical_coated, habitat_temps, label="Vertical configuration")

plt.xlabel("Heat flux (W)")
plt.ylabel("Habitat internal temperature (C)")
plt.title("Internal temperature vs heat loss")

plt.legend()
plt.show()
```

Plot a breakdown of heat fluxes where the internal temperature is 20 degrees

```python

plt.figure(figsize=(8,6))
plt.title("Heat loss breakdown in Mars Direct-style Habitats")
tools.plot_power_balance([
    reports_vertical_coated[i] for i in range(len(reports_vertical_coated)) 
    if (configs_vertical_coated[i]["T_habitat"] == 293)], 
labels=["Vertical"])

```
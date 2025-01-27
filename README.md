# NM-project
_Repository for the project of Mathematical Modeling and Simulation, held by prof. Marchetti (M.Sc. in Quantitative and Computational Biology, 2024-2025)._

[MOESM4_RSSA_2.m](https://github.com/Sara-Baldinelli/NM-project/blob/main/code/MOESM4_RSSA_2.m) implements the RSSA simulation and [histograms.m](https://github.com/Sara-Baldinelli/NM-project/blob/main/code/histograms.m) generates histograms for RSSA and RTC comparison.

### MOESM4_RSSA_2.m example usage:
```
[T, Dynamics] = MOESM4_RSSA_2(4000, 40, 40, -50, 0.1, 0.0001); 
```
_(Maximum Simulation Time, Calcium channels, Potassium channels, Voltage, delta, step)_

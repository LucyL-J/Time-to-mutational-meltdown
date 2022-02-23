# Time-to-mutational-meltdown
Documentation of the simulation and the computational analysis of the manuscript "The time to mutational meltdown".

The simulation of the population dynamics and mutation accumulation is done in the script `simulate_mutational_meltdown.jl`. First, the parameters are initialised and a founder population created, followed by repeated mutation, reproduction and population-size control steps. Possible outputs are the extinction time, the number or the distribution of mutations every generation, the beginning and ending of pre-ratchet, ratchet and meltdown phase and the times until the fittest class is lost.

We run the simulations for the parameter regime described in the manuscript using the script `launcher_script.jl`. Here, we are interested in the times of pre-ratchet, ratchet and meltdown phase and the times until loss of the fittest class. The outputs are a data frame with the simulation output and another data frame with the parameters used for this specific simulation. They can be identified via the seed from the random number generator.

The analysis is done in the jupyter notebook `Computational_analysis.ipynb`. This notebook is self-explanatory.

![Flow chart](flowchart.png)


The simulations are done in the programming language [julia](https://julialang.org).
To replicate the data in the manuscript, install julia, start it and run
```
include("launcher_script.jl")
ratchet_speeds_report()
times_phases_report()
```
To start the notebook, run
```
using IJulia
notebook()
```

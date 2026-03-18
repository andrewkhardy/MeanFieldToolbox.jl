# MeanFieldToolkit.jl

[![Build Status](https://https://github.com/Toronto-Condensed-Matter-Theory/MeanFieldToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://https://github.com/Toronto-Condensed-Matter-Theory/MeanFieldToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)

MeanFieldToolkit.jl is a Julia package meant for solving generalized self-consistent mean-field equations on a lattice. It is an updated version of [MeanFieldToolkit.jl](https://github.com/Anjishnubose/MeanFieldToolkit.jl)

Currently supported :
* Lattice implementation is done through [TightBindingToolkit.jl](https://https://github.com/Toronto-Condensed-Matter-TheoryTightBindingToolkit.jl). Any custom lattice in $d=1,2,3$ is supported.
* User can input any two-site interaction in the form of arrays, and their corresponding mean-field equations. Simple four-fermion interactions are already built in (such as Hubbard, Spin-Spin interactions etc.).
* Can track any hopping and pairing order parameters.
* Self-consistentcy solver is implemented using [FixedPointToolkit.jl](https://https://github.com/Toronto-Condensed-Matter-TheoryFixedPointToolkit.jl). One can customize the solver, the tolerance of convergence, the maximum number of iterations and so on.
* Can checkpoint and save results into JLD2 files, and resume iterations from reading such files.
* Can plot results of order parameters, and the mean-field ground state energy as a function of iterations.

# Documentation

For further details, please refer to the [Documentation](https://Toronto-Condensed-Matter-Theory.github.io/MeanFieldToolkit.jl/dev/).

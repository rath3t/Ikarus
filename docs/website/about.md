# About

## What is Ikarus?

Ikarus is a project developed by the Institute for Structural Mechanics, University of Stuttgart.
It strives to develop an easy-to-read and easy-to-use framework to perform finite element analysis.
It acts as a front-end code that utilizes several modules from [DUNE](https://dune-project.org/)
and has some other dependencies.
<figure class="inline end" markdown>
![Ikarus Logo](auxiliaryImages/BigLogo_transparent.png)
</figure>
Ikarus endeavors to not *reinvent the wheel*
but to just provide a generic interface that makes use of the existing modules, resulting in a more expressive,
modular library to solve partial differential equations.
The main idea is to provide users with a platform to rapidly prototype their ideas,
thereby enabling researchers across the world to quickly test their various theories.
One of the major principles that we follow is *"take what you need,"*
where the users can pick different functionalities as per their requirements to try out their own examples.
We also provide certain [examples](02_examples/index.md#examples) that can help a user get familiar with Ikarus.
Ikarus itself is written in C++, but it also comes with [Python bindings](https://pypi.org/project/pyikarus/),
which can then aid in combining the functionalities of Ikarus with different Python packages.
It uses C++ template (meta) programming to achieve a general and easy-to-extend library.
The idea of interface segregation between a grid space and an ansatz space from DUNE is heavily exploited here.
Some of the features that we provide are

- Element library (linear elastic, nonlinear elastic, Kirchhoff-Love Shell, truss)
- Material library
    - St. Venant-Kirchhoff, Neo-Hookean
    - Hyperelastic materials constructed using various deviatoric and volumetric functions.
        - Included deviatoric functions are Arruda-Boyce, Blatz-Ko, Gent, Invariant-based, and Ogden models.
        - The invariant-based model is a general model and can be used to create other material models,
        for example, Mooney-Rivlin and Yeoh models.
- Nonlinear Solvers (Newton-Raphson Method, Trust region Method)
- Assemblers for sparse and dense matrices
- Control routines (load control, displacement control, arc-length control)
- Handle Dirichlet boundary conditions
- Observer patterns to log user-desired information as an output

## Dependencies

!!! bug
    This section will be updated.

## Ikarus Core Developers

The Ikarus Core Developers develop, maintain and look after Ikarus.

- [Alexander Müller](https://www.ibb.uni-stuttgart.de/en/institute/team/Mueller-00006/)
- [Tarun Kumar Mitruka Vinod Kumar Mitruka](https://www.ibb.uni-stuttgart.de/en/institute/team/Vinod-Kumar-Mitruka/)
- [Henrik Jakob](https://www.ibb.uni-stuttgart.de/en/institute/team/Jakob-00003/)

## Talks

- **Jakob, H.**, Vinod Kumar Mitruka, T. K. M., Müller, A.
  *Ikarus and dune-iga: Easy-To-Use C++ Libraries With Python Bindings for Structural Analysis Within DUNE*.
  9th European Congress on Computational Methods in Applied Sciences and Engineering (ECCOMAS),
  Lisbon, Portugal, June 03–07, 2024.
- **Jakob, H.**, Vinod Kumar Mitruka, T. K. M., Müller, A.
  *Ikarus and dune-iga: Easy-To-Use C++ Libraries With Python Bindings.*
  FE im Schnee 2024,
  Hirschegg, Austria, March 17–20, 2024.
- **Müller, A.**, Jakob, H., Vinod Kumar Mitruka, T. K. M., Sander, O., Bischoff, M.
  *Ikarus with dune-iga: An open-source C++ library for Isogeometric Analysis within the Dune framework.*
  11th International Conference on Isogeometric Analysis 2023 (IGA 2023),
  Lyon, France, June 18–21, 2023.

## Datasets

Different versions of Ikarus are released via the Data Repository of the University of Stuttgart ([DaRUS](https://darus.uni-stuttgart.de/)).

- [Ikarus v0.4](https://doi.org/10.18419/darus-3889)
- [Ikarus v0.3](https://doi.org/10.18419/darus-3303)

## How to cite us?

Here is the BibTex entry for our latest dataset:

```text
@data{Ikarusv04_2024,
author    = {Müller, Alexander and Vinod Kumar Mitruka, Tarun Kumar Mitruka and Jakob, Henrik},
publisher = {DaRUS},
title     = {Ikarus v0.4},
year      = {2024},
version   = {V1},
doi       = {10.18419/darus-3889}
}
```

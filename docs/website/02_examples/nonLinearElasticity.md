---
status: new
---

# Non-linear Elasticity for 2D solids

## Description

In `iks006_nonlinear2DSolid.cpp`, an automatic differentiation-based implementation is used to perform a non-linear analysis on a 2D block.
Various methods to obtain a 2D grid via Dune are also shown in the commented section at
the beginning. Python is used to provide a Neumann boundary condition, providing a demonstration for the usage of
Python-based code within the Ikarus framework. The load control method is chosen as the desired control routine, and
Newton-Raphson (or trust region methods) are used to solve the non-linear problem itself.

## Code highlights

This example uses two macros, `gridType` and `solverType`, that are to be set as desired before executing the example.
The `gridType` can be set to 0, 1, or 2, denoting an `ALUGrid`, a `YaspGrid`, and a `NURBSGrid` respectively.
The `solverType` can be set to either 0 or 1 for `Newton-Raphson` and `Trust region` methods respectively.

When `ALUGrid` is chosen, the mesh file `auxiliaryFiles/unstructuredTrianglesfine.msh` is read using `Dune::GmshReader`
and is then globally refined once. The `YaspGrid` created a square block of length 1 with 10 elements in either direction.
The `NURBSGrid` also creates a square block of length 1 with a polynomial degree of 2 in both directions.
The block is subjected to a Neumann load on the right edge ($x=1$) that is incorporated using the dune-python interface as shown below:

```cpp
Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
std::string lambdaNeumannVertices = std::string("lambda x: ( x[0]>0.999 )");

Python::start();
Python::Reference main = Python::import("__main__");
Python::run("import math");
Python::runStream() << std::endl << "import sys" << std::endl << "import os" << std::endl;
auto pythonNeumannVertices = Python::make_function<bool>(Python::evaluate(lambdaNeumannVertices));

for (auto &&vertex : vertices(gridView)) {
  bool isNeumann                          = pythonNeumannVertices(vertex.geometry().corner(0));
  neumannVertices[indexSet.index(vertex)] = isNeumann;
}
BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);
```

After the basis is defined, the non-linear elastic finite element is created as shown below:

```cpp
auto volumeLoad = [](auto &globalCoord, auto &lamb) {
  Eigen::Vector2d fext;
  fext.setZero();
  return fext;
};

auto neumannBoundaryLoad = [](auto &globalCoord, auto &lamb) {
  Eigen::Vector2d fext;
  fext.setZero();
  fext[1] = lamb / 40;
  return fext;
};

auto matParameter = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});

Ikarus::StVenantKirchhoff matSVK(matParameter);
auto reducedMat = planeStress(matSVK);

std::vector<Ikarus::NonLinearElastic<decltype(basis), decltype(reducedMat)>>> fes;
for (auto &element : elements(gridView))
  fes.emplace_back(*basis, element, reducedMat, &neumannBoundary, neumannBoundaryLoad, volumeLoad);
```

The functors `volumeLoad` and `neumannBoundaryLoad` are used to obtain the external volume and surface loads acting on a particular
position.
We use a Saint Venant–Kirchhoff material model, which we transform to a plane stress material law for our two-dimensional simulation.
The line $y=0$ is clamped by applying the Dirichlet boundary condition expressed below:

```cpp
auto basisP = std::make_shared<const decltype(basis)>(basis);
Ikarus::DirichletValues dirichletValues(basisP->flat());

dirichletValues.fixBoundaryDOFs([&](auto &dirichletFlags, auto &&localIndex, auto &&localView, auto &&intersection) {
  if (std::abs(intersection.geometry().center()[1]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
});
```

The finite element requirements are defined by using the affordance `#!cpp Ikarus::AffordanceCollections::elastoStatics`.
This is then used to create functors to get the stiffness matrix, residual vector, and energy value using a sparse assembler.
A non-linear operator and the linear solver used by the `solverType` are defined as:

```cpp
auto nonLinOp = Ikarus::makeNonLinearOperator(functions(energyFunction, residualFunction, KFunction), d);
auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::sd_UmfPackLU);
```

An object for the Newton-Raphson method or the trust region method can then be defined as

```cpp
#if solverType == 0
  auto nr = Ikarus::makeNewtonRaphson(nonLinOp.subOperator<1, 2>(), std::move(linSolver));
#endif
#if solverType == 1
  auto nr = Ikarus::makeTrustRegion(nonLinOp);
  nr->setup({.verbosity = 1,
             .maxiter   = 30,
             .grad_tol  = 1e-8,
             .corr_tol  = 1e-8,
             .useRand   = false,
             .rho_reg   = 1e6,
             .Delta0    = 1});
#endif
```

All the available output messages are subscribed to be displayed by using the following commands:

```cpp
auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
nr->subscribeAll(nonLinearSolverObserver);
```

The [load control](../01_framework/controlRoutines.md#load-control) method is finally used as the path-following technique to solve this
non-linear problem.
It also subscribes to all the available information being written to the `vtkWriter`.
Output files are written for the deformed configuration at every load step that can be visualized using Paraview.
The load control method is executed by the following commands:

```cpp
auto lc = Ikarus::LoadControl(nr, 20, {0, 2000});
lc.subscribeAll(vtkWriter);
lc.run();
```

For postprocessing purposes we now write our results in a different VTK File. First we take a look at the stresses, in this case
the second Piola-Kirchhoff stress tensor. This quantity can be computed with the elements' `calculateAt()` function.
But here we will be using a so-called `ResultFunction`. This is a helper function that gathers the results over the
whole grid and can be used to generate data for a `VTKWriter`. To create this function, we can use the following helper

```cpp
auto stressFunction = Ikarus::makeResultFunction<ResultType::PK2Stress>(assembler);
```

To add this as vertex data to the VTK file we can do the following:

```cpp
Dune::VTKWriter resultWriter(gridView);
resultWriter.addVertexData(stressFunction);
```

As this functionality only writes out the results of the `calculateAt()` function as is, we can use `ResultEvaluators`
to further process our results. For example, we can use `ResultEvaluators::VonMises` to compute the Von Mises stress:

```cpp
  auto vonMisesFunction
      = Ikarus::makeResultFunction<ResultType::PK2Stress>(assembler, ResultEvaluators::VonMises{});
```

## Takeaways

- Grid types, finite element discretizations, and solver types are independent entities that are used to solve the problem at hand and
  can be switched easily to compare various formulations.
- The dune-python interface can be used to read external codes written in [Python](https://www.python.org/).
- A geometrically non-linear elastic finite element can be used from the Ikarus library.
- The Newton-Raphson and trust regions methods can be used as non-linear solvers.
- The load control method is used here as the path-following technique.
- Stress results can be written to a VTK file with a `ResultFunction`.

---
status: new
---

# Control routines

## Load control

A load control object is constructed as follows:

```cpp
auto lc = Ikarus::LoadControl(nonlinearSolver, numLoadSteps, {loadFactorStartValue, loadFactorEndValue});
```

- `nonlinearSolver` is a nonlinear solver, e.g., Newton-Raphson method, trust-region method, etc.
- `numLoadSteps` is the number of load steps.
- `loadFactorStartValue` is the value of the load factor at the beginning of the simulation.
- `loadFactorEndValue` is the value of the load factor at the end of the simulation.

The load control is started with the ``run()`` method, i.e., for the above-mentioned example:

```cpp
lc.run();
```

## Obtaining information from control routines

The load control is an observable object, i.e. one can subscribe to the messages of the load control method.
To read further on the implementation of observer patterns in Ikarus, see [here](observer.md).

The following messages are available:

```cpp
enum class ControlMessages {
  BEGIN,
  CONTROL_STARTED,
  CONTROL_ENDED,
  STEP_STARTED,
  STEP_ENDED,
  SOLUTION_CHANGED,
  END
};
```

## Path-following techniques

A general routine based on the standard arc-length method is included, which uses a scalar subsidiary function to impose
a constraint on the non-linear system of equations. The previously mentioned LoadControl method can also be recreated
using this technique. For more details on the standard arc-length method, see, among others, the works
of Wempner[@WEMPNER19711581], Crisfield[@CRISFIELD198155], Ramm[@wunderlich_strategies_1981] and
Riks[@riks_application_1972] among others. A path-following object is constructed as follows:

```cpp
auto alc = Ikarus::PathFollowing(nr, load_steps, stepSize, pft);
```

where `#!cpp nr` is a Newton-Raphson solver which considers a scalar subsidiary function and is defined by

```cpp
auto nr = Ikarus::makeNewtonRaphsonWithSubsidiaryFunction(nonLinOp, std::move(linSolver));
```

and `#!cpp pft` is the desired path-following technique. Three different path-following techniques are included, namely

- Standard arc-length method
- Load control method (as a subsidiary function under this generalized implementation)
- Displacement control method (uses a vector of indices which are all controlled by a same `#!cpp stepSize`).

These can be invoked by defining

```cpp
auto pft = Ikarus::StandardArcLength{};
auto pft = Ikarus::LoadControlWithSubsidiaryFunction{};
auto pft = Ikarus::DisplacementControl{controlledIndices};
```

!!! note
    The default path-following type is the `#!cpp Ikarus::StandardArcLength{}`.
    In the current implementation, it is assumed that the external forces are given by
    $F_{ext} = F_{ext}^0\lambda$ such that
    $$
    -\frac{\partial \mathbf{R}}{\partial \lambda} = F_{ext}^0
    $$
    An implementation for a general non-linear $F_{ext} = F_{ext}^0\left(\mathbf{D},\lambda\right)$ is an [open task](../03_contribution/openTask.md#control-routines---addons).

In order to create an own implementation for the scalar subsidiary function, the user has to create a `#!cpp struct`
with the following three member functions:

```cpp
void evaluateSubsidiaryFunction(SubsidiaryArgs& args) const;
void initialPrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args);
void intermediatePrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args);
```

For each Newton-Raphson iteration, the function `#!cpp evaluateSubsidiaryFunction(SubsidiaryArgs& args)` is used to evaluate the
subsidiary function and
its derivatives with respect to the displacement $\mathbf{D}$ and the load factor $\lambda$. The other two functions
are used to specify a prediction for $\mathbf{D}$ and $\lambda$ for the initial step and for
all the other intermediate subsequent `#!cpp load_steps`, respectively.

`#!cpp SubsidiaryArgs` is a `#!cpp struct` which is defined as

```cpp
struct SubsidiaryArgs {
  double stepSize; // (1)!
  Eigen::VectorX<double> DD; // (2)!
  double Dlambda{}; // (3)!
  double f{}; // (4)!
  Eigen::VectorX<double> dfdDD; // (5)!
  double dfdDlambda{}; // (6)!
  int currentStep; // (7)!
};
```

1. User-desired step size
2. Vector of displacement increments
3. Increment in the load factor
4. Scalar value evaluated from the subsidiary function
5. Derivative of the subsidiary function with respect to the displacement increment
6. Derivative of the subsidiary function with respect to the load factor increment
7. Current load step number

An example for the standard arc-length method is shown below:

```cpp
struct StandardArcLength {
    void evaluateSubsidiaryFunction(SubsidiaryArgs& args) const {
      if (psi) {
        const auto root = sqrt(args.DD.squaredNorm() + psi.value() * psi.value() * args.Dlambda * args.Dlambda);
        args.f          = root - args.stepSize;
        args.dfdDD      = args.DD / root;
        args.dfdDlambda = (psi.value() * psi.value() * args.Dlambda) / root;
      } else
        DUNE_THROW(Dune::InvalidStateException,
                   "You have to call initialPrediction first. Otherwise psi is not defined");
    }

    template <typename NonLinearOperator>
    void initialPrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      auto linearSolver
          = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);  // for the linear predictor step

      nonLinearOperator.lastParameter() = 1.0;  // lambda =1.0

      nonLinearOperator.template update<0>();
      const auto& R = nonLinearOperator.value();
      const auto& K = nonLinearOperator.derivative();

      linearSolver.factorize(K);
      linearSolver.solve(args.DD, -R);

      const auto DD2 = args.DD.squaredNorm();

      psi    = sqrt(DD2);
      auto s = sqrt(psi.value() * psi.value() + DD2);

      args.DD      = args.DD * args.stepSize / s;
      args.Dlambda = args.stepSize / s;

      nonLinearOperator.firstParameter() = args.DD;
      nonLinearOperator.lastParameter()  = args.Dlambda;
    }

    template <typename NonLinearOperator>
    void intermediatePrediction(NonLinearOperator& nonLinearOperator, SubsidiaryArgs& args) {
      nonLinearOperator.firstParameter() += args.DD;
      nonLinearOperator.lastParameter() += args.Dlambda;
    }

    std::string name = "Arc length";

  private:
    std::optional<double> psi;
  };
```

## Adaptive step-sizing for the path-following techniques

### Interface

The general interface for adaptive step-sizing is represented by the following concept.

```cpp
namespace Ikarus::Concepts {
    template <typename AdaptiveStepSizing, typename NonLinearSolverInformation, typename SubsidiaryArgs,
              typename NonLinearOperator>
    concept AdaptiveStepSizingStrategy = requires(AdaptiveStepSizing adaptiveSS, NonLinearSolverInformation info,
                                                  SubsidiaryArgs args, NonLinearOperator nop) {
      { adaptiveSS(info, args, nop) } -> std::same_as<void>; // (1)!
      { adaptiveSS.targetIterations() } -> std::same_as<int>; // (2)!
      { adaptiveSS.setTargetIterations(std::declval<int>()) } -> std::same_as<void>; // (3)!
    };
}
```

1. `#!cpp operator()` is overloaded such that the step size is modified.
2. Function that returns `#!cpp targetIterations` to be achieved. See [Iteration-based](#iteration-based), for example.
3. Function used to set `#!cpp targetIterations`.

For implementation details, refer to `ikarus/controlroutines/adaptivestepsizing.hh`.

### Implementations

#### No Operation

By default, `AdaptiveStepSizing::NoOp` is used with a path-following technique.
`AdaptiveStepSizing::NoOp` uses the step size provided by the user and doesn't modify them while using `PathFollowing`.
`NoOp` here stands for [No Operation](https://en.wikipedia.org/wiki/NOP_(code)).

#### Iteration-based

Instead of using a constant step size, the step size can be automatically adapted for efficient computations.
The `AdaptiveStepSizing::IterationBased` is implemented according to Ramm[@wunderlich_strategies_1981].
The step size can be scaled as shown below:

$$
\hat{s}_{k+1} = \left(\frac{\hat{i}}{i_k}\right)^{1/2}\hat{s}_k
$$

Here, $\hat{s}_{k+1}$ and $\hat{s}_{k}$ are the step sizes at $k+1$-th and $k$-th iteration.
Here, $\hat{i}$ is the desired number of iterations and $i_k$ is the number of iterations used in the previous step.

# Newton-Raphson method

## Description

`iks005_newtonRaphson.cpp` shows a basic example of the [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton's_method) to solve
a non-linear set of equations.
A function that shows the algorithm explicitly is provided, and another function which is implemented in Ikarus is
demonstrated. The function that depicts the Ikarus implementation uses a
[non-linear operator](../01_framework/nonlinearOperator.md) to
perform the Newton-Raphson iterations. A logger can also be subscribed to in order to observe the residual norms,
for instance.

## Code highlights

Here, the `#!cpp main()` function uses two functions, namely `#!cpp void newtonRaphsonVeryBasicExample();` and `#!cpp void newtonRaphsonBasicExampleWithLogger();`
that demonstrate the implementation of the Newton-Raphson scheme and also show the method to subscribe to loggers for information, respectively.
The function, which is solved in this example, and its derivative are mentioned below:

```cpp
auto f(double &x) { return 0.5 * x * x + x - 2; }
auto df(double &x) { return x + 1; }
```

First, the `#!cpp void newtonRaphsonVeryBasicExample();` function is described.
The settings for the Newton-Raphson method are defined as:

```cpp
double x               = 13; // (1)!
const double eps       = 1e-10; // (2)!
const int maxIter      = 20; // (3)!
const double xExpected = std::sqrt(5.0) - 1.0; // (4)!
```

1. Starting point for the Newton-Raphson method to find the closest root
2. Tolerance level below which the iterations stop
3. Maximum number of iterations
4. Expected value (analytical solution)

Functors are created then to evaluate the function to be solved for and its derivative. These are then passed to the non-linear operator
as shown below:

```cpp
auto fvLambda  = [&](auto &&x) { return f(x); };
auto dfvLambda = [&](auto &&x) { return df(x); };
Ikarus::NonLinearOperator nonLinOp(Ikarus::functions(fvLambda, dfvLambda), x);
```

The standard implementation of the Newton-Raphson method is illustrated in this function, which also uses `nonLinOp`.

```cpp
int iterCount = 1;
typename Ikarus::NonLinearOperator::Domain x;
auto f = nonLinOp(x);
auto df = derivative(nonLinOp)(x);
while (abs(nonLinOp.value()) > eps and iterCount <= maxIter) {

  x -=f / df;
  f = nonLinOp(x);
  df = derivative(nonLinOp)(x);
  iterCount++;

  std::cout << "nonlinearOperator, value(): " << nonLinOp(x) << "\n";
  std::cout << "nonlinearOperator, x: " << x << "\n";
}
```

One could also use the existing functionality in Ikarus to obtain a similar solution from the Newton-Raphson scheme, as depicted below:

```cpp
Ikarus::NewtonRaphson nr(nonLinOp);
nr.setup({eps, maxIter});
const auto solverInfo = nr.solve(x);

std::cout << "success: " << solverInfo.success << "\n";
std::cout << "iterations: " << solverInfo.iterations << "\n";
std::cout << "residuum: " << solverInfo.residualNorm << "\n";
std::cout << "solution: " << x << "\n";
std::cout << "expected solution: " << xExpected << "\n";
```

Further details on the non-linear solver can be found [here](../01_framework/solvers.md#nonlinear-solver).

Instead of using the multiple `#!cpp std::cout` statements, one can simply subscribe for the desired information
using the functionalities from the [observer module](../01_framework/observer.md).
This is exemplified by the function `#!cpp void newtonRaphsonBasicExampleWithLogger();` and the `#!cpp class OurFirstObserver`.
It is also possible to subscribe to the existing non-linear solver messages mentioned [here](../01_framework/observer.md#messages).

```cpp
class OurFirstObserver : public IObserver<NonLinearSolverMessages> {
 public:
  void updateImpl(NonLinearSolverMessages message) override {
    if (message == NonLinearSolverMessages::ITERATION_STARTED) std::cout << "Iteration started.\n";
  }
};

void newtonRaphsonBasicExampleWithLogger() {
  double x               = 13;
  const double eps       = 1e-10;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;

  auto fvLambda  = [&](auto &&x) { return f(x); };
  auto dfvLambda = [&](auto &&x) { return df(x); };
  Ikarus::NonLinearOperator nonLinOp(Ikarus::functions(fvLambda, dfvLambda), Ikarus::parameter(x));

  Ikarus::NewtonRaphson nr(nonLinOp);
  nr.setup({eps, maxIter});

  auto ourSimpleObserver = std::make_shared<OurFirstObserver>();
  nr.subscribe(NonLinearSolverMessages::ITERATION_STARTED, ourSimpleObserver);

  const auto solverInfo = nr.solve(x);
  if (solverInfo.success)
    std::cout << "solution: " << x << "\n";
  else
    std::cout << "The Newton-Raphson procedure failed to converge" << std::endl;
}
```

## Takeaways

- Functors for the function and its derivative can be used to create a simple non-linear operator.
- A `NewtonRaphson` object can be created using the non-linear operator.
- The settings for the Newton-Raphson scheme can be modified by using the `#!cpp setup()` function.
- Nonlinear solver messages can be subscribed to print the desired quantities.

# Nonlinear operator

The class ``NonLinearOperator`` consists of a collection of a function and its derivatives, including their dependence
on parameters.

Let us assume a function `f(x)` and its derivative `df(x)`.
Then, a ``NonLinearOperator`` can be constructed as follows:

```cpp
double x               = 13;
auto fvLambda  = [&](auto&& x) { return f(x); };
auto dfvLambda = [&](auto&& x) { return df(x); };

auto nonLinOp = Ikarus::NonLinearOperator(functions(fvLambda, dfvLambda), x);
```

!!! note
    It is assumed that the second function is the derivative of the first function, the third function is the derivative of the second
    function (2nd derivative of the first function), and so on.

``functions(...)`` and ``parameter(...)`` are helper functions. They are necessary to distinguish which argument is a function and which
argument is a parameter.

``nonLinOp`` provides the following features:

```cpp
void updateAll() // (1)!
void update<n>() // (2)!
auto& value() // (3)!
auto& derivative() // (4)!
auto& secondDerivative() // (5)!
auto& nthDerivative<n>() // (6)!
auto& firstParameter() // (7)!
auto& secondParameter() // (8)!
auto& nthParameter<n>() // (9)!
auto& lastParameter() // (10)!
auto subOperator<n,m,...>() // (11)!
```

1. Evaluates all functions.
2. Evaluates the n-th function in ``functions(...)`` . Counting starts from 0, as always in C++.
3. Returns the result of the function evaluation.
4. Returns the result of the evaluation of the first derivative (if the function for the first derivative is passed to the nonlinear
   operator during construction).
5. Returns the result of the evaluation of the second derivative (if the function for the second derivative is passed to the nonlinear
   operator during construction).
6. Returns the result of the evaluation of the n-th derivative (if the function for the n-th derivative is passed to the nonlinear
   operator during construction).
7. Returns the value of the first parameter.
8. Returns the value of the second parameter (if available).
9. Returns the value of the n-th parameter (if available).
10. Returns the value of the last parameter.
11. Creates an `Ikarus::NonLinearOperator` with a subset of the derivatives. For example, let us consider a nonlinear operator with
    (function, first derivative, second derivative). ``subOperator<0,1>()`` then returns a nonlinear operator with
    (function, first derivative).

<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Finite elements

Several disciplines associate finite elements with different meanings.
In Ikarus, finite elements have two different objectives.
The first one is to provide an evaluation of the scalars, vectors, and matrices. 
These are associated with an algebraic representation of discrete energies, weak forms, or bilinear forms.
These algebraic objects are usually constructed using some combination of [local functions](localFunctions.md) and 
parameters stemming from the underlying physical problem, e.g., load factor, Young's modulus, or viscosity.

The second task of finite elements is to evaluate derived results in the element parameter space, e.g., stresses or geometric quantities.
This leads to the following interface for the finite elements:
## Interface
Local functions provide the following interface
```cpp
ScalarType evaluateScalar(const FErequirements& req);
void evaluateVector(const FErequirements& req, VectorType& b);
void evaluateMatrix(const FErequirements& req, MatrixType& A);
void calculateLocalSystem(const FErequirements& req, MatrixType& A, VectorType& b);
void calculateAt(const Resultrequirements& req, const Eigen::Vector<double, Traits::mydim>& local,
                     ResultTypeMap<ScalarType>& result);
void globalIndices(std::vector<GlobalIndex>& globalIndices);
```

Please refer to the [FE requirements](feRequirements.md) to learn more about the finite element requirements and result requirements. 
The first four methods receive an object of type `FErequirements`. This object is responsible for passing different types of information needed for the local evaluation of the local linear algebra objects.
The first method, `evaluateScalar`, simply returns by value because it is cheaper to return a `double`, for example when evaluating energy.
The other methods, `evaluateVector`, `evaluateMatrix`, and `calculateLocalSystem`, receive one or two additional output arguments where the results are to be written.
This interface is needed to circumvent the dynamic memory allocation, that is required if these methods return by value.

The method `calculateAt` is responsible for evaluating several results, and it receives the `ResultRequirements` object.
These results are stored inside the output argument `result`, which is of the type `ResultTypeMap`.
Additionally, there is the argument `local`, which contains the element coordinates where the results are to be evaluated.


A typical `calculateAt` method is implemented as shown below: 

```cpp
typename ResultTypeMap<double>::ResultArray res;
if(req.isResultRequested( ResultType::gradientNormOfMagnetization)) {
  res.resize(1,1);
  res(0,0)=...;
  result.insertOrAssignResult(ResultType::gradientNormOfMagnetization,res);
}
if(req.isResultRequested( ResultType::BField)) {
  res.setZero(3,1);
  res=...;
  result.insertOrAssignResult(ResultType::BField,res);
}
if(req.isResultRequested( ResultType::cauchyStress)) {
  res.setZero(3,3);
  res = ...;
  result.insertOrAssignResult(ResultType::cauchyStress,res);
}
```
!!! note "`ResultTypeMap<double>::ResultArray`"
    `#!cpp ResultTypeMap<double>::ResultArray` is an object of type `#!cpp Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,3,3>`.
    Thus, the maximum size of `result` is limited to a 3x3 matrix. This is used to circumvent dynamic memory allocations again.


The last method is the `globalIndices`. It is used to write a finite element's global indices to the output parameter `globalIndices`.
This information originates from a `basis` object. See existing implementations for details.

## Linear and Non-linear Elasticity
* To be added

## Enhanced Assumed Strain Elements
The Enhanced Assumed Strain (EAS) elements are a class of finite elements that helps to avoid the locking phenomenon.
They are obtained by re-parametrizing the Hu-Washizu principle and enforcing an orthogonality condition. 
This results in an extension of the standard pure displacement formulation with an enhanced strain field ($\tilde\epsilon$) 
as an additional independent variable.The locking characteristics of the pure displacement formulations can be eliminated with an appropriate choice of ansatz space for $\tilde\epsilon$. For further theoretical aspects, the readers are referred to [@simo_class_1990] 
and [@andelfinger_eas-elements_1993]. The EAS formulation is currently implemented for the linear-elastic case, but 
it could be extended to the non-linear regime. The currently implemented EAS elements are the following:

* Q1E4
* Q1E5
* Q1E7
* H1E9
* H1E21

The notation used here is described as follows: the first alphabet stands for a Quadrilateral (Q) or a Hexahedral (H) element.
The second index denotes the order of the element. E stands for the EAS element, and the number following that denotes the 
number of EAS parameters used to enhance the strain field. The only difference amongst various EAS formulations arises 
from the matrix, $\mathbf{M}$ which is used to approximate the enhanced strain field. An example for the calculation of the 
matrix $\mathbf{M}$ for a Q1E4 element is shown below:
```cpp
template <typename Geometry>
struct EASQ1E4 {
  static constexpr int strainSize         = 3;
  static constexpr int enhancedStrainSize = 4;

  EASQ1E4() = default;
  explicit EASQ1E4(const Geometry& geometry)
      : geometry{std::make_unique<Geometry>(geometry)}, T0InverseTransformed{calcTransformationMatrix2D(geometry)} {}

  auto calcM(const Dune::FieldVector<double, 2>& quadPos) const {
    Eigen::Matrix<double, strainSize, enhancedStrainSize> M;
    M.setZero(strainSize, enhancedStrainSize);
    const double xi   = quadPos[0];
    const double eta  = quadPos[1];
    M(0, 0)           = 2 * xi - 1.0;
    M(1, 1)           = 2 * eta - 1.0;
    M(2, 2)           = 2 * xi - 1.0;
    M(2, 3)           = 2 * eta - 1.0;
    const double detJ = geometry->integrationElement(quadPos);
    M                 = T0InverseTransformed / detJ * M;
    return M;
  }

  std::unique_ptr<Geometry> geometry;
  Eigen::Matrix3d T0InverseTransformed;
};
```
It is to be noted that the ansatz spaces for the matrix $\mathbf{M}$ are to be modified such that they fulfill the orthogonality 
condition in the $\left[0,1\right]$ element domain used in DUNE, in contrast to the $\left[-1,1\right]$ usually found in 
literature.

In order to add a new EAS element, the following additions are to be made: 

1. Create a `#!cpp struct` to calculate the matrix $\mathbf{M}$ as shown above exemplarily for the Q1E4 element.
2. Add the new variant to the corresponding list of 2D and 3D variants as shown below:
```cpp
template <typename Geometry>
using EAS2DVariant = std::variant<EASQ1E4<Geometry>, EASQ1E5<Geometry>, EASQ1E7<Geometry>>;
template <typename Geometry>
using EAS3DVariant = std::variant<EASH1E9<Geometry>, EASH1E21<Geometry>>;
```
3. Finally, add the new EAS variant with an appropriate switch statement (as shown below) to automatically call the 
desired functions
```cpp
void setEASType(int numberOfEASParameters) {
    if constexpr (Traits::mydim == 2) {
      switch (numberOfEASParameters) {
        case 0:
          onlyDisplacementBase = true;
          break;
        case 4:
          easVariant = EASQ1E4(DisplacementBasedElement::getLocalView().element().geometry());
          break;
        case 5:
          easVariant = EASQ1E5(DisplacementBasedElement::getLocalView().element().geometry());
          break;
        case 7:
          easVariant = EASQ1E7(DisplacementBasedElement::getLocalView().element().geometry());
          break;
        default:
          DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available for the 2D case.");
          break;
      }
    } else if constexpr (Traits::mydim == 3) {
      switch (numberOfEASParameters) {
        case 0:
          onlyDisplacementBase = true;
          break;
        case 9:
          easVariant = EASH1E9(DisplacementBasedElement::getLocalView().element().geometry());
          break;
        case 21:
          easVariant = EASH1E21(DisplacementBasedElement::getLocalView().element().geometry());
          break;
        default:
          DUNE_THROW(Dune::NotImplemented, "The given EAS parameters are not available for the 3D case.");
          break;
      }
    }
  }
```

If the number of EAS parameters is set to zero, the pure displacement formulation is then utilised for analysis.

\bibliography
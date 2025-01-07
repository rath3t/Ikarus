// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file resultfunction.hh
 * \brief Ikarus Result Function for Stress and other finite element results
 * \ingroup io
 *
 */

#pragma once

#include <type_traits>

#include <dune/vtk/function.hh>
#include <dune/vtk/vtkwriter.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/concepts.hh>

namespace Ikarus {
namespace Impl {
  struct DefaultUserFunction
  {
  };
} // namespace Impl

/**
 * \brief Wrapper to evaluate results for a vtkwriter.
 * \details
 * Usage:
 * \code
 * // Usage with Dune::Vtk::VtkWriter
 * auto resultFunction = Ikarus::makeResultVtkFunction<resType>(assembler);
 * vtkwriter.addPointData(resultFunction);
 *
 * // Usage with the native Dune::VTKWriter
 * auto resultFunction = Ikarus::makeResultFunction<resType>(assembler);
 * vtkWriter.addVertexData(resultFunction);
 * \endcode
 * \ingroup io
 * \relates
 * \tparam AS underlying assembler (provides the finite elements and the requested results)
 * \tparam RT requested result type
 * \tparam UserFunction Type of the user-defined function for custom result evaluation (default is
DefaultUserFunction)
 */
template <typename AS, template <typename, int, int> class RT, typename UserFunction = Impl::DefaultUserFunction>
requires(Concepts::FlatAssembler<AS> and Concepts::ResultType<RT>)
class ResultFunction : public Dune::VTKFunction<typename AS::GridView>
{
public:
  using Assembler         = AS;
  using GridView          = typename Assembler::GridView;
  using FERequirementType = typename Assembler::FERequirement;
  using FEContainer       = typename Assembler::FEContainer;
  using FiniteElement     = typename std::remove_cvref_t<FEContainer>::value_type;

  using ctype                  = typename GridView::ctype;
  constexpr static int griddim = GridView::dimension;
  using Entity                 = typename GridView::template Codim<0>::Entity;

  /**
   * \brief Evaluate the component at a given entity and local coordinates.
   * \details
   * This function is required by the Dune::VTKFunction interface.
   *
   * \param comp Stress component index
   * \param e Entity on which to evaluate the stress
   * \param local Local coordinates within the entity
   * \return Stress component value
   */
  double evaluate(int comp, const Entity& e, const Dune::FieldVector<ctype, griddim>& local) const override {
    const auto index = gridView().indexSet().index(e);
    return evaluateComponent(index, local, comp);
  }

  /**
   * \brief Get the number of components.
   * \details
   * This function is required by the Dune::VTKFunction interface.
   *
   * \return Number of stress components
   */
  [[nodiscard]] int ncomps() const override {
    if constexpr (std::is_same_v<UserFunction, Impl::DefaultUserFunction>) {
      Dune::FieldVector<ctype, griddim> val(0.0);
      auto sigma = finiteElements().at(0).template calculateAt<RT>(requirement(), val).asVec();
      return static_cast<int>(sigma.size());
    } else
      return userFunction_.ncomps();
  }

  /**
   * \brief Get the name of the result type.
   * \details
   * This function is required by the Dune::VTKFunction interface.
   *
   * \return String representing the name of the result type
   */
  [[nodiscard]] constexpr std::string name() const override {
    if constexpr (requires { userFunction_.name(); })
      return userFunction_.name();
    else
      return toString<RT>();
  }

  /**
   * \brief Get the precision used for this result
   * \details
   * This function is part of the Dune::VTKFunction interface.
   * This has no affect when the ResultFunction is used as part of a the dune-vtk module
   *
   * \return Precision (i.e. float64 or float32)
   */
  Dune::VTK::Precision precision() const override { return prec_; }

  /**
   * \brief Constructor for ResultFunction.
   * \details
   * Constructs a ResultFunction object with given finite elements, ferequirements
   *
   * \param assembler shared pointer to the underlying assembler (provides the finite elements and the requested
   * results)
   * \param prec (optional) specify the used precision (only has an effect when using resultfunciton with
   * Dune::VTK::VTKWriter)
   */
  ResultFunction(std::shared_ptr<Assembler> assembler, Dune::VTK::Precision prec = Dune::VTK::Precision::float64)
      : assembler_(assembler),
        prec_{prec},
        userFunction_{UserFunction{}} {}

private:
  double evaluateComponent(int eleID, const Dune::FieldVector<ctype, griddim>& local, int comp) const {
    auto result = finiteElements().at(eleID).template calculateAt<RT>(requirement(), local).asVec();

    if constexpr (!std::is_same_v<UserFunction, Impl::DefaultUserFunction>)
      return userFunction_(result, comp);
    else
      return result(comp);
  }

  const FEContainer& finiteElements() const { return assembler_->finiteElements(); }
  const FERequirementType& requirement() const { return assembler_->requirement(); }
  const GridView& gridView() const { return assembler_->gridView(); }

  std::shared_ptr<Assembler> assembler_;

  Dune::VTK::Precision prec_;
  [[no_unique_address]] std::string name_{};
  UserFunction userFunction_;
};

/**
 * \brief Function to create a ResultFunction as a shared_ptr
 * \details
 * Constructs a ResultFunction object with given finite elements, ferequirements as shared_ptr to be used with
 * the native Dune VTKWriter
 *
 * \param assembler shared pointer to the underlying assembler (provides the finite elements and the requested results)
 * \tparam AS type of the assembler
 * \param prec (optional) specify the used precision
 * \tparam RT requested result type
 * \tparam UserFunction Type of the user-defined function for custom result evaluation (default is DefaultUserFunction)
 */
template <template <typename, int, int> class RT, typename UserFunction = Impl::DefaultUserFunction, typename AS>
auto makeResultFunction(std::shared_ptr<AS> assembler, Dune::VTK::Precision prec = Dune::VTK::Precision::float64) {
  return std::make_shared<ResultFunction<AS, RT, UserFunction>>(assembler, prec);
}

/**
 * \brief Function to create a ResultFunction as a gridfunction that can be used with dune-vtk
 * \details
 * Constructs a ResultFunction object with given finite elements, ferequirements as a VTK::Function to be used with
 * dune-vtk It is possible to construct a localFunction from this as follows
 * \code
 * auto localResultFunction = localFunction(vtkResultFunction);
 * localResultFunction.bind(element);
 * \endcode
 * \param assembler shared pointer to the underlying assembler (provides the finite elements and the requested results)
 * \tparam AS  type of the assembler
 * \tparam RT requested result type
 * \tparam UserFunction Type of the user-defined function for custom result evaluation (default is
 * DefaultUserFunction)
 */
template <template <typename, int, int> class RT, typename UserFunction = Impl::DefaultUserFunction, typename AS>
auto makeResultVtkFunction(std::shared_ptr<AS> assembler) {
  return Dune::Vtk::Function<typename AS::GridView>(std::make_shared<ResultFunction<AS, RT, UserFunction>>(assembler));
}

} // namespace Ikarus

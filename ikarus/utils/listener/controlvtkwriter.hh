// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlvtkwriter.hh
 * \brief Observer implementation for writing vtk files when notified
 */

#pragma once
#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <Eigen/Core>

#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/listener/listener.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

/**
 * \brief ControlSubsamplingVertexVTKWriter class for writing VTK files with subsampling based on control messages.
 *
 * \details It inherits from the Listener class and is specifically designed for handling SOLUTION_CHANGED messages.
 *
 * \tparam B The type of the grid basis.
 */
template <typename B>
class ControlSubsamplingVertexVTKWriter : public Listener
{
  using Basis                     = B;
  static constexpr int components = Basis::LocalView::Tree::degree() == 0 ? 1 : Basis::LocalView::Tree::degree();

public:
  /**
   * \brief Constructor for ControlSubsamplingVertexVTKWriter.
   *
   * Initializes the VTK writer with the provided basis, solution, and refinement levels.
   *
   * \param basis The grid basis.
   * \param sol The solution vector.
   * \param refinementLevels The refinement levels for subsampling.
   */
  ControlSubsamplingVertexVTKWriter(const Basis& basis, const Eigen::VectorXd& sol, int refinementLevels = 0)
      : basis_{&basis},
        vtkWriter_(basis.gridView(), Dune::refinementLevels(refinementLevels)),
        solution_{&sol} {}

  template <typename BC>
  ControlSubsamplingVertexVTKWriter& subscribeTo(BC& bc) {
    this->subscribe(bc, [&](ControlMessages message) { this->updateImpl(message); });
    return *this;
  }

  /**
   * \brief Set field information for the VTK file.
   *
   * \param name The name of the field.
   * \param type The type of the field.
   * \param size The size of the field.
   * \param prec The precision of the field.
   * \return The field information.
   */
  auto setFieldInfo(std::string&& name, Dune::VTK::FieldInfo::Type type, std::size_t size,
                    Dune::VTK::Precision prec = Dune::VTK::Precision::float32) {
    fieldInfo_      = Dune::VTK::FieldInfo(std::move(name), type, size, prec);
    isFieldInfoSet_ = true;
  }

  /**
   * \brief Set the file name prefix for VTK files.
   *
   * \param p_name The file name prefix.
   */
  auto setFileNamePrefix(std::string&& name) { prefixString_ = std::move(name); }

  /**
   * \brief Implementation of the update method.
   *
   * This method is called upon receiving a SOLUTION_CHANGED control message.
   * It writes VTK files with subsampling based on the provided field information.
   *
   * \param message The received control message.
   */
  void updateImpl(ControlMessages message) {
    assert(isFieldInfoSet_ && "You need to call setFieldInfo first!");
    switch (message) {
      case ControlMessages::SOLUTION_CHANGED: {
        auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, components>>(*basis_,
                                                                                                            *solution_);
        vtkWriter_.addVertexData(disp, fieldInfo_);
        vtkWriter_.write(prefixString_ + std::to_string(step_++));
      } break;
      default:
        break; // default: do nothing when notified
    }
  }

private:
  const Basis* basis_;
  Dune::SubsamplingVTKWriter<typename Basis::GridView> vtkWriter_;
  const Eigen::VectorXd* solution_;
  int step_{0};
  Dune::VTK::FieldInfo fieldInfo_{"Default", Dune::VTK::FieldInfo::Type::scalar, 1};
  std::string prefixString_{};
  bool isFieldInfoSet_{false};
};
} // namespace Ikarus
#pragma GCC diagnostic pop

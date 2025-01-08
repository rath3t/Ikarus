// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "../pythonhelpers.hh"

#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>

#include <spdlog/spdlog.h>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/python/finiteelements/material.hh>

void addMaterialsSubModule(pybind11::module& m) {
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;
  using namespace Eigen;

  auto materials = m.def_submodule("materials", "This is the submodule for materials in Ikarus");

  ENUM_BINDINGS_WITH_MODULE(StrainTags, materials);
  ENUM_BINDINGS_WITH_MODULE(StressTags, materials);
  ENUM_BINDINGS_WITH_MODULE(TangentModuliTags, materials);

  pybind11::class_<Materials::LinearElasticity> linElastic(materials, "LinearElasticity");
  Ikarus::Python::registerLinearElasticity(materials, linElastic);

  pybind11::class_<Materials::StVenantKirchhoff> svk(materials, "StVenantKirchhoff");
  Ikarus::Python::registerStVenantKirchhoff(materials, svk);

  pybind11::class_<Materials::NeoHooke> nh(materials, "NeoHooke");
  Ikarus::Python::registerNeoHooke(materials, nh);

  /**
   * \brief Transform strain from one type to another.
   *
   * This function transforms one strain component matrix from one type to another, based on the provided strain tags
   * \ingroup  materials
   * \param from Type of the source strain tag.
   * \param to Type of the target strain tag.
   * \param E Eigen matrix representing the input strain (can be in Voigt notation).
   * \return The transformed strain matrix.
   */
  materials.def(
      "transformStrain",
      [](StrainTags from, StrainTags to, Eigen::MatrixXd E) -> Eigen::MatrixXd {
        auto callTransformStrain =
            []<StrainTags from_, StrainTags to_>(const Eigen::MatrixXd& eRaw) -> Eigen::MatrixXd {
          if (eRaw.cols() == 1) {
            Eigen::Vector<double, 6> E = eRaw;
            return transformStrain<from_, to_>(E);
          } else {
            Eigen::Matrix<double, 3, 3> E = eRaw;
            return transformStrain<from_, to_>(E);
          }
        };
        if (!((E.rows() == 3 && E.cols() == 3) || (E.rows() == 6 && E.cols() == 1)))
          DUNE_THROW(Dune::IOError,
                     "Strain conversions are only implemented for matrices of dimension 3 or Voigt vectors of size 6");

        if (from == StrainTags::linear) {
          spdlog::warn("No useful transformation available for linear strains");
          return E;
        }
        if (from == to)
          return E;

        if (to == StrainTags::greenLagrangian) {
          switch (from) {
            case StrainTags::deformationGradient:
              return callTransformStrain
                  .template operator()<StrainTags::deformationGradient, StrainTags::greenLagrangian>(E);
            case Ikarus::StrainTags::displacementGradient:
              return callTransformStrain
                  .template operator()<StrainTags::displacementGradient, StrainTags::greenLagrangian>(E);
            case Ikarus::StrainTags::rightCauchyGreenTensor:
              return callTransformStrain
                  .template operator()<StrainTags::rightCauchyGreenTensor, StrainTags::greenLagrangian>(E);
            default:
              __builtin_unreachable();
          }
        } else if (to == StrainTags::deformationGradient) {
          switch (from) {
            case StrainTags::greenLagrangian:
              return callTransformStrain
                  .template operator()<StrainTags::greenLagrangian, StrainTags::deformationGradient>(E);
            case Ikarus::StrainTags::displacementGradient:
              return callTransformStrain
                  .template operator()<StrainTags::displacementGradient, StrainTags::deformationGradient>(E);
            case Ikarus::StrainTags::rightCauchyGreenTensor:
              return callTransformStrain
                  .template operator()<StrainTags::rightCauchyGreenTensor, StrainTags::deformationGradient>(E);
            default:
              __builtin_unreachable();
          }
        } else if (to == StrainTags::rightCauchyGreenTensor) {
          switch (from) {
            case StrainTags::greenLagrangian:
              return callTransformStrain
                  .template operator()<StrainTags::greenLagrangian, StrainTags::rightCauchyGreenTensor>(E);
            case Ikarus::StrainTags::displacementGradient:
              return callTransformStrain
                  .template operator()<StrainTags::displacementGradient, StrainTags::rightCauchyGreenTensor>(E);
            case Ikarus::StrainTags::deformationGradient:
              return callTransformStrain
                  .template operator()<StrainTags::deformationGradient, StrainTags::rightCauchyGreenTensor>(E);
            default:
              __builtin_unreachable();
          }
        }
        __builtin_unreachable();
      },
      py::arg("to"), py::arg("from"), py::arg("strain"));
}
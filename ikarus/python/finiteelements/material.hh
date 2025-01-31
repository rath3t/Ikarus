// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file material.hh
 * \brief Python bindings for materials
 */

#pragma once

#include <dune/common/classname.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <ikarus/finiteelements/mechanics/materials.hh>
#if ENABLE_MUESLI
  #include <ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh>
#endif
#include <ikarus/utils/concepts.hh>

#define MAKE_MaterialFunction(clsName, materialName, functionname, vecSize)                                    \
  clsName.def(                                                                                                 \
      #functionname,                                                                                           \
      [](materialName& self, StrainTags straintag, Eigen::Ref<const Eigen::Vector<double, vecSize>> eVoigt_) { \
        if constexpr (not(materialName::strainTag == StrainTags::linear)) {                                    \
          Eigen::Vector<double, vecSize> eVoigt = eVoigt_;                                                     \
          if (straintag == StrainTags::rightCauchyGreenTensor)                                                 \
            return self.template functionname<StrainTags::rightCauchyGreenTensor>(eVoigt);                     \
          else if (straintag == StrainTags::greenLagrangian)                                                   \
            return self.template functionname<StrainTags::greenLagrangian>(eVoigt);                            \
          else if (straintag == StrainTags::linear)                                                            \
            DUNE_THROW(Dune::MathError, "Passing linear strain to " + std::string(#materialName) +             \
                                            " does not makes sense use LinearElastic class");                  \
          else if (straintag == StrainTags::displacementGradient)                                              \
            DUNE_THROW(Dune::MathError,                                                                        \
                       "Passing displacementGradient strain in 6d Voigt notation does not make any sense!");   \
          else if (straintag == StrainTags::deformationGradient)                                               \
            DUNE_THROW(Dune::MathError,                                                                        \
                       "Passing deformationGradient strain in 6d Voigt notation does not make any sense!");    \
          else                                                                                                 \
            DUNE_THROW(Dune::MathError, toString(straintag) + "is not a valid strain tag.");                   \
        } else {                                                                                               \
          Eigen::Vector<double, vecSize> eVoigt = eVoigt_; /* Linear elastic path */                           \
          if (straintag == StrainTags::linear)                                                                 \
            return self.template functionname<StrainTags::linear>(eVoigt);                                     \
          else                                                                                                 \
            DUNE_THROW(Dune::MathError, "Linear elastic material only accepts linear strains!");               \
        }                                                                                                      \
        __builtin_unreachable();                                                                               \
      },                                                                                                       \
      "StrainName"_a, "strainVector"_a);

namespace Ikarus::Python {

namespace Impl {
  // Function to extract and convert parameters using a conversion strategy
  template <typename T>
  LamesFirstParameterAndShearModulus convertMaterialParameters(const pybind11::kwargs& kwargs,
                                                               const std::string& param1, const std::string& param2) {
    auto converter =
        convertLameConstants(T{kwargs[param1.c_str()].cast<double>(), kwargs[param2.c_str()].cast<double>()});

    double lambda = converter.toLamesFirstParameter();
    double mu     = converter.toShearModulus();
    return {lambda, mu};
  }

  Ikarus::LamesFirstParameterAndShearModulus extractMaterialParameters(const pybind11::kwargs& kwargs) {
    // clang-format off
    static const std::map<std::array<std::string, 2>, std::function<LamesFirstParameterAndShearModulus(const pybind11::kwargs&)>> conversionMap = {
      {{"E", "nu"}, [](const auto& kw){ return convertMaterialParameters<YoungsModulusAndPoissonsRatio>(kw, "E", "nu"); }},
      {{"E", "mu"}, [](const auto& kw){ return convertMaterialParameters<YoungsModulusAndShearModulus>(kw, "E", "mu"); }},
      {{"E", "K"}, [](const auto& kw){ return convertMaterialParameters<YoungsModulusAndBulkModulus>(kw, "E", "K"); }},
      {{"E", "Lambda"}, [](const auto& kw){ return convertMaterialParameters<YoungsModulusAndLamesFirstParameter>(kw, "E", "Lambda"); }},
      {{"K", "Lambda"}, [](const auto& kw){ return convertMaterialParameters<BulkModulusAndLamesFirstParameter>(kw, "K", "Lambda"); }},
      {{"Lambda", "mu"}, [](const auto& kw){ return LamesFirstParameterAndShearModulus{kw["Lambda"].template cast<double>(), kw["mu"].template cast<double>()}; }}
    };
    // clang-format on

    if (kwargs.size() != 2)
      DUNE_THROW(Dune::IOError, "The number of material parameters passed to the material should be 2");

    for (const auto& [materialParameters, parameterConverter] : conversionMap) {
      const auto [firstPar, secondPar] = materialParameters;
      if (kwargs.contains(firstPar) && kwargs.contains(secondPar)) {
        return parameterConverter(kwargs);
      }
    }

    DUNE_THROW(Dune::IOError,
               "No suitable combination of material parameters found, valid combinations are: (E, nu), (E, mu), (E, "
               "K), (E, Lambda), (K, Lambda), (Lambda, nu)");
  }
} // namespace Impl
template <class Material, size_t vecSize, bool registerConstructor, class... options>
void registerMaterial(pybind11::handle scope, pybind11::class_<Material, options...> cls) {
  using pybind11::operator""_a;
  namespace py = pybind11;

  if constexpr (registerConstructor) {
    cls.def(pybind11::init([](const py::kwargs& kwargs) {
      auto matParameter = Impl::extractMaterialParameters(kwargs);
      return new Material(matParameter);
    }));
  }

  cls.def_property_readonly("name", [](Material& self) { return self.name(); });

  std::string materialname = Material::name();

  MAKE_MaterialFunction(cls, Material, storedEnergy, vecSize);
  MAKE_MaterialFunction(cls, Material, stresses, vecSize);
  MAKE_MaterialFunction(cls, Material, tangentModuli, vecSize);

  using PlaneStressClass = decltype(Materials::planeStress(std::declval<Material>()));
  auto includes          = Dune::Python::IncludeFiles{"ikarus/finiteelements/mechanics/materials.hh"};
  auto pS                = Dune::Python::insertClass<PlaneStressClass>(
                scope, std::string("PlaneStress_") + materialname,
                Dune::Python::GenerateTypeName(
                    "Ikarus::Materials::VanishingStress<std::array<Ikarus::Materials::MatrixIndexPair, 3ul >"
                                   "{"
                                   "{Ikarus::Materials::MatrixIndexPair{2ul, 1ul}, Ikarus::Materials::MatrixIndexPair{2ul, 0ul},"
                                   "Ikarus::Materials::MatrixIndexPair{2ul, 2ul}}}," +
                    Dune::className<Material>() + ">"),
                includes)
                .first;
  MAKE_MaterialFunction(pS, PlaneStressClass, storedEnergy, 3);
  MAKE_MaterialFunction(pS, PlaneStressClass, stresses, 3);
  MAKE_MaterialFunction(pS, PlaneStressClass, tangentModuli, 3);
  MAKE_MaterialFunction(pS, PlaneStressClass, storedEnergy, 6);
  MAKE_MaterialFunction(pS, PlaneStressClass, stresses, 6);
  MAKE_MaterialFunction(pS, PlaneStressClass, tangentModuli, 6);

  cls.def("asPlaneStress", [](Material& self) {
    return Materials::planeStress(self);
  }); /* no keep_alive since planeStress copies the material */

  using PlaneStrainClass = decltype(Materials::planeStrain(std::declval<Material>()));
  auto pStrain           = Dune::Python::insertClass<PlaneStrainClass>(
                     scope, std::string("PlaneStrain_") + materialname,
                     Dune::Python::GenerateTypeName(
                         "Ikarus::Materials::VanishingStrain<std::array<Ikarus::Materials::MatrixIndexPair, "
                                   "3ul>{{Ikarus::Materials::MatrixIndexPair{2ul, 1ul},"
                                   "Ikarus::Materials::MatrixIndexPair{2ul,0ul}, Ikarus::Materials::MatrixIndexPair{"
                                   "2ul, 2ul}}}," +
                         Dune::className<Material>() + ">"),
                     includes)
                     .first;
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, storedEnergy, 3);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, stresses, 3);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, tangentModuli, 3);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, storedEnergy, 6);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, stresses, 6);
  MAKE_MaterialFunction(pStrain, PlaneStrainClass, tangentModuli, 6);

  cls.def("asPlaneStrain", [](Material& self) {
    return Materials::planeStrain(self);
  }); /* no keep_alive since planeStrain copies the material */
  using ShellMaterialClass = decltype(Materials::shellMaterial(std::declval<Material>()));
  auto shellmaterial       = Dune::Python::insertClass<ShellMaterialClass>(
                           scope, std::string("Shell_") + materialname,
                           Dune::Python::GenerateTypeName(
                               "Ikarus::Materials::VanishingStress<std::array<Ikarus::Materials::MatrixIndexPair,"
                                     "1ul>{{Ikarus::Materials::MatrixIndexPair{2ul, 2ul}}}," +
                               Dune::className<Material>() + ">"),
                           includes)
                           .first;
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, storedEnergy, 5);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, stresses, 5);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, tangentModuli, 5);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, storedEnergy, 6);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, stresses, 6);
  MAKE_MaterialFunction(shellmaterial, ShellMaterialClass, tangentModuli, 6);

  cls.def("asShellMaterial", [](Material& self) {
    return Materials::shellMaterial(self);
  }); /* no keep_alive since shellMaterial copies the material */
  using BeamMaterialClass = decltype(Materials::beamMaterial(std::declval<Material>()));
  auto beammaterial       = Dune::Python::insertClass<BeamMaterialClass>(
                          scope, std::string("Beam_") + materialname,
                          Dune::Python::GenerateTypeName(
                              "Ikarus::Materials::VanishingStress<std::array<Ikarus::Materials::MatrixIndexPair, "
                                    "2ul>{{Materials::MatrixIndexPair{1, 1},Ikarus::Materials::MatrixIndexPair{2ul, 2ul}}}," +
                              Dune::className<Material>() + ">"),
                          includes)
                          .first;
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, storedEnergy, 4);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, stresses, 4);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, tangentModuli, 4);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, storedEnergy, 6);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, stresses, 6);
  MAKE_MaterialFunction(beammaterial, BeamMaterialClass, tangentModuli, 6);
  cls.def("asBeamMaterial", [](Material& self) {
    return Materials::beamMaterial(self);
  }); /* no keep_alive since beamMaterial copies the material */
}

#define MAKE_MATERIAL_REGISTRY_FUNCTION(name, vecSize)                                      \
  template <class Material, class... options>                                               \
  void register##name(pybind11::handle scope, pybind11::class_<Material, options...> cls) { \
    Ikarus::Python::registerMaterial<Material, vecSize, true>(scope, cls);                  \
  }

MAKE_MATERIAL_REGISTERY_FUNCTION(LinearElasticity, 6);
MAKE_MATERIAL_REGISTERY_FUNCTION(StVenantKirchhoff, 6);
MAKE_MATERIAL_REGISTERY_FUNCTION(NeoHooke, 6);

#if ENABLE_MUESLI

template <class MuesliMaterial, class... options>
void registerMuesliMaterial(pybind11::handle scope, pybind11::class_<MuesliMaterial, options...> cls) {
  Ikarus::Python::registerMaterial<MuesliMaterial, 6, false>(scope, cls);

  cls.def(pybind11::init([](const pybind11::kwargs& kwargs) {
    if constexpr (std::same_as<typename MuesliMaterial::MaterialModel, muesli::neohookeanMaterial> or
                  std::same_as<typename MuesliMaterial::MaterialModel, muesli::svkMaterial> or
                  std::same_as<typename MuesliMaterial::MaterialModel, muesli::elasticIsotropicMaterial>) {
      auto matParameter     = Impl::extractMaterialParameters(kwargs);
      auto muesliParameters = Ikarus::Materials::propertiesFromIkarusMaterialParameters(matParameter);
      return new MuesliMaterial(muesliParameters);
    } else if constexpr (std::same_as<typename MuesliMaterial::MaterialModel, muesli::arrudaboyceMaterial>) {
      bool compressible = kwargs.contains("compressible") ? kwargs["compressible"].cast<bool>() : true;
      return Materials::makeMuesliArrudaBoyce(kwargs["C1"].cast<double>(), kwargs["lambda_m"].cast<double>(),
                                              kwargs["K"].cast<double>(), true);
    } else if constexpr (std::same_as<typename MuesliMaterial::MaterialModel, muesli::yeohMaterial>) {
      bool compressible = kwargs.contains("compressible") ? kwargs["compressible"].cast<bool>() : true;
      return Materials::makeMuesliYeoh(kwargs["C"].cast<std::array<double, 3>>(), kwargs["K"].cast<double>(), true);
    } else if constexpr (std::same_as<typename MuesliMaterial::MaterialModel, muesli::mooneyMaterial>) {
      bool incompressible = kwargs.contains("incompressible") ? kwargs["incompressible"].cast<bool>() : false;
      return Materials::makeMooneyRivlin(kwargs["alpha"].cast<std::array<double, 3>>(), incompressible);
    } else {
      DUNE_THROW(Dune::NotImplemented, "No known constructor for the specified Muesli material mode");
    }
  }));

  cls.def("printDescription", [](MuesliMaterial& self) { self.material().print(std::cout); });
}
#endif

} // namespace Ikarus::Python

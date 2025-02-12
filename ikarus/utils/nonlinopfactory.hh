// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinopfactory.hh
 * \brief Contains the generic NonLinearOperatorFactory class.
 */

#pragma once

#include <utility>

#include <ikarus/assembler/dirichletbcenforcement.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>

namespace Ikarus {

        template<class Signature>
    struct DerivativeTraitsDense
    {
      typedef Dune::Functions::InvalidRange Range;
    };

template<typename K>
struct DerivativeTraitsDense< K(K) >
{
  typedef double Range;
};

template<typename K, Eigen::Index n>
struct DerivativeTraitsDense<K(Eigen::Vector<K,n>)>
{
  typedef Eigen::Vector<K,n> Range;
};

template<typename K,  Eigen::Index n,  Eigen::Index m>
struct DerivativeTraitsDense<Eigen::Vector<K,m>(Eigen::Vector<K,n>)>
{
  typedef Eigen::Matrix<K,m,n> Range;
};

template <typename K,  Eigen::Index n,FESolutions sol, FEParameter para, typename SV, typename PM>
struct DerivativeTraitsDense<Eigen::Vector<K,n>(FERequirements<sol,para,SV,PM>)>
{
  typedef Eigen::Matrix<K,n,n> Range;
};

struct DerivativeTraitsDenseD
{
  template<typename D>
  using Traits= DerivativeTraitsDense<D>;
};


      template<class Signature>
    struct DerivativeTraitsSparse
    {
      typedef Dune::Functions::InvalidRange Range;
    };

template<typename K>
struct DerivativeTraitsSparse< K(K) >
{
  typedef double Range;
};

template<typename K, Eigen::Index n>
struct DerivativeTraitsSparse<K(Eigen::Vector<K,n>)>
{
  typedef Eigen::Vector<K,n> Range;
};

template<typename K,  Eigen::Index n,  Eigen::Index m>
struct DerivativeTraitsSparse<Eigen::Vector<K,m>(Eigen::Vector<K,n>)>
{
  typedef Eigen::SparseMatrix<K> Range;
};

template <typename K,  Eigen::Index n,FESolutions sol, FEParameter para, typename SV, typename PM>
struct DerivativeTraitsSparse<Eigen::Vector<K,n>(FERequirements<sol,para,SV,PM>)>
{
  typedef Eigen::SparseMatrix<K> Range;
};

struct DerivativeTraitsSparseD
{
  template<typename D>
  using Traits= DerivativeTraitsSparse<D>;
};

struct NonLinearOperatorFactory
{
      template<class Signature>
    struct DerivativeTraitsDense
    {
      typedef Dune::Functions::InvalidRange Range;
    };



  template <typename Assembler, typename... Affordances>
  static auto op(Assembler&& as, typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement& req,
                 AffordanceCollection<Affordances...> affordances, DBCOption dbcOption) {
    auto assemblerPtr = [as]() {
      if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                    traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value)
        return as;
      else
        return std::make_shared<std::remove_cvref_t<Assembler>>(std::forward<Assembler>(as));
    }();

    using FERequirement             = typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement;
    [[maybe_unused]] auto KFunction = [dbcOption, assembler = assemblerPtr, affordances](FERequirement& req) -> auto& {
      return assembler->matrix(req, affordances.matrixAffordance(), dbcOption);
    };

    [[maybe_unused]] auto residualFunction = [dbcOption, assembler = assemblerPtr, affordances](FERequirement& req) -> auto& {
      return assembler->vector(req, affordances.vectorAffordance(), dbcOption);
    };



    assert(req.populated() && " Before you calls this method you have to pass populated fe requirements");
    using DerivativeTraitsDummy = std::conditional_t<traits::EigenSparseMatrix<decltype(KFunction(req))>, DerivativeTraitsSparseD, DerivativeTraitsDenseD>;
    if constexpr (affordances.hasScalarAffordance) {
          [[maybe_unused]] auto energyFunction = [assembler = assemblerPtr, affordances](FERequirement& req) -> auto& {

      return assembler->scalar(req, affordances.scalarAffordance());
    };

      using EnergyFunctionSignature = Dune::Functions::SignatureTraits<decltype(energyFunction)>;
      auto sigTag =Dune::Functions::SignatureTag<EnergyFunctionSignature,DerivativeTraitsDummy::template Traits>();
      return Dune::Functions::makeDifferentiableFunctionFromCallables(sigTag,std::move(energyFunction), std::move(residualFunction), std::move(KFunction));
    } else
     {
            using ResidualFunctionSignature = Dune::Functions::SignatureTraits<decltype(residualFunction)>;

            auto sigTag =Dune::Functions::SignatureTag<ResidualFunctionSignature,DerivativeTraitsDummy::template Traits>();

       return Dune::Functions::DifferentiableFunctionFromCallables(sigTag,std::move(residualFunction), std::move(KFunction));}
  }

  template <typename Assembler>
  static auto op(Assembler&& as, DBCOption dbcOption) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to a fe requirement and an affordance collection before you can call "
                 "this method");
    };
    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (as->boundToRequirement() and as->boundToAffordanceCollection()) {
        return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), dbcOption);
      } else {
        ex();
      }
    } else {
      if (as->boundToRequirement() and as->boundToAffordanceCollection()) {
        return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), dbcOption);
      } else {
        ex();
      }
    }
    __builtin_unreachable();
  }

  template <typename Assembler>
  static auto op(Assembler&& as) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to a fe requirement to an affordance collection and to an "
                 "DBCOption before you can call "
                 "this method");
    };
    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (not as->bound())
        ex();
      return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), as->dBCOption());
    } else {
      if (not as.bound())
        ex();
      return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), as.dBCOption());
    }
  }

  template <typename Assembler, typename... Affordances>
  static auto op(Assembler&& as, AffordanceCollection<Affordances...> affordances,
                 DBCOption dbcOption = DBCOption::Full) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to a fe requirement before you can call "
                 "this method");
    };

    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (not as->boundToRequirement())
        ex();
      return op(std::forward<Assembler>(as), as->requirement(), affordances, dbcOption);
    } else {
      if (not as.boundToRequirement())
        ex();
      return op(std::forward<Assembler>(as), as.requirement(), affordances, dbcOption);
    }
  }

  template <typename Assembler>
  static auto op(Assembler&& as, typename traits::remove_pointer_t<std::remove_cvref_t<Assembler>>::FERequirement& req,
                 DBCOption dbcOption) {
    auto ex = []() {
      DUNE_THROW(Dune::InvalidStateException,
                 "Assembler has to be bound to an affordance collection before you can call "
                 "this method");
    };

    if constexpr (std::is_pointer_v<std::remove_cvref_t<Assembler>> or
                  traits::isSharedPtr<std::remove_cvref_t<Assembler>>::value) {
      if (not as->boundToAffordanceCollection())
        ex();
      return op(std::forward<Assembler>(as), as->requirement(), as->affordanceCollection(), dbcOption);
    } else {
      if (not as.boundToAffordanceCollection())
        ex();
      return op(std::forward<Assembler>(as), as.requirement(), as.affordanceCollection(), dbcOption);
    }
  }
};

} // namespace Ikarus

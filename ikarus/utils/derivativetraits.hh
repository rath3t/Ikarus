// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file derivativetraits.hh
 * \brief Contains derivative traits for common vbalue and derivative relations used for makeDifferentiableFunctionFromCallables
 */

#pragma once


#include <dune/functions/common/differentiablefunctionfromcallables.hh>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Ikarus {

  enum class FEParameter;
  enum class FESolutions;
  template <FESolutions sol, FEParameter para, typename SV, typename PM>
  class FERequirements;

        template<class Signature>
    struct DerivativeTraitsDense
    {
      typedef Dune::Functions::InvalidRange Range;
    };

template<typename K> requires std::floating_point<K>
struct DerivativeTraitsDense< K(K) >
{
  typedef double Range;
};

template<typename K, int n>
struct DerivativeTraitsDense<K(Eigen::Vector<K,n>)>
{
  typedef Eigen::Vector<K,n> Range;
};

template <typename K, FESolutions sol, FEParameter para, typename SV, typename PM>
struct DerivativeTraitsDense<K(FERequirements<sol,para,SV,PM>)>
{
  typedef FERequirements<sol,para,SV,PM>::SolutionVectorType Range;
};

template<typename K,  int n,  int m>
struct DerivativeTraitsDense<Eigen::Vector<K,m>(Eigen::Vector<K,n> )>
{
  typedef Eigen::Matrix<K,m,n> Range;
};

template <typename K,  int n,FESolutions sol, FEParameter para, typename SV, typename PM>
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
    struct DerivativeTraitsSparse: DerivativeTraitsDense<Signature>
    {
    };

template<typename K,  int n,  int m>
struct DerivativeTraitsSparse<Eigen::Vector<K,m>(Eigen::Vector<K,n>)>
{
  typedef Eigen::SparseMatrix<K> Range;
};

template <typename K,  int n,FESolutions sol, FEParameter para, typename SV, typename PM>
struct DerivativeTraitsSparse<Eigen::Vector<K,n>(FERequirements<sol,para,SV,PM>)>
{
  typedef Eigen::SparseMatrix<K> Range;
};

struct DerivativeTraitsSparseD
{
  template<typename D>
  using Traits= DerivativeTraitsSparse<D>;
};


} // namespace Ikarus

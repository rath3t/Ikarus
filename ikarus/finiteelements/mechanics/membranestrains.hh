// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file membranestrains.hh
 * \brief Implementation of  membrane strain for shells
 */

#pragma once
#include <ranges>

#include <dune/common/fvector.hh>

#include <Eigen/Core>
namespace Ikarus {

struct DefaultMembraneStrain
{
  /**
   * \brief Compute the strain vector at a given integration point.
   *
   * \param gpPos The position of the integration point.
   * \param geo The geometry object providing the transposed Jacobian.
   * \param uFunction The function representing the displacement field.
   *
   * \tparam GEO The type of the geometry object.
   *
   * \return The strain vector at the given integration point.
   */
  template <typename GEO>
  static auto value(const Dune::FieldVector<double, 2>& gpPos, const GEO& geo,
                    const auto& uFunction) -> Eigen::Vector3<typename std::remove_cvref_t<decltype(uFunction)>::ctype> {
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    Eigen::Vector3<ScalarType> epsV;
    const auto J = Dune::toEigen(geo.jacobianTransposed(gpPos));
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
        uFunction.evaluateDerivative(gpPos, // Here the gpIndex could be passed
                                     Dune::wrt(spatialAll), Dune::on(Dune::DerivativeDirections::referenceElement)));
    const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

    epsV << J.row(0).dot(gradu.col(0)) + 0.5 * gradu.col(0).squaredNorm(),
        J.row(1).dot(gradu.col(1)) + 0.5 * gradu.col(1).squaredNorm(), j.row(0).dot(j.row(1)) - J.row(0).dot(J.row(1));
    return epsV;
  }

  /**
   * \brief Compute the strain-displacement matrix for a given node and integration point.
   *
   * \param gpPos The position of the integration point.
   * \param jcur The Jacobian matrix.
   * \param dNAtGp The derivative of the shape functions at the integration point.
   * \param geo The geometry object of the finite element.
   * \param uFunction The function representing the displacement field.
   * \param localBasis The local basis object.
   * \param node The FE node index.
   *
   * \tparam GEO The type of the geometry object.
   * \tparam ST The scalar type for the matrix elements.
   *
   * \return The strain-displacement matrix for the given node and integration point.
   */
  template <typename GEO, typename ST>
  static auto derivative(const Dune::FieldVector<double, 2>& gpPos, const Eigen::Matrix<ST, 2, 3>& jcur,
                         const auto& dNAtGp, const GEO& geo, const auto& uFunction, const auto& localBasis,
                         const int node) {
    Eigen::Matrix<ST, 3, 3> bop;
    bop.row(0) = jcur.row(0) * dNAtGp(node, 0);
    bop.row(1) = jcur.row(1) * dNAtGp(node, 1);
    bop.row(2) = jcur.row(0) * dNAtGp(node, 1) + jcur.row(1) * dNAtGp(node, 0);

    return bop;
  }

  /**
   * \brief Compute the second derivative of the membrane strain for a given node pair and integration point.
   * \details This function computes the geometric tangent stiffness for a given node pair at a given integration
   * point.
   *
   * \param gpPos The position of the integration point.
   * \param dNAtGp The derivative of the shape functions at the integration point.
   * \param geo The geometry object.
   * \param uFunction The function representing the displacement field.
   * \param localBasis The local basis object.
   * \param S The strain vector.
   * \param I The index of the first node.
   * \param J The index of the second node.
   *
   * \tparam GEO The type of the geometry object.
   * \tparam ST The scalar type for the matrix elements.
   *
   * \return The second derivative of the membrane strain.
   */
  template <typename GEO, typename ST>
  static auto secondDerivative(const Dune::FieldVector<double, 2>& gpPos, const auto& dNAtGp, const GEO& geo,
                               const auto& uFunction, const auto& localBasis, const Eigen::Vector3<ST>& S, int I,
                               int J) {
    const auto& dN1i           = dNAtGp(I, 0);
    const auto& dN1j           = dNAtGp(J, 0);
    const auto& dN2i           = dNAtGp(I, 1);
    const auto& dN2j           = dNAtGp(J, 1);
    const ST NS                = dN1i * dN1j * S[0] + dN2i * dN2j * S[1] + (dN1i * dN2j + dN2i * dN1j) * S[2];
    Eigen::Matrix<ST, 3, 3> kg = Eigen::Matrix<double, 3, 3>::Identity() * NS;
    return kg;
  }
};

} // namespace Ikarus

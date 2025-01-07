// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file scalarwrapper.hh
 * \brief Provides a wrapper for scalar types to support passing by reference in Python bindings.
 *
 * Since Python does not support passing python float by reference to a double&, we have to wrap everything.
 * See also: https://pybind11.readthedocs.io/en/stable/faq.html#limitations-involving-reference-arguments
 */

#pragma once

#include <dune/common/referencehelper.hh>

/**
 * \class ScalarWrapper
 * \brief A wrapper class for scalar types to facilitate reference passing in Python bindings.
 *
 * This class resolves the limitations involving reference arguments in Python, enabling
 * seamless integration with C++ functions expecting reference types.
 *
 * \tparam T The type of the scalar to be wrapped.
 */
template <typename T>
struct ScalarWrapper
{
  /**
   * \brief Constructs a ScalarWrapper with the given value.
   * \param val The initial value of the wrapped scalar.
   */
  ScalarWrapper(T val)
      : val{val} {}

  using RawScalarType = Dune::ResolveRef_t<T>;

  /**
   * \brief Implicitly converts the wrapper to a reference of the raw scalar type.
   * \return A reference to the wrapped scalar value.
   */
  operator RawScalarType&() const { return Dune::resolveRef(val); }

  /**
   * \brief Gets the wrapped scalar value as a constant reference.
   * \return A constant reference to the wrapped scalar value.
   */
  const RawScalarType& value() const { return Dune::resolveRef(val); }

  /**
   * \brief Gets the wrapped scalar value as a reference.
   * \return A reference to the wrapped scalar value.
   */
  RawScalarType& value() { return Dune::resolveRef(val); }

  /**
   * \brief Adds the values of two ScalarWrapper instances. This returns the raw type since this makes makes no since
   * when T is of reference_wrapper type \param v Another ScalarWrapper instance. \return The result of the addition as
   * a new scalar value.
   */
  RawScalarType operator+(const ScalarWrapper& v) const {
    return RawScalarType{Dune::resolveRef(val) + Dune::resolveRef(v.val)};
  }
  /**
   * \brief Subtracts the value of another ScalarWrapper from this instance.
   * \param v Another ScalarWrapper instance.
   * \return The result of the subtraction as a new scalar value.
   */
  RawScalarType operator-(const ScalarWrapper& v) const {
    return RawScalarType{Dune::resolveRef(val) - Dune::resolveRef(v.val)};
  }

  /**
   * \brief Negates the wrapped scalar value.
   * \return The negated scalar value.
   */
  RawScalarType operator-() const { return RawScalarType{-Dune::resolveRef(val)}; }

  /**
   * \brief Multiplies the wrapped scalar value by another value.
   * \param value The value to multiply with.
   * \return The result of the multiplication.
   */
  RawScalarType operator*(RawScalarType value) const { return RawScalarType{val * value}; }

  /**
   * \brief Adds another ScalarWrapper's value to this instance.
   * \param v Another ScalarWrapper instance.
   * \return A reference to this instance after the addition.
   */
  ScalarWrapper& operator+=(const ScalarWrapper& v) {
    Dune::resolveRef(val) += Dune::resolveRef(v.val);
    return *this;
  }

  /**
   * \brief Subtracts another ScalarWrapper's value from this instance.
   * \param v Another ScalarWrapper instance.
   * \return A reference to this instance after the subtraction.
   */
  ScalarWrapper& operator-=(const ScalarWrapper& v) {
    Dune::resolveRef(val) -= Dune::resolveRef(v.val);
    return *this;
  }

  /**
   * \brief Multiplies the wrapped scalar value by another value and assigns the result.
   * \param v The value to multiply with.
   * \return A reference to this instance after the multiplication.
   */
  ScalarWrapper& operator*=(RawScalarType v) {
    Dune::resolveRef(val) *= v;
    return *this;
  }

  /**
   * \brief Divides the wrapped scalar value by another value and assigns the result.
   * \param v The value to divide by.
   * \return A reference to this instance after the division.
   */
  ScalarWrapper& operator/=(RawScalarType v) {
    Dune::resolveRef(val) /= v;
    return *this;
  }

private:
  T val;

  /**
   * \brief Multiplies a scalar value by a ScalarWrapper's value.
   * \param f A scalar value.
   * \param v A ScalarWrapper instance.
   * \return The result of the multiplication.
   */
  friend RawScalarType operator*(RawScalarType f, const ScalarWrapper& v) {
    return RawScalarType{f * Dune::resolveRef(v.val)};
  }
};
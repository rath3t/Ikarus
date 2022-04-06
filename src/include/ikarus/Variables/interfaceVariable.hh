//
// Created by Alex on 25.05.2021.
//

#pragma once

#include <memory>
#include <span>

#include <Eigen/Core>

#include <ikarus/variables/variablePolicies.hh>
#include <ikarus/utils/linearAlgebraTypedefs.hh>

namespace Ikarus {
  enum class VariableTags;

  class IVariable {
  public:
    template <Concepts::Variable VAR>
    explicit IVariable(const VAR &vo) : variableImpl{std::make_unique<VarImpl<VAR> >(vo)} {}

    ~IVariable() = default;
    IVariable(const IVariable &other) : variableImpl{other.variableImpl->clone()} {}
    IVariable &operator=(IVariable &&) noexcept = default;
    IVariable(IVariable &&) noexcept            = default;

    IVariable &operator=(const IVariable &other) {
      IVariable tmp(other);  // Temporary-swap idiom
      std::swap(variableImpl, tmp.variableImpl);
      return *this;
    }

    using UpdateType     = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 8, 1>;
    using CoordinateType = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 8, 1>;

    double &operator[](int i);

  private:
    struct VarBase {
      virtual ~VarBase() = default;

      [[nodiscard]] virtual int do_valueSize() const                              = 0;
      [[nodiscard]] virtual int do_correctionSize() const                         = 0;
      [[nodiscard]] virtual bool do_equalComparison(const IVariable &other) const = 0;
      [[nodiscard]] virtual bool do_lessComparison(const IVariable &other) const  = 0;
      virtual void do_assignAdd(const UpdateType &other)                          = 0;
      virtual void do_setValue(const UpdateType &other)                           = 0;
      [[nodiscard]] virtual CoordinateType do_getValue() const                    = 0;
      [[nodiscard]] virtual int do_getTag() const                                 = 0;
      [[nodiscard]] virtual const double &operator[](int i) const                 = 0;
      [[nodiscard]] virtual double &operator[](int i)                             = 0;
      [[nodiscard]] virtual std::unique_ptr<VarBase> clone() const                = 0;
    };

    template <typename VAR>
    struct VarImpl : public VarBase {
      explicit VarImpl(VAR voarg) : vo{voarg} {};

      [[nodiscard]] int do_valueSize() const final { return VAR::valueSize; }
      [[nodiscard]] int do_correctionSize() const final { return VAR::correctionSize; }
      [[nodiscard]] int do_getTag() const final { return vo.getTag(); };
      [[nodiscard]] CoordinateType do_getValue() const final { return vo.getValue(); };
      void do_assignAdd(const UpdateType &other) final { vo.update(other); }
      void do_setValue(const UpdateType &other) final { return vo.setValue(other); }
      [[nodiscard]] bool do_equalComparison(const IVariable &other) const final {
        return (this->do_getTag() == getTag(other));
      };

      [[nodiscard]] bool do_lessComparison(const IVariable &other) const final {
        return (this->do_getTag() < getTag(other));
      };

      [[nodiscard]] const double &operator[](int i) const final { return vo[i]; };
      [[nodiscard]] double &operator[](int i) final { return vo[i]; };

      [[nodiscard]] std::unique_ptr<VarBase> clone() const final { return std::make_unique<VarImpl>(*this); }

      VAR vo;
    };

    std::unique_ptr<VarBase> variableImpl;  // Pimpl idiom / Bridge Design Pattern

    friend IVariable &operator+=(IVariable &vo, const UpdateType &correction);
    friend IVariable operator+(IVariable &vo, const UpdateType &correction);
    friend IVariable &operator-=(IVariable &vo, const UpdateType &correction);
    friend IVariable operator-(IVariable &vo, const UpdateType &correction);

    friend void setValue(IVariable &vo, const UpdateType &value);
    friend CoordinateType getValue(const IVariable &vo);
    friend int valueSize(const IVariable &vo);
    friend int correctionSize(const IVariable &vo);
    friend bool operator==(const IVariable &var, const IVariable &other);
    friend bool operator<(const IVariable &var, const IVariable &other);
    friend int getTag(const IVariable &var);
    friend std::ostream &operator<<(std::ostream &s, const IVariable &var);
  };

  IVariable &operator+=(IVariable &vo, const IVariable::UpdateType &correction);
  IVariable &operator+=(IVariable *vo, const IVariable::UpdateType &correction);
  IVariable operator+(IVariable &vo, const IVariable::UpdateType &correction);
  IVariable operator+(IVariable *vo, const IVariable::UpdateType &correction);
  IVariable &operator-=(IVariable &vo, const IVariable::UpdateType &correction);
  IVariable &operator-=(IVariable *vo, const IVariable::UpdateType &correction);
  IVariable operator-(IVariable &vo, const IVariable::UpdateType &correction);
  IVariable operator-(IVariable *vo, const IVariable::UpdateType &correction);
  void setValue(IVariable &vo, const IVariable::UpdateType &value);
  IVariable::CoordinateType getValue(const IVariable &vo);
  IVariable::CoordinateType getValue(const IVariable *vo);
  int valueSize(const IVariable &vo);
  int correctionSize(const IVariable &vo);
  bool operator==(const IVariable &var, const IVariable &other);
  bool operator<(const IVariable &var, const IVariable &other);
  int getTag(const IVariable &var);
  std::string getName(const IVariable &var);
  std::ostream &operator<<(std::ostream &s, const IVariable &var);
  std::ostream &operator<<(std::ostream &s, const IVariable *var);
  size_t valueSize(std::span<const IVariable> varSpan);
  size_t correctionSize(std::span<const IVariable> varSpan);
  void update(std::span<IVariable> varSpan, const Eigen::VectorXd &correction);
  bool isType(const IVariable &vo, Ikarus::VariableTags tag);
  bool isType(IVariable *vo, Ikarus::VariableTags tag);
}  // namespace Ikarus
// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"
#include "tests/src/testcommon.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/febase.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverstate.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename PreFE, typename FE>
class DummySkill;

struct DummySkillPre
{
  template <typename PreFE, typename FE>
  using Skill = DummySkill<PreFE, FE>;
};

MAKE_ENUM(UpdateMessages, INCREMENT, RESET)

template <typename PreFE, typename FE>
class DummySkill : public ResultTypeBase<ResultTypes::linearStress>
{
public:
  using Traits       = PreFE::Traits;
  using BasisHandler = typename Traits::BasisHandler;
  using FlatBasis    = typename Traits::FlatBasis;
  using Requirement =
      FERequirementsFactory<FESolutions::displacement, FEParameter::loadfactor, Traits::useEigenRef>::type;
  using LocalView = typename Traits::LocalView;
  using Geometry  = typename Traits::Geometry;
  using GridView  = typename Traits::GridView;
  using Element   = typename Traits::Element;
  using Pre       = DummySkillPre;

  explicit DummySkill(const Pre& pre) {}

protected:
  // This returns a tuple functions to be registered
  template <typename MT, typename BC>
  void subscribeToImpl(BC& bc) {
    if constexpr (std::same_as<MT, NonLinearSolverMessages>) {
      using NLSState = typename BC::State;
      underlying().subscribe(bc, [&](NonLinearSolverMessages message, NLSState& state) {
        this->updateState(message, state.solution, state.correction);
      });
    } else if constexpr (std::same_as<MT, UpdateMessages>) {
      underlying().subscribe(bc, [&](UpdateMessages message, int val) { this->updateState(message, val); });
      underlying().subscribe(bc, [&](UpdateMessages message) { this->updateState(message); });
    } else if constexpr (std::same_as<MT, ControlMessages>) {
      underlying().subscribe(bc, [&](ControlMessages message) { this->updateState(message); });
    } else
      static_assert(Dune::AlwaysFalse<MT>::value, "No registration for MT");
  }

  template <class ScalarType>
  auto calculateScalarImpl(const Requirement& par, ScalarAffordance affo,
                           const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx =
                               std::nullopt) const -> ScalarType {
    return counter_;
  }

  template <typename ScalarType>
  void calculateMatrixImpl(
      const Requirement& par, const MatrixAffordance& affordance, typename Traits::template MatrixType<> K,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    K.setConstant(counter_);
  }

  template <typename ScalarType>
  void calculateVectorImpl(
      const Requirement& par, VectorAffordance affordance, typename Traits::template VectorType<ScalarType> force,
      const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx = std::nullopt) const {
    force.setConstant(counter_);
  }

  template <template <typename, int, int> class RT>
  requires(supportsResultType<RT>())
  auto calculateAtImpl(const Requirement& req, [[maybe_unused]] const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<0>) const {}

  void updateState(NonLinearSolverMessages message, const Eigen::VectorXd& vec,
                   const std::remove_reference_t<typename Traits::template VectorType<>>& correction) {
    // We are hijacking the NLSolverMessages here to get either increment the counter or reset it to zero
    if (message == NonLinearSolverMessages::FINISHED_SUCESSFULLY)
      counter_ = 0;
    else
      counter_ = counter_ + 1;
  }

  void updateState(UpdateMessages message, double val) {
    // We are using a second set of UpdateMessages
    counter_ = val;
  }
  void updateState(UpdateMessages message) {
    if (message == UpdateMessages::RESET)
      counter_ = 0;
  }
  void updateState(ControlMessages message) {
    if (message == ControlMessages::CONTROL_STARTED)
      counter_ = 10;
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }

  int counter_{0};
};

inline auto dummySkill() {
  DummySkillPre pre{};

  return pre;
}
using NRStateDummy = NonlinearSolverState<const Eigen::VectorXd&, const Eigen::VectorXd&>;
struct DummyBroadcaster : public Broadcasters<void(NonLinearSolverMessages, NRStateDummy& state),
                                              void(UpdateMessages, int), void(UpdateMessages)>
{
  using State = NRStateDummy;

  void emitMessage(NonLinearSolverMessages message) {
    auto state = State(x_, dx_);
    this->notify(message, state);
  }
  void emitMessage(UpdateMessages message, int val) { this->notify(message, val); }
  void emitMessage(UpdateMessages message) { this->notify(message); }

  DummyBroadcaster(size_t size)
      : size_(size) {
    x_  = Eigen::VectorXd::Zero(size);
    dx_ = Eigen::VectorXd::Zero(size);
  }

private:
  size_t size_;
  Eigen::VectorXd x_;
  Eigen::VectorXd dx_;
};

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  constexpr double h = 1.0;
  constexpr double L = 10.0;
  auto grid          = createGrid<Grids::OneDFoamGridIn2D>();
  auto gridView      = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  auto sk      = skills(dummySkill());
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  DirichletValues dirichletValues(basis.flat());
  auto sparseFlatAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  /// Create non-linear operator
  double lambda = 0.0;
  Eigen::VectorXd d;
  size_t size = basis.flat().size();
  d.setZero(size);
  auto req = FEType::Requirement();
  req.insertGlobalSolution(d).insertParameter(lambda);

  sparseFlatAssembler->bind(req);
  sparseFlatAssembler->bind(Ikarus::AffordanceCollections::elastoStatics);
  sparseFlatAssembler->bind(Ikarus::DBCOption::Full);

  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);
  Eigen::VectorXd F = Eigen::VectorXd::Zero(8);

  auto& fe = fes.front();

  auto checkMatrixAndVector = [&](int counter, const std::string& messageIfFailed) {
    fe.calculateMatrixImpl<double>(req, MatrixAffordance::stiffness, K);
    fe.calculateVectorImpl<double>(req, VectorAffordance::forces, F);

    t.check(K.sum() == counter * 8 * 8) << messageIfFailed;
    t.check(F.sum() == counter * 8) << messageIfFailed;
  };

  // Counter is zero, we haven't updated yet
  checkMatrixAndVector(0, testLocation());

  DummyBroadcaster broadcaster(size);
  for (auto& fe : fes) {
    fe.subscribeTo<NonLinearSolverMessages>(broadcaster);
    fe.subscribeTo<UpdateMessages>(broadcaster);
  }

  broadcaster.emitMessage(NonLinearSolverMessages::SOLUTION_CHANGED);
  checkMatrixAndVector(1, testLocation());

  // wenn we emit NonLinearSolverMessages::FINISHED_SUCESSFULLY, counter should reset
  broadcaster.emitMessage(NonLinearSolverMessages::FINISHED_SUCESSFULLY);
  checkMatrixAndVector(0, testLocation());

  broadcaster.emitMessage(UpdateMessages::INCREMENT, 2);
  checkMatrixAndVector(2, testLocation());

  broadcaster.emitMessage(UpdateMessages::RESET);
  checkMatrixAndVector(0, testLocation());

  // Now check with lc
  auto linSolver = LinearSolver(SolverTypeTag::d_LDLT);
  NewtonRaphsonConfig nrConfig({}, linSolver);
  NonlinearSolverFactory nrFactory(nrConfig);
  auto nr       = nrFactory.create(sparseFlatAssembler);
  auto nonLinOp = Ikarus::NonLinearOperatorFactory::op(sparseFlatAssembler);
  auto lc       = ControlRoutineFactory::create(LoadControlConfig{1, 0.0, 1.0}, nr, sparseFlatAssembler);

  lc.notify(Ikarus::ControlMessages::CONTROL_STARTED);
  checkMatrixAndVector(10, testLocation());

  // Set to 2, then unsubscribe a unrelated listener, set to 6, unsubscribe from all, then set to 4, but it should still
  // be 2
  broadcaster.emitMessage(UpdateMessages::INCREMENT, 2);
  fe.unSubscribeLast();
  // This should not deregister the following message listener method
  broadcaster.emitMessage(UpdateMessages::INCREMENT, 6);
  checkMatrixAndVector(6, testLocation());
  // This also deregisters this listener method
  fe.unSubscribeAll();
  broadcaster.emitMessage(UpdateMessages::INCREMENT, 4);
  checkMatrixAndVector(6, testLocation());

  // Using lower-level api
  int i       = 0;
  auto callee = [&]() { i++; };
  auto token  = fe.subscribe(broadcaster, [&](UpdateMessages message) { callee(); });

  broadcaster.emitMessage(UpdateMessages::INCREMENT);
  t.check(i == 1) << testLocation();

  fe.unSubscribe(std::move(token));
  broadcaster.emitMessage(UpdateMessages::INCREMENT);
  t.check(i == 1) << testLocation();

  static_assert(Concepts::PointerOrSmartPointer<std::shared_ptr<int>>);

  return t.exit();
}

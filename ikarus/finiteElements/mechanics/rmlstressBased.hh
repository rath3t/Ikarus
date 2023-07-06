// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "externalload.hh"
#include "fesettings.hh"

#include <utility>

//#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/tensorproductquadrule.hh>
namespace Ikarus {


  template<typename ResultantBasedShell>
  struct StressBasedShellRM : private ResultantBasedShell
  {
    template <typename ScalarType>
    using KinematicVariables =typename ResultantBasedShell::KinematicVariables<ScalarType>;
    Dune::QuadratureRule<double,3> rule;
    static constexpr int midSurfaceDim         = 3;
    static constexpr int directorDim           = 3;
    static constexpr int directorCorrectionDim = directorDim - 1;
    template <typename... Args>
    explicit StressBasedShellRM(Args&&... args) : ResultantBasedShell{std::forward<Args>(args)...},geo{*this->geo_} {

      const auto& twoDRule = Dune::QuadratureRules<double,2>::rule(Dune::GeometryTypes::quadrilateral,this->order);
      const auto& oneDRule = Dune::QuadratureRules<double,1>::rule(Dune::GeometryTypes::line,2);

      numberOfThicknessIntegrationPoints = oneDRule.size();
      rule = Ikarus::tensorProductQuadrature(twoDRule,oneDRule);
//      fESettings = ResultantBasedShell::getfFESettings();
//      order = ResultantBasedShell::getOrder();
//      numberOfNodes = ResultantBasedShell::getNumberOfNodes();
//      membraneStrain = ResultantBasedShell::getMembraneStrain();
//      std::cout<<"Constructor: "<<rule.size()<<std::endl;

    }

    int numberOfThicknessIntegrationPoints;
//    FESettings fESettings;
//    size_t numberOfNodes{0};
//    int order{};
//    mutable MembraneStrain membraneStrain;

    using GlobalIndex = typename ResultantBasedShell::GlobalIndex;
    using ResultantBasedShell::size ;
    using ResultantBasedShell::globalFlatIndices ;
    using ResultantBasedShell::localView ;
    using ResultantBasedShell::ndofEmbedded ;
    using Basis = typename ResultantBasedShell::Basis;
    using FlatBasis = typename Basis::FlatBasis;
    using FERequirementType = typename ResultantBasedShell::FERequirementType;
    using ResultRequirementsType = ResultRequirements<FERequirementType>;
    using LocalView = typename FlatBasis::LocalView;
    using Element = typename LocalView::Element;
    using Geometry = typename Element::Geometry;
    const Geometry& geo;
    using GridView = typename FlatBasis::GridView;
    static constexpr int useEigenRef = ResultantBasedShell::useEigenRef;
    using Traits = TraitsFromLocalView<LocalView, useEigenRef>;
   protected:
    static constexpr int myDim = Traits::mydim;
    static constexpr int worlddim = Traits::worlddim;
    static constexpr bool isEuclidean = ResultantBasedShell::isEuclidean;
    using ResultantBasedShell::localViewBlocked ;
   public:
    inline void calculateMatrix(const FERequirementType &par, typename Traits::template MatrixType<> K) const {
      calculateMatrixImpl<double>(par, K);
    }

    inline void calculateVector(const FERequirementType &par, typename Traits::template VectorType<> force) const {
      calculateVectorImpl<double>(par, force);
    }
    inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }
    template<typename ScalarType>
    auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = this->createFunctions(par,dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      ScalarType energy = 0.0;
      this->membraneStrain.pre(geo,displacementFunction);

      const auto &thickness_ = this->fESettings.template request<double>("thickness");
      KinematicVariables<ScalarType> kin{};
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};

        kin.t           = directorFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));

        auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,gpIndex2D,geo,displacementFunction,kin);

//        const auto
//            [C, epsV, kappaV, j, J, h,H, a3N, a3] = this->computeMaterialAndStrains(gp2DPos, gpIndex2D, geo, uFunction);

        const auto G = this->calc3DMetric(kin,zeta);
//        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto C3D = this->materialTangent(Ginv);
        Eigen::Vector3<ScalarType> rhoV;

        rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
        const auto strainsV= (epsV+ zeta*kappaV+zeta*zeta*rhoV).eval();
        Eigen::Vector<ScalarType,5> strains;
        strains<< strainsV,gammaV;
        const ScalarType energyVal = 0.5*strains.dot(C3D*strains);
        energy += (energyVal)*geo.integrationElement(gp2DPos)*gp.weight()*2*thickness_/2.0;
        ++gpIndex;
      }

      if (this->volumeLoad) {
        for (const auto& [gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
          const auto u                                       = displacementFunction.evaluate(gpIndex);
          const Eigen::Vector<double, Traits::worlddim> fExt = this->volumeLoad(geo.global(gp.position()), lambda)[0];
          energy -= u.dot(fExt) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e., neumann boundary
      if (not this->neumannBoundary and not this->neumannBoundaryLoad) return energy;
      forEachInterSectionIntegrationPoint(this->localView().element(),this->neumannBoundary,this->order,
                                          [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                            const auto neumannValue
                                                = this->neumannBoundaryLoad(globalPos, lambda)[0];
                                            const auto u = displacementFunction.evaluate(quadPos);
                                            energy -= neumannValue.dot(u) * intElement;
                                          });
      return energy;
    }

    template<typename ScalarType>
    void calculateVectorImpl(const FERequirementType &par, typename Traits::template VectorType<ScalarType> force,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      force.setZero();
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = this->createFunctions(par,dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      const int midSurfaceDofs = this->numNodeMidSurface * midSurfaceDim;
      const auto &thickness_ = this->fESettings.template request<double>("thickness");
      this->membraneStrain.pre(geo,displacementFunction);
      KinematicVariables<ScalarType> kin{};
      // Internal forces
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceI;
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorI;
      std::cout<<"calculateVectorImpl"<<std::endl;
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};
        kin.t           = directorFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        const Eigen::Matrix<ScalarType,2,3> jE = kin.a1anda2.transpose();
        auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,gpIndex2D,geo,displacementFunction,kin);

        const auto G = this->calc3DMetric(kin,zeta);
        //        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto g1Andg2 = (kin.a1anda2+ zeta* kin.td1Andtd2).eval();

        const auto C3D = this->materialTangent(Ginv);
        Eigen::Vector3<ScalarType> rhoV;
        rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
        const auto strainsV= (epsV+ zeta*kappaV+zeta*zeta*rhoV).eval();
        Eigen::Vector<ScalarType,5> strains;
        strains<< strainsV,gammaV;
        const auto S = (C3D*strains).eval();

        const auto &Nd = this->localBasisMidSurface.evaluateJacobian(gpIndex2D);
        for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
          std::cout<<"bopMidSurfaceIForces"<<std::endl;
          bopMidSurfaceI= bopMembrane(kin, gpIndex2D, gp2DPos,jE, Nd, displacementFunction, i, zeta);

          force.template segment<3>(3*i) +=
              bopMidSurfaceI.transpose()*S*geo.integrationElement(gp2DPos)*gp.weight()*2*thickness_/2.0;
          // the first two fixes the change of the integration mapping from 0..1 to -1..1,
          // and the h/2 factor is the factor for the correct thickness
        }
        for (int i = 0; i < this->numNodeDirector; ++i) {
          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
          std::cout<<"bopDirectorIForces"<<std::endl;

          bopDirectorI =bopDirector(kin, g1Andg2, gpIndex2D, directorFunction, i, zeta);
          force.template segment<directorCorrectionDim>(indexI) += bopDirectorI.transpose() * S*geo.integrationElement(gp2DPos)*gp.weight()*2*thickness_/2.0;
        }

        ++gpIndex;
      }
      using namespace Dune::Indices;
      //External forces volume forces over the domain
      if (this->volumeLoad) {
        for (const auto& [gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
          const auto [fext,mext] = this->volumeLoad(geo.global(gp.position()), lambda);
          for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
            const auto udCi = displacementFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<worlddim>(worlddim * i)
                -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
          }
          for (size_t i = 0; i < this->numNodeDirector; ++i) {
            const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
            const auto tdCi = directorFunction.evaluateDerivative(gpIndex, Dune::wrt(coeff(_1,i)));
            force.template segment<directorCorrectionDim>(indexI)
                -= (tdCi.transpose() * mext) * geo.integrationElement(gp.position()) * gp.weight();
          }
        }
      }

      // External forces, boundary forces, i.e., at the Neumann boundary
      if (not this->neumannBoundary) return;
      forEachInterSectionIntegrationPoint(this->localView().element(),this->neumannBoundary,this->order,
                                          [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                            const auto [fext,mext] = this->neumannBoundaryLoad(globalPos, lambda);

                                            for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
                                              const auto udCi = displacementFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

                                              // Value of the Neumann data at the current position
                                              force.template segment<worlddim>(worlddim * i) -= udCi * fext * intElement;
                                            }
                                            for (size_t i = 0; i < this->numNodeDirector; ++i) {
                                              const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
                                              const auto tdCi = directorFunction.evaluateDerivative(quadPos, Dune::wrt(coeff(_1,i)));
                                              force.template segment<directorCorrectionDim>(indexI)
                                                  -= (tdCi.transpose() * mext) * intElement;
                                            }
                                          });
    }
    template<typename ScalarType>
    auto bopMembrane(const KinematicVariables<ScalarType>& kin, const auto& gpIndex2D,const auto& gp2DPos,const auto& jE, const auto& Nd, const auto& displacementFunction,int i,double zeta) const
    {
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceI;
      const auto bopIMembrane   = this->membraneStrain.derivative(gp2DPos,jE, Nd, geo,displacementFunction,this->localBasisMidSurface, i);
      const auto bopIBending   = this->boperatorMidSurfaceBending(kin,gpIndex2D, i,displacementFunction);
      const auto bopIShear   = this->boperatorMidSurfaceShear(kin,gpIndex2D, i,displacementFunction);
      bopMidSurfaceI<< bopIMembrane+zeta*bopIBending, bopIShear;
      return bopMidSurfaceI;
    }
    template<typename ScalarType>
    auto bopDirector(const KinematicVariables<ScalarType>& kin, const auto& g1Andg2, const auto& gpIndex2D,const auto& directorFunction,int i,double zeta)const
    {
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorI;

      const auto bopBendingI   = this->boperatorDirectorBending(g1Andg2, gpIndex2D, i, directorFunction);
      const auto bopShearI   = this->boperatorDirectorShear(kin, gpIndex2D, i, directorFunction);
      bopDirectorI<<zeta*bopBendingI,bopShearI;
      return bopDirectorI;

    }
    template<typename ScalarType>
    void calculateMatrixImpl(const FERequirementType &par, typename Traits::template MatrixType<ScalarType> K,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      K.setZero();
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = this->createFunctions(par,dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      this->membraneStrain.pre(geo,displacementFunction);
      const int midSurfaceDofs = this->numNodeMidSurface * midSurfaceDim;
      KinematicVariables<ScalarType> kin{};


      const auto &thickness_ = this->fESettings.template request<double>("thickness");
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceI;
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorI;
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceJ;
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorJ;
      // Internal forces
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};
        kin.t           = directorFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        const Eigen::Matrix<ScalarType,2,3> jE = kin.a1anda2.transpose();
        auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,gpIndex2D,geo,displacementFunction,kin);

        const auto G = this->calc3DMetric(kin,zeta);
        //        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto g1Andg2 = (kin.a1anda2+ zeta* kin.td1Andtd2).eval();

        const auto C3D = this->materialTangent(Ginv);
        Eigen::Vector3<ScalarType> rhoV;
        rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
        const auto strainsV= (epsV+ zeta*kappaV+zeta*zeta*rhoV).eval();
        Eigen::Vector<ScalarType,5> strains;
        strains<< strainsV,gammaV;
        const auto S = (C3D*strains).eval();
        Eigen::Vector<ScalarType,8> S8;
        Eigen::Vector<ScalarType,3> SSec;
        SSec= zeta*zeta*S.template segment<3>(0);
        S8.template segment<3>(0)<< S.template segment<3>(0);
        S8.template segment<3>(3)<< zeta*S.template segment<3>(0);
        S8.template segment<2>(6)<< S.template segment<2>(3);
        const auto intElement =geo.integrationElement(gp2DPos) * gp.weight() * 2 * thickness_ / 2.0;

        const auto &Nd = this->localBasisMidSurface.evaluateJacobian(gpIndex2D);
        const auto &N = this->localBasisMidSurface.evaluateFunction(gpIndex2D);
        const auto &dNdirector = this->localBasisDirector.evaluateJacobian(gpIndex2D);
        const auto &Ndirector = this->localBasisDirector.evaluateFunction(gpIndex2D);
        for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
          const auto indexI = midSurfaceDim * i;
          bopMidSurfaceI= bopMembrane(kin, gpIndex2D, gp2DPos,jE, Nd, displacementFunction, i, zeta);

          for (size_t j = i; j < this->numNodeMidSurface; ++j) {
            const auto indexJ = midSurfaceDim * j;
            std::cout<<"bopMidSurfaceIStiffness"<<std::endl;

            bopMidSurfaceI= bopMembrane(kin, gpIndex2D, gp2DPos,jE, Nd, displacementFunction, j, zeta);

            K.template block<3, 3>(indexI, indexJ) += (bopMidSurfaceI.transpose()*C3D*bopMidSurfaceJ)*intElement;

            Eigen::Matrix<ScalarType, 3, 3> kgMembraneIJ = this->membraneStrain.secondDerivative(gp2DPos,Nd,geo,displacementFunction,this->localBasisMidSurface, S.template segment<3>(0).eval(), i, j);
            K.template block<3, 3>(indexI, indexJ) += kgMembraneIJ*intElement;

          }
          for (int j = 0; j < this->numNodeDirector; ++j) {
            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;
            const auto bopBendingJ   = this->boperatorDirectorBending(g1Andg2, gpIndex2D, j, directorFunction);
            const auto bopShearJ   = this->boperatorDirectorShear(kin, gpIndex2D, j, directorFunction);
            std::cout<<"bopDirectorIStiffness"<<std::endl;
            bopDirectorJ =bopDirector(kin, g1Andg2, gpIndex2D, directorFunction, j, zeta);

            Eigen::Matrix<ScalarType, 3, 2> kg= this->kgMidSurfaceDirectorBending(kin,N,Nd,gpIndex2D,i,j,displacementFunction,directorFunction,S8);
            Eigen::Matrix<ScalarType, 3, 2> kg2= this->kgMidSurfaceDirectorShear(kin,N,Nd,gpIndex2D,i,j,displacementFunction,directorFunction,S8);

            K.template block<midSurfaceDim, directorCorrectionDim>(indexI, indexJ)
                += bopMidSurfaceI.transpose() * C3D * bopDirectorJ * intElement;

            K.template block<midSurfaceDim, directorCorrectionDim>(indexI, indexJ)
                += (kg+kg2) * intElement;
          }
        }
        for (int i = 0; i < this->numNodeDirector; ++i) {
          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
          bopDirectorI =bopDirector(kin, g1Andg2, gpIndex2D, directorFunction, i, zeta);

          for (int j = i; j < this->numNodeDirector; ++j) {
            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;

            bopDirectorJ =bopDirector(kin, g1Andg2, gpIndex2D, directorFunction, j, zeta);
            std::cout<<"bopDirectorJ"<<std::endl;
            std::cout<<bopDirectorJ<<std::endl;
            Eigen::Matrix<ScalarType, 2, 2> kgBending= this->kgDirectorDirectorBending(kin,Ndirector,dNdirector,gpIndex2D,i,j,displacementFunction,directorFunction,S8);
            std::cout<<"Before"<<std::endl;
            std::cout<<"kgBending"<<std::endl;
            std::cout<<kgBending<<std::endl;
//            Eigen::Matrix<ScalarType, 2, 2> kgBending2= this->kgSecondDirectorDirectorBending(kin,Ndirector,dNdirector,gpIndex2D,i,j,displacementFunction,directorFunction,SSec);
//            std::cout<<"kgBending2"<<std::endl;
//            std::cout<<kgBending2<<std::endl;
            Eigen::Matrix<ScalarType, 2, 2> kgShear= this->kgDirectorDirectorShear(kin,Ndirector,dNdirector,gpIndex2D,i,j,displacementFunction,directorFunction,S8);
            std::cout<<"kgShear"<<std::endl;
            std::cout<<kgShear<<std::endl;

            K.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ)
                += bopDirectorI.transpose() * C3D * bopDirectorJ * intElement;

            K.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ)
                += (kgBending+kgShear)* intElement;
          }
        }
        ++gpIndex;
      }
      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
          }
  };

}  // namespace Ikarus

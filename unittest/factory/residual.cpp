///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2021, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "residual.hpp"
#include "pinocchio/multibody/geometry.hpp"
#include "crocoddyl/multibody/residuals/state.hpp"
#include "crocoddyl/core/residuals/control.hpp"
#include "crocoddyl/multibody/residuals/com-position.hpp"
#include "crocoddyl/multibody/residuals/centroidal-momentum.hpp"
#include "crocoddyl/multibody/residuals/frame-placement.hpp"
#include "crocoddyl/multibody/residuals/frame-rotation.hpp"
#include "crocoddyl/multibody/residuals/frame-translation.hpp"
#include "crocoddyl/multibody/residuals/frame-velocity.hpp"
#include "crocoddyl/multibody/residuals/control-gravity.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance-sqr.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance-exp.hpp"
#include "crocoddyl/multibody/residuals/fly-high-sqr.hpp"
#include "crocoddyl/multibody/residuals/fly-high-exp.hpp"
#include "crocoddyl/core/utils/exception.hpp"

namespace crocoddyl {
namespace unittest {

const std::vector<ResidualModelTypes::Type> ResidualModelTypes::all(ResidualModelTypes::init_all());

std::ostream& operator<<(std::ostream& os, ResidualModelTypes::Type type) {
  switch (type) {
    case ResidualModelTypes::ResidualModelState:
      os << "ResidualModelState";
      break;
    case ResidualModelTypes::ResidualModelControl:
      os << "ResidualModelControl";
      break;
    case ResidualModelTypes::ResidualModelCoMPosition:
      os << "ResidualModelCoMPosition";
      break;
    case ResidualModelTypes::ResidualModelCentroidalMomentum:
      os << "ResidualModelCentroidalMomentum";
      break;
    case ResidualModelTypes::ResidualModelFramePlacement:
      os << "ResidualModelFramePlacement";
      break;
    case ResidualModelTypes::ResidualModelFrameRotation:
      os << "ResidualModelFrameRotation";
      break;
    case ResidualModelTypes::ResidualModelFrameTranslation:
      os << "ResidualModelFrameTranslation";
      break;
    case ResidualModelTypes::ResidualModelFrameVelocity:
      os << "ResidualModelFrameVelocity";
      break;
    case ResidualModelTypes::ResidualModelControlGrav:
      os << "ResidualModelControlGrav";
      break;
    case ResidualModelTypes::ResidualModelObstacleAvoidance:
      os << "ResidualModelObstacleAvoidance";
      break;
    case ResidualModelTypes::ResidualModelObstacleAvoidanceSqr:
      os << "ResidualModelObstacleAvoidanceSqr";
      break;
    case ResidualModelTypes::ResidualModelObstacleAvoidanceExp:
      os << "ResidualModelObstacleAvoidanceExp";
      break;
    case ResidualModelTypes::ResidualModelFlyHighSqr:
      os << "ResidualModelFlyHighSqr";
      break;
    case ResidualModelTypes::ResidualModelFlyHighExp:
      os << "ResidualModelFlyHighExp";
      break;
    case ResidualModelTypes::NbResidualModelTypes:
      os << "NbResidualModelTypes";
      break;
    default:
      break;
  }
  return os;
}

ResidualModelFactory::ResidualModelFactory() {}
ResidualModelFactory::~ResidualModelFactory() {}

boost::shared_ptr<crocoddyl::ResidualModelAbstract> ResidualModelFactory::create(
    ResidualModelTypes::Type residual_type, StateModelTypes::Type state_type, std::size_t nu) const {
  StateModelFactory state_factory;
  boost::shared_ptr<crocoddyl::ResidualModelAbstract> residual;
  boost::shared_ptr<crocoddyl::StateMultibody> state =
      boost::static_pointer_cast<crocoddyl::StateMultibody>(state_factory.create(state_type));
  pinocchio::FrameIndex frame_index = state->get_pinocchio()->frames.size() - 1;
  // Fixed frame on lf foot
  // pinocchio::FrameIndex frame_index = state->get_pinocchio()->getFrameId("LF_FOOT");
  pinocchio::SE3 frame_SE3 = pinocchio::SE3::Random();
  // New part for the obstacle avoidance task
  pinocchio::FrameIndex universeId = state->get_pinocchio()->getFrameId("universe");
  boost::shared_ptr<pinocchio::GeometryModel> geomModel(new pinocchio::GeometryModel());

  pinocchio::SE3 sphere_pos(pinocchio::SE3::Identity());
  //sphere_pos.translation() = pinocchio::SE3::LinearType(2., 1., 0.5);
  sphere_pos.translation() = state->get_pinocchio()->frames[frame_index].placement.translation();
  pinocchio::SE3 box_pos(pinocchio::SE3::Identity());
  box_pos.translation() = pinocchio::SE3::LinearType(0., 0., 0);

  pinocchio::GeometryObject sphere("foot", frame_index, state->get_pinocchio()->frames[frame_index].parent, boost::shared_ptr<hpp::fcl::Sphere>(new hpp::fcl::Sphere(0)), sphere_pos);
  pinocchio::GeometryObject box("obstacle", universeId, state->get_pinocchio()->frames[universeId].parent, boost::shared_ptr<hpp::fcl::Box>(new hpp::fcl::Box(0.1, 0.1, 0.1)), box_pos);
  geomModel->addGeometryObject(sphere);
  geomModel->addGeometryObject(box);
  geomModel->addAllCollisionPairs();
  const double beta = 1e-2;

  if (nu == std::numeric_limits<std::size_t>::max()) {
    nu = state->get_nv();
  }
  switch (residual_type) {
    case ResidualModelTypes::ResidualModelState:
      residual = boost::make_shared<crocoddyl::ResidualModelState>(state, state->rand(), nu);
      break;
    case ResidualModelTypes::ResidualModelControl:
      residual = boost::make_shared<crocoddyl::ResidualModelControl>(state, Eigen::VectorXd::Random(nu));
      break;
    case ResidualModelTypes::ResidualModelCoMPosition:
      residual = boost::make_shared<crocoddyl::ResidualModelCoMPosition>(state, Eigen::Vector3d::Random(), nu);
      break;
    case ResidualModelTypes::ResidualModelCentroidalMomentum:
      residual = boost::make_shared<crocoddyl::ResidualModelCentroidalMomentum>(state, Vector6d::Random(), nu);
      break;
    case ResidualModelTypes::ResidualModelFramePlacement:
      residual = boost::make_shared<crocoddyl::ResidualModelFramePlacement>(state, frame_index, frame_SE3, nu);
      break;
    case ResidualModelTypes::ResidualModelFrameRotation:
      residual =
          boost::make_shared<crocoddyl::ResidualModelFrameRotation>(state, frame_index, frame_SE3.rotation(), nu);
      break;
    case ResidualModelTypes::ResidualModelFrameTranslation:
      residual = boost::make_shared<crocoddyl::ResidualModelFrameTranslation>(state, frame_index,
                                                                              frame_SE3.translation(), nu);
      break;
    case ResidualModelTypes::ResidualModelFrameVelocity:
      residual = boost::make_shared<crocoddyl::ResidualModelFrameVelocity>(
          state, frame_index, pinocchio::Motion::Random(), static_cast<pinocchio::ReferenceFrame>(rand() % 2),
          nu);  // the code cannot test LOCAL_WORLD_ALIGNED
      break;
    case ResidualModelTypes::ResidualModelControlGrav:
      residual = boost::make_shared<crocoddyl::ResidualModelControlGrav>(state, nu);
      break;
    case ResidualModelTypes::ResidualModelObstacleAvoidance:
      residual = boost::make_shared<crocoddyl::ResidualModelObstacleAvoidance>(
          state, nu, geomModel, 0, frame_index, pinocchio::LOCAL_WORLD_ALIGNED, beta
          );
      break;
    case ResidualModelTypes::ResidualModelObstacleAvoidanceSqr:
      residual = boost::make_shared<crocoddyl::ResidualModelObstacleAvoidanceSqr>(
          state, nu, geomModel, 0, frame_index, pinocchio::LOCAL_WORLD_ALIGNED, beta
          );
      break;
    case ResidualModelTypes::ResidualModelObstacleAvoidanceExp:
      residual = boost::make_shared<crocoddyl::ResidualModelObstacleAvoidanceExp>(
          state, nu, geomModel, 0, frame_index, pinocchio::LOCAL_WORLD_ALIGNED, beta
          );
    case ResidualModelTypes::ResidualModelFlyHighSqr:
      residual = boost::make_shared<crocoddyl::ResidualModelFlyHighSqr>(
          state, nu, frame_index, beta
          );
      break;
    case ResidualModelTypes::ResidualModelFlyHighExp:
      residual = boost::make_shared<crocoddyl::ResidualModelFlyHighExp>(
          state, nu, frame_index, beta
          );
      break;
    default:
      throw_pretty(__FILE__ ": Wrong ResidualModelTypes::Type given");
      break;
  }
  return residual;
}

boost::shared_ptr<crocoddyl::ResidualModelAbstract> create_random_residual(StateModelTypes::Type state_type) {
  static bool once = true;
  if (once) {
    srand((unsigned)time(NULL));
    once = false;
  }

  ResidualModelFactory factory;
  ResidualModelTypes::Type rand_type =
      static_cast<ResidualModelTypes::Type>(rand() % ResidualModelTypes::NbResidualModelTypes);
  return factory.create(rand_type, state_type);
}

}  // namespace unittest
}  // namespace crocoddyl

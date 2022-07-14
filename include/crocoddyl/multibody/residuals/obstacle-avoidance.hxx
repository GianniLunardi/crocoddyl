// Author: Gianni Lunardi, University of Trento 2022

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/utils/stop-watch.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance.hpp"


namespace crocoddyl {

template <typename Scalar>
ResidualModelObstacleAvoidanceTpl<Scalar>::ResidualModelObstacleAvoidanceTpl(boost::shared_ptr<StateMultibody> state,
                                                                             const std::size_t nu,
                                                                             boost::shared_ptr<GeometryModel> geom_model,
                                                                             const pinocchio::PairIndex pair_id,
                                                                             const pinocchio::FrameIndex frame_id,
                                                                             const pinocchio::ReferenceFrame type,
                                                                             const double beta)
    : Base(state, 2, nu, true, true, false),
      pin_model_(*state->get_pinocchio()),
      geom_model_(geom_model),
      pair_id_(pair_id),
      frame_id_(frame_id),
      type_(type),
      beta_(beta)  {}

template <typename Scalar>
ResidualModelObstacleAvoidanceTpl<Scalar>::~ResidualModelObstacleAvoidanceTpl()  {}

template <typename Scalar>
void ResidualModelObstacleAvoidanceTpl<Scalar>::calc(const boost::shared_ptr<ResidualDataAbstract> &data,
                                                     const Eigen::Ref<const VectorXs> &,
                                                     const Eigen::Ref<const VectorXs> &) {
    START_PROFILER("ResidualModelObstacleAvoidance::calc");
    Data* d = static_cast<Data*>(data.get());

    // Compute the distance for the collision pair
    START_PROFILER("ResidualModelObstacleAvoidance::calc.update");
    pinocchio::updateGeometryPlacements(pin_model_, *d->pinocchio, *geom_model_.get(), d->geometry);
    STOP_PROFILER("ResidualModelObstacleAvoidance::calc.update");

    START_PROFILER("ResidualModelObstacleAvoidance::calc.computeDistance");
    pinocchio::computeDistance(*geom_model_.get(), d->geometry, pair_id_);
    STOP_PROFILER("ResidualModelObstacleAvoidance::calc.computeDistance");
    // Compute the velocity of the foot
    d->v = (pinocchio::getFrameVelocity(pin_model_, *d->pinocchio, frame_id_, type_)).toVector();
    // Compute the residual
    if(d->geometry.distanceResults[pair_id_].min_distance > 0) {
        // Square root of the distance btw the obstacles
        d->dist_sqrt = std::sqrt(d->geometry.distanceResults[pair_id_].min_distance);
        d->r[0] = d->v[0] / (d->dist_sqrt + beta_);
        d->r[1] = d->v[1] / (d->dist_sqrt + beta_);
    }
    else {
        d->r[0] = d->v[0] / beta_;
        d->r[1] = d->v[1] / beta_;
    }
    STOP_PROFILER("ResidualModelObstacleAvoidance::calc");
}

template <typename Scalar>
void ResidualModelObstacleAvoidanceTpl<Scalar>::calcDiff(const boost::shared_ptr<ResidualDataAbstract> &data,
                                                         const Eigen::Ref<const VectorXs> &,
                                                         const Eigen::Ref<const VectorXs> &) {
    START_PROFILER("ResidualModelObstacleAvoidance::calcDiff");
    Data* d = static_cast<Data*>(data.get());
    // const Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    START_PROFILER("ResidualModelObstacleAvoidance::calcDiff.getFrame");
    // Compute the frame Jacobian, consider the first three rows (derivatives wrt the position of the frame)
    pinocchio::getFrameJacobian(pin_model_, *d->pinocchio, frame_id_, type_, d->J);

    // Compute the derivatives of the foot's velocity
    pinocchio::getFrameVelocityDerivatives(pin_model_, *d->pinocchio, frame_id_, type_,
                                           d->dv_dx.leftCols(nv_), d->dv_dx.rightCols(nv_));
    STOP_PROFILER("ResidualModelObstacleAvoidance::calcDiff.getFrame");
    Vector3s omega = d->v.tail(3);
    // Add the missing term in the top left corner (LOCAL_WORLD_ALIGNED)
    d->dv_dx.topLeftCorner(3, nv_).noalias() += pinocchio::skew(omega) * d->J.topRows(3);
    // Compute the residual Jacobian, considering the two cases
    if(d->geometry.distanceResults[pair_id_].min_distance > 0) {
        // 1) min dist > 0 --> compute the complete derivatives of the residual
        // Normal: take the opposite direction wrt the one compute by hpp-fcl
        d->dist_dx.head(nv_) = -d->geometry.distanceResults[pair_id_].normal.transpose() * d->J.topRows(3);
        // Chain rule
        d->Rx = 1 / (d->dist_sqrt + beta_) * d->dv_dx.topRows(2)
                - 1 / (2 * d->dist_sqrt * (d->dist_sqrt + beta_) * (d->dist_sqrt + beta_))* d->v.head(2) * d->dist_dx.transpose();
    }
    else {
        d->Rx = 1 / beta_ * d->dv_dx.topRows(2);
    }
    STOP_PROFILER("ResidualModelObstacleAvoidance::calcDiff");
}

template <typename Scalar>
boost::shared_ptr<ResidualDataAbstractTpl<Scalar> > ResidualModelObstacleAvoidanceTpl<Scalar>::createData(
        DataCollectorAbstract *const data) {
    return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this, data);
}

template <typename Scalar>
const pinocchio::GeometryModel &ResidualModelObstacleAvoidanceTpl<Scalar>::get_geometry() const {
    return *geom_model_.get();
}

template <typename Scalar>
pinocchio::PairIndex ResidualModelObstacleAvoidanceTpl<Scalar>::get_pair_id() const {
    return pair_id_;
}

template <typename Scalar>
pinocchio::FrameIndex ResidualModelObstacleAvoidanceTpl<Scalar>::get_frame_id() const {
    return frame_id_;
}

template <typename Scalar>
pinocchio::ReferenceFrame ResidualModelObstacleAvoidanceTpl<Scalar>::get_type() const {
    return type_;
}

template <typename Scalar>
void ResidualModelObstacleAvoidanceTpl<Scalar>::set_type(const pinocchio::ReferenceFrame type) {
    type_ = type;
}


}   // namespace crocoddyl

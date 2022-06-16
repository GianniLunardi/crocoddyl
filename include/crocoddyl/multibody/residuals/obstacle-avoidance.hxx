// Author: Gianni Lunardi, University of Trento 2022

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
//#include <armadillo>
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
                                                     const Eigen::Ref<const VectorXs> &x,
                                                     const Eigen::Ref<const VectorXs> &) {
    START_PROFILER("ResidualModelObstacleAvoidance::calc");
    Data* d = static_cast<Data*>(data.get());

    // Compute the distance for the collision pair
    pinocchio::updateGeometryPlacements(pin_model_, *d->pinocchio, *geom_model_.get(), d->geometry, x.head(nq_));
    pinocchio::computeDistance(*geom_model_.get(), d->geometry, pair_id_);
    // Compute the distance between the two nearest points
    d->p_diff = d->geometry.distanceResults[pair_id_].nearest_points[0] -
                d->geometry.distanceResults[pair_id_].nearest_points[1];
    d->dist = d->p_diff.norm();
    // Save the square root of the distance, to avoid its computation again
    d->dist_sqrt = std::sqrt(d->dist);
    // Compute the velocity of the foot
    d->v = (pinocchio::getFrameVelocity(pin_model_, *d->pinocchio, frame_id_, type_)).toVector();
    // Compute the residual
    if(d->geometry.distanceResults[pair_id_].min_distance > 0) {
        d->r[0] = d->v[0] * d->v[0] / (d->dist_sqrt + beta_);
        d->r[1] = d->v[1] * d->v[1] / (d->dist_sqrt + beta_);
    }
    else {
        d->r[0] = d->v[0] * d->v[0] / beta_;
        d->r[1] = d->v[1] * d->v[1] / beta_;
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

    // IMPORTANT: translate the derivatives of the linear velocity from the WORLD to LOCAL_WORLD_ALIGNED reference frame
    pinocchio::getFrameVelocityDerivatives(pin_model_, *d->pinocchio, frame_id_, pinocchio::WORLD,
                                           d->dv0_dx.leftCols(nv_), d->dv0_dx.rightCols(nv_));
    STOP_PROFILER("ResidualModelObstacleAvoidance::calcDiff.getFrame");
    Vector3s p = d->pinocchio->oMf[frame_id_].translation();
    Vector3s omega = d->v.tail(3);
    d->dv_dx.topLeftCorner(3, nv_) = d->dv0_dx.topLeftCorner(3, nv_);
    d->dv_dx.topLeftCorner(3, nv_).noalias() -= pinocchio::skew(p) * d->dv0_dx.bottomLeftCorner(3, nv_);
    d->dv_dx.topLeftCorner(3, nv_).noalias() += pinocchio::skew(omega) * d->J.topRows(3);
    // Compute the residual Jacobian, considering the two cases
    if(d->geometry.distanceResults[pair_id_].min_distance > 0) {
        // 1) min dist > 0 --> compute the complete derivatives of the residual
        Vector3s norm_der = d->p_diff / d->dist;
        d->dist_der.head(nv_) = norm_der.transpose() * d->J.topRows(3);
        Vector2s vel_square;
        vel_square << d->v[0] * d->v[0], d->v[1] * d->v[1];
        // Chain rule
        d->Rx = 2 / (d->dist_sqrt + beta_) * d->v.head(2).asDiagonal() * d->dv_dx.topRows(2)
                - 1 / (2 * d->dist_sqrt * (d->dist_sqrt + beta_) * (d->dist_sqrt + beta_))* vel_square * d->dist_der.transpose();
    }
    else {
        d->Rx = 2 / beta_ * d->v.head(2).asDiagonal() * d->dv_dx.topRows(2);
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

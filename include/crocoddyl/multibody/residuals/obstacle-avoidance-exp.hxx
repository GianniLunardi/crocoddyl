// Author: Gianni Lunardi, University of Trento 2022

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/utils/stop-watch.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance-exp.hpp"


namespace crocoddyl {

template <typename Scalar>
ResidualModelObstacleAvoidanceExpTpl<Scalar>::ResidualModelObstacleAvoidanceExpTpl(boost::shared_ptr<StateMultibody> state,
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
ResidualModelObstacleAvoidanceExpTpl<Scalar>::~ResidualModelObstacleAvoidanceExpTpl()  {}

template <typename Scalar>
void ResidualModelObstacleAvoidanceExpTpl<Scalar>::calc(const boost::shared_ptr<ResidualDataAbstract> &data,
                                                        const Eigen::Ref<const VectorXs> &,
                                                        const Eigen::Ref<const VectorXs> &) {
    START_PROFILER("ResidualModelObstacleAvoidanceExp::calc");
    Data* d = static_cast<Data*>(data.get());

    // Compute the distance for the collision pair
    pinocchio::updateGeometryPlacements(pin_model_, *d->pinocchio, *geom_model_.get(), d->geometry);
    pinocchio::computeDistance(*geom_model_.get(), d->geometry, pair_id_);
    d->dist_exp = std::exp(beta_ * d->geometry.distanceResults[pair_id_].min_distance);
    // Compute the velocity of the foot
    d->v = (pinocchio::getFrameVelocity(pin_model_, *d->pinocchio, frame_id_, type_)).toVector();
    // Compute the residual
    d->r[0] = d->v[0] / d->dist_exp;
    d->r[1] = d->v[1] / d->dist_exp;
    STOP_PROFILER("ResidualModelObstacleAvoidanceExp::calc");
}

template <typename Scalar>
void ResidualModelObstacleAvoidanceExpTpl<Scalar>::calcDiff(const boost::shared_ptr<ResidualDataAbstract> &data,
                                                            const Eigen::Ref<const VectorXs> &,
                                                            const Eigen::Ref<const VectorXs> &) {
    START_PROFILER("ResidualModelObstacleAvoidanceExp::calcDiff");
    Data* d = static_cast<Data*>(data.get());
    const Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

    // Compute the frame Jacobian, consider the first three rows (derivatives wrt the position of the frame)
    pinocchio::getFrameJacobian(pin_model_, *d->pinocchio, frame_id_, type_, d->J);
    // Compute the derivatives of the foot's velocity
    pinocchio::getFrameVelocityDerivatives(pin_model_, *d->pinocchio, frame_id_, type_,
                                           d->dv_dx.leftCols(nv_), d->dv_dx.rightCols(nv_));
    Vector3s omega = d->v.tail(3);
    // Add the missing term in the top left corner (LOCAL_WORLD_ALIGNED)
    d->dv_dx.topLeftCorner(3, nv_).noalias() += pinocchio::skew(omega) * d->J.topRows(3);
    // Normal: take the opposite direction wrt the one compute by hpp-fcl
    d->dist_dx.head(nv_) = -d->geometry.distanceResults[pair_id_].normal.transpose() * d->J.topRows(3);
    // Compute the residual Jacobian
    d->Rx = 1 / d->dist_exp * d->dv_dx.topRows(2) - beta_ * d->r * d->dist_dx.transpose();
    STOP_PROFILER("ResidualModelObstacleAvoidanceExp::calcDiff");
}

template <typename Scalar>
boost::shared_ptr<ResidualDataAbstractTpl<Scalar> > ResidualModelObstacleAvoidanceExpTpl<Scalar>::createData(
        DataCollectorAbstract *const data) {
    return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this, data);
}

template <typename Scalar>
const pinocchio::GeometryModel &ResidualModelObstacleAvoidanceExpTpl<Scalar>::get_geometry() const {
    return *geom_model_.get();
}

template <typename Scalar>
pinocchio::PairIndex ResidualModelObstacleAvoidanceExpTpl<Scalar>::get_pair_id() const {
    return pair_id_;
}

template <typename Scalar>
pinocchio::FrameIndex ResidualModelObstacleAvoidanceExpTpl<Scalar>::get_frame_id() const {
    return frame_id_;
}

template <typename Scalar>
pinocchio::ReferenceFrame ResidualModelObstacleAvoidanceExpTpl<Scalar>::get_type() const {
    return type_;
}

template <typename Scalar>
void ResidualModelObstacleAvoidanceExpTpl<Scalar>::set_type(const pinocchio::ReferenceFrame type) {
    type_ = type;
}


}   // namespace crocoddyl
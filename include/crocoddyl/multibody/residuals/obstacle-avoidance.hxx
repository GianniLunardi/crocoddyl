// Author: Gianni Lunardi, University of Trento 2022

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <armadillo>
#include "crocoddyl/core/utils/exception.hpp"
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
    Data* d = static_cast<Data*>(data.get());
    const Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

    pinocchio::forwardKinematics(pin_model_, *d->pinocchio, x.head(nq_), x.tail(nv_));
    // Compute the distance for the collision pair
    pinocchio::updateGeometryPlacements(pin_model_, *d->pinocchio, *geom_model_.get(), d->geometry, x.head(nq_));
    pinocchio::computeDistance(*geom_model_.get(), d->geometry, pair_id_);
    // Compute the distance between the two nearest points
    d->p_diff = d->geometry.distanceResults[pair_id_].nearest_points[0] -
                d->geometry.distanceResults[pair_id_].nearest_points[1];
    d->dist = d->p_diff.norm();
    //std::cout << "p_diff = \n" << d->p_diff.format(HeavyFmt) << std::endl;
    // Compute the velocity of the foot
    d->v = (pinocchio::getFrameVelocity(pin_model_, *d->pinocchio, frame_id_, type_)).toVector();
    // Compute the residual
    if(d->geometry.distanceResults[pair_id_].min_distance > 0) {
        d->r[0] = std::pow(d->v[0], 2) / (std::sqrt(d->dist) + beta_);
        d->r[1] = std::pow(d->v[1], 2) / (std::sqrt(d->dist) + beta_);
    }
    else {
        d->r[0] = std::pow(d->v[0], 2) / beta_;
        d->r[1] = std::pow(d->v[1], 2) / beta_;
    }
    //std::cout << "r = \n" << d->r.format(HeavyFmt) << std::endl;
}

template <typename Scalar>
void ResidualModelObstacleAvoidanceTpl<Scalar>::calcDiff(const boost::shared_ptr<ResidualDataAbstract> &data,
                                                         const Eigen::Ref<const VectorXs> &,
                                                         const Eigen::Ref<const VectorXs> &) {
    Data* d = static_cast<Data*>(data.get());
    const Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

    // Compute the frame Jacobian, consider the first three rows (derivatives wrt the position of the frame)
    pinocchio::getFrameJacobian(pin_model_, *d->pinocchio, frame_id_, type_, d->J);
    Matrix3xs p_diff_der = d->J.template topRows<3>();
    //std::cout << "p_diff_der = \n" << p_diff_der.format(HeavyFmt) << std::endl;

    // Compute the derivatives of the foot's velocity
    pinocchio::getFrameVelocityDerivatives(pin_model_, *d->pinocchio, frame_id_, type_,
                                           d->dv_dx.leftCols(nv_), d->dv_dx.rightCols(nv_));

    // IMPORTANT: translate the derivatives of the linear velocity from the WORLD to LOCAL_WORLD_ALIGNED reference frame
    Matrix6xs dv0_dx(Matrix6xs::Zero(6, state_->get_ndx()));
    pinocchio::getFrameVelocityDerivatives(pin_model_, *d->pinocchio, frame_id_, pinocchio::WORLD,
                                           dv0_dx.leftCols(nv_), dv0_dx.rightCols(nv_));
    //std::cout << "dv0_dx = \n" << dv0_dx.format(HeavyFmt) << std::endl;

    Vector3s p = d->pinocchio->oMf[frame_id_].translation();
    Vector3s omega = d->v.tail(3);
    //std::cout << "p = " << p.format(HeavyFmt) << ", omega = " << omega.format(HeavyFmt) << std::endl;
    Matrix3s p_skew = pinocchio::skew(p);
    Matrix3s omega_skew = pinocchio::skew(omega);
//    std::cout << "prod 1:\n" << (p_skew * dv0_dx.bottomLeftCorner(3, nv_)).format(HeavyFmt) << "\n" <<
//                 "prod 2:\n" << (omega_skew * p_diff_der).format(HeavyFmt) << std::endl;
    d->dv_dx.topLeftCorner(3, nv_) = dv0_dx.topLeftCorner(3, nv_) - p_skew * dv0_dx.bottomLeftCorner(3, nv_) +
                                     omega_skew * p_diff_der;
    //std::cout << "dv_dx = \n" << d->dv_dx.format(HeavyFmt) << std::endl;
    // Compute the residual Jacobian, considering the two cases
    if(d->geometry.distanceResults[pair_id_].min_distance > 0) {
        // 1) min dist > 0 --> compute the complete derivatives of the residual
        VectorXs norm_der = d->p_diff / d->dist;
        //std::cout << "norm_der = \n" << norm_der.format(HeavyFmt) << std::endl;
        VectorXs dist_der_dq = norm_der.transpose() * p_diff_der;
        VectorXs dist_der_dv(VectorXs::Zero(nv_, 1));
        VectorXs dist_der(dist_der_dq.size() + dist_der_dv.size());
        dist_der << dist_der_dq, dist_der_dv;
        //std::cout << "dist_der = \n" << dist_der.format(HeavyFmt) << std::endl;

        // Chain rule
        for(size_t i = 0; i < nr_; i++) {
            for(size_t j = 0; j < ndx_; j++) {
                d->Rx(i,j) = - (std::pow(d->v[i], 2) * dist_der[j]) / (2 * std::pow((beta_ + std::sqrt(d->dist)), 2) * std::sqrt(d->dist)) +
                                (2 * d->v[i] * d->dv_dx(i,j)) / (beta_ + std::sqrt(d->dist));
            }
        }
    }
    else {
        for(size_t i = 0; i < nr_; i++) {
            for(size_t j = 0; j < ndx_; j++) {
                d->Rx(i,j) = (2 * d->v[i] * d->dv_dx(i,j)) / beta_;
            }
        }
    }
    //std::cout << "Rx = \n" << d->Rx.format(HeavyFmt) << std::endl;
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

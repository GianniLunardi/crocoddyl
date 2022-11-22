// Author: Gianni Lunardi, University of Trento 2022

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/utils/stop-watch.hpp"
#include "crocoddyl/multibody/residuals/fly-high-sqr.hpp"


namespace crocoddyl {

template <typename Scalar>
ResidualModelFlyHighSqrTpl<Scalar>::ResidualModelFlyHighSqrTpl(boost::shared_ptr<StateMultibody> state,
                                                               const std::size_t nu,
                                                               const pinocchio::FrameIndex frame_id,
                                                               const double beta)
    : Base(state, 2, nu, true, true, false),
      pin_model_(*state->get_pinocchio()),
      frame_id_(frame_id),
      beta_(beta) {}

template <typename Scalar>
ResidualModelFlyHighSqrTpl<Scalar>::~ResidualModelFlyHighSqrTpl() {}

template <typename Scalar>
void ResidualModelFlyHighSqrTpl<Scalar>::calc(const boost::shared_ptr<ResidualDataAbstract> &data,
                                              const Eigen::Ref<const VectorXs> &,
                                              const Eigen::Ref<const VectorXs> &) {
    START_PROFILER("ResidualModelFlyHighSqr::calc");
    Data* d = static_cast<Data*>(data.get());

    pinocchio::updateFramePlacement(pin_model_, *d->pinocchio, frame_id_);
    // Compute the velocity of the foot in the LWA reference frame
    d->v = pinocchio::getFrameVelocity(pin_model_,
                                       *d->pinocchio,
                                       frame_id_,
                                       pinocchio::LOCAL_WORLD_ALIGNED).toVector();

    // Compute the height from the ground
    d->h = d->pinocchio->oMf[frame_id_].translation()[2];
    if(d->h > 0) {
        d->h_sqrt = std::sqrt(d->h);
        d->r[0] = d->v[0] * d->v[0] / (d->h_sqrt + beta_);
        d->r[1] = d->v[1] * d->v[1] / (d->h_sqrt + beta_);
    }
    else {
        d->r[0] = d->v[0] * d->v[0] / beta_;
        d->r[1] = d->v[1] * d->v[1] / beta_;
    }
    STOP_PROFILER("ResidualModelFlyHighSqr::calc");
}

template <typename Scalar>
void ResidualModelFlyHighSqrTpl<Scalar>::calcDiff(const boost::shared_ptr<ResidualDataAbstract> &data,
                                                  const Eigen::Ref<const VectorXs> &,
                                                  const Eigen::Ref<const VectorXs> &) {
    START_PROFILER("ResidualModelFlyHighSqr::calcDiff");
    Data* d = static_cast<Data*>(data.get());
    // const Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

    pinocchio::getFrameVelocityDerivatives(pin_model_, *d->pinocchio, frame_id_, pinocchio::LOCAL,
                                           d->l_dv_dx.leftCols(nv_), d->l_dv_dx.rightCols(nv_));
    Vector3s v = pinocchio::getFrameVelocity(pin_model_, *d->pinocchio,
                                             frame_id_, pinocchio::LOCAL).linear();
    Matrix3s R = d->pinocchio->oMf[frame_id_].rotation();
    // Compute the derivatives of the velocity in the LWA frame
    Matrix3xs vxJ = pinocchio::skew(-v) * d->l_dv_dx.bottomRightCorner(3, nv_);
    vxJ += d->l_dv_dx.topLeftCorner(3, nv_);

    d->dv_dx.leftCols(nv_) = R * vxJ;
    d->dv_dx.rightCols(nv_) = R * d->l_dv_dx.topRightCorner(3, nv_);
    // Compute the residual Jacobian, considering the two cases
    d->Rx = 2 / (d->h_sqrt + beta_) * d->v.head(2).asDiagonal() * d->dv_dx.topRows(2);
    if(d->h > 0) {
        Vector2s vel_square;
        vel_square << d->v[0] * d->v[0], d->v[1] * d->v[1];
        d->Rx.leftCols(nv_) += - 1 / (2 * d->h_sqrt * (d->h_sqrt + beta_) * (d->h_sqrt + beta_)) * vel_square *
                               d->dv_dx.rightCols(nv_).row(2);
    }
    STOP_PROFILER("ResidualModelFlyHighSqr::calcDiff");
}

template <typename Scalar>
boost::shared_ptr<ResidualDataAbstractTpl<Scalar> > ResidualModelFlyHighSqrTpl<Scalar>::createData(
        DataCollectorAbstract *const data) {
    return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this, data);
}

template <typename Scalar>
pinocchio::FrameIndex ResidualModelFlyHighSqrTpl<Scalar>::get_frame_id() const {
    return frame_id_;
}

template <typename Scalar>
Scalar ResidualModelFlyHighSqrTpl<Scalar>::get_beta() const {
    return beta_;
}


}   // namespace crocoddyl
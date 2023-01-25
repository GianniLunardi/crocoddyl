///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2023, University of Edinburgh, LAAS-CNRS
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/numdiff/actuation.hpp"

namespace crocoddyl {

template <typename Scalar>
ActuationModelNumDiffTpl<Scalar>::ActuationModelNumDiffTpl(boost::shared_ptr<Base> model)
    : Base(model->get_state(), model->get_nu()),
      model_(model),
      e_jac_(std::sqrt(2.0 * std::numeric_limits<Scalar>::epsilon())) {}

template <typename Scalar>
ActuationModelNumDiffTpl<Scalar>::~ActuationModelNumDiffTpl() {}

template <typename Scalar>
void ActuationModelNumDiffTpl<Scalar>::calc(const boost::shared_ptr<ActuationDataAbstract>& data,
                                            const Eigen::Ref<const VectorXs>& x, const Eigen::Ref<const VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != model_->get_state()->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(model_->get_state()->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }
  Data* d = static_cast<Data*>(data.get());
  model_->calc(d->data_0, x, u);
  data->tau = d->data_0->tau;
}

template <typename Scalar>
void ActuationModelNumDiffTpl<Scalar>::calc(const boost::shared_ptr<ActuationDataAbstract>& data,
                                            const Eigen::Ref<const VectorXs>& x) {
  if (static_cast<std::size_t>(x.size()) != model_->get_state()->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(model_->get_state()->get_nx()) + ")");
  }
  Data* d = static_cast<Data*>(data.get());
  model_->calc(d->data_0, x);
  data->tau = d->data_0->tau;
}

template <typename Scalar>
void ActuationModelNumDiffTpl<Scalar>::calcDiff(const boost::shared_ptr<ActuationDataAbstract>& data,
                                                const Eigen::Ref<const VectorXs>& x,
                                                const Eigen::Ref<const VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != model_->get_state()->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(model_->get_state()->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }
  Data* d = static_cast<Data*>(data.get());
  const VectorXs& tau0 = d->data_0->tau;
  d->du.setZero();

  // Computing the d actuation(x,u) / dx
  model_->get_state()->diff(model_->get_state()->zero(), x, d->dx);
  d->x_norm = d->dx.norm();
  d->dx.setZero();
  d->xh_jac = e_jac_ * std::max(1., d->x_norm);
  for (std::size_t ix = 0; ix < model_->get_state()->get_ndx(); ++ix) {
    d->dx(ix) = d->xh_jac;
    model_->get_state()->integrate(x, d->dx, d->xp);
    model_->calc(d->data_x[ix], d->xp, u);
    d->dtau_dx.col(ix) = (d->data_x[ix]->tau - tau0) / d->xh_jac;
    d->dx(ix) = 0.;
  }

  // Computing the d actuation(x,u) / du
  d->uh_jac = e_jac_ * std::max(1., u.norm());
  for (unsigned iu = 0; iu < model_->get_nu(); ++iu) {
    d->du(iu) = d->uh_jac;
    model_->calc(d->data_u[iu], x, u + d->du);
    d->dtau_du.col(iu) = (d->data_u[iu]->tau - tau0) / d->uh_jac;
    d->du(iu) = 0.;
  }
}

template <typename Scalar>
void ActuationModelNumDiffTpl<Scalar>::calcDiff(const boost::shared_ptr<ActuationDataAbstract>& data,
                                                const Eigen::Ref<const VectorXs>& x) {
  if (static_cast<std::size_t>(x.size()) != model_->get_state()->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(model_->get_state()->get_nx()) + ")");
  }
  Data* d = static_cast<Data*>(data.get());
  const VectorXs& tau0 = d->data_0->tau;
  d->dx.setZero();

  // Computing the d actuation(x,u) / dx
  model_->get_state()->diff(model_->get_state()->zero(), x, d->dx);
  d->x_norm = d->dx.norm();
  d->dx.setZero();
  d->xh_jac = e_jac_ * std::max(1., d->x_norm);
  for (std::size_t ix = 0; ix < model_->get_state()->get_ndx(); ++ix) {
    d->dx(ix) = d->xh_jac;
    model_->get_state()->integrate(x, d->dx, d->xp);
    model_->calc(d->data_x[ix], d->xp);
    d->dtau_dx.col(ix) = (d->data_x[ix]->tau - tau0) / d->xh_jac;
    d->dx(ix) = 0.;
  }
}

template <typename Scalar>
void ActuationModelNumDiffTpl<Scalar>::commands(const boost::shared_ptr<ActuationDataAbstract>& data,
                                                const Eigen::Ref<const VectorXs>& x,
                                                const Eigen::Ref<const VectorXs>& tau) {
  if (static_cast<std::size_t>(x.size()) != model_->get_state()->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(model_->get_state()->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(tau.size()) != model_->get_state()->get_nv()) {
    throw_pretty("Invalid argument: "
                 << "tau has wrong dimension (it should be " + std::to_string(model_->get_state()->get_nv()) + ")");
  }
  Data* d = static_cast<Data*>(data.get());

  model_->torqueTransform(d->data_x[0], x, tau);
  data->u.noalias() = d->data_x[0]->Mtau * tau;
}

template <typename Scalar>
boost::shared_ptr<ActuationDataAbstractTpl<Scalar> > ActuationModelNumDiffTpl<Scalar>::createData() {
  return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this);
}

template <typename Scalar>
const boost::shared_ptr<ActuationModelAbstractTpl<Scalar> >& ActuationModelNumDiffTpl<Scalar>::get_model() const {
  return model_;
}

template <typename Scalar>
const Scalar ActuationModelNumDiffTpl<Scalar>::get_disturbance() const {
  return e_jac_;
}

template <typename Scalar>
void ActuationModelNumDiffTpl<Scalar>::set_disturbance(const Scalar disturbance) {
  if (disturbance < 0.) {
    throw_pretty("Invalid argument: "
                 << "Disturbance constant is positive");
  }
  e_jac_ = disturbance;
}

}  // namespace crocoddyl

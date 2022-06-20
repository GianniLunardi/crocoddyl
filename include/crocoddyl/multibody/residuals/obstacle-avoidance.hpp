// Author: Gianni Lunardi, University of Trento 2022

#ifndef CROCODDYL_MULTIBODY_RESIDUALS_OBSTACLE_AVOIDANCE_HPP_
#define CROCODDYL_MULTIBODY_RESIDUALS_OBSTACLE_AVOIDANCE_HPP_

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/motion.hpp>
#include <pinocchio/multibody/fcl.hpp>
#include <pinocchio/multibody/geometry.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/geometry.hpp>
#include <Eigen/src/Core/DenseBase.h>

#include "crocoddyl/multibody/fwd.hpp"
#include "crocoddyl/core/residual-base.hpp"
#include "crocoddyl/multibody/states/multibody.hpp"
#include "crocoddyl/multibody/data/multibody.hpp"
#include "crocoddyl/core/utils/exception.hpp"

namespace crocoddyl {

/**

BRIEF FOR DOXYGEN

*/


template <typename _Scalar>
class ResidualModelObstacleAvoidanceTpl: public ResidualModelAbstractTpl<_Scalar> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    typedef MathBaseTpl<Scalar> MathBase;
    typedef ResidualModelAbstractTpl<Scalar> Base;
    typedef ResidualDataObstacleAvoidanceTpl<Scalar> Data;
    typedef StateMultibodyTpl<Scalar> StateMultibody;
    typedef ResidualDataAbstractTpl<Scalar> ResidualDataAbstract;
    typedef DataCollectorAbstractTpl<Scalar> DataCollectorAbstract;

    typedef pinocchio::GeometryModel GeometryModel;
    typedef typename MathBase::Vector2s Vector2s;
    typedef typename MathBase::Vector3s Vector3s;
    typedef typename MathBase::VectorXs VectorXs;
    typedef typename MathBase::Matrix3s Matrix3s;
    typedef typename MathBase::Matrix6xs Matrix6xs;

    ResidualModelObstacleAvoidanceTpl(boost::shared_ptr<StateMultibody> state,
                                      const std::size_t nu,
                                      boost::shared_ptr<GeometryModel> geom_model,
                                      const pinocchio::PairIndex pair_id,
                                      const pinocchio::FrameIndex frame_id,
                                      const pinocchio::ReferenceFrame type,
                                      const double beta);

    /*
    Brief
    */

    virtual ~ResidualModelObstacleAvoidanceTpl();

    virtual void calc(const boost::shared_ptr<ResidualDataAbstract> &data,
                      const Eigen::Ref<const VectorXs> &x, const Eigen::Ref<const VectorXs> &u);

    virtual void calcDiff(const boost::shared_ptr<ResidualDataAbstract> &data,
                          const Eigen::Ref<const VectorXs> &x, const Eigen::Ref<const VectorXs> &u);

    virtual boost::shared_ptr<ResidualDataAbstract> createData(DataCollectorAbstract *const data);

    const pinocchio::GeometryModel &get_geometry() const;

    pinocchio::PairIndex get_pair_id() const;

    pinocchio::FrameIndex get_frame_id() const;

    pinocchio::ReferenceFrame get_type() const;

    void set_type(const pinocchio::ReferenceFrame type);

  protected:
    using Base::nu_;
    using Base::nr_;
    using Base::state_;
    using Base::q_dependent_;
    using Base::v_dependent_;
    using Base::unone_;
    const std::size_t nq_ = state_->get_nq();
    const std::size_t nv_ = state_->get_nv();
    const std::size_t ndx_ = state_->get_ndx();


  private:
    typename StateMultibody::PinocchioModel pin_model_;
    boost::shared_ptr<pinocchio::GeometryModel> geom_model_;
    pinocchio::PairIndex pair_id_;
    pinocchio::FrameIndex frame_id_;
    pinocchio::ReferenceFrame type_;
    const double beta_;

};

template <typename _Scalar>
struct ResidualDataObstacleAvoidanceTpl : public ResidualDataAbstractTpl<_Scalar> {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    typedef MathBaseTpl<Scalar> MathBase;
    typedef ResidualDataAbstractTpl<Scalar> Base;
    typedef StateMultibodyTpl<Scalar> StateMultibody;
    typedef DataCollectorAbstractTpl<Scalar> DataCollectorAbstract;

    typedef typename MathBase::MatrixXs MatrixXs;
    typedef typename MathBase::Matrix6xs Matrix6xs;
    typedef typename MathBase::Vector3s Vector3s;
    typedef typename MathBase::Vector6s Vector6s;
    typedef typename MathBase::VectorXs VectorXs;

    template <template <typename Scalar> class Model>
    ResidualDataObstacleAvoidanceTpl(Model <Scalar> *const model, DataCollectorAbstract *const data)
        :   Base(model, data),
            geometry(pinocchio::GeometryData(model->get_geometry())),
            dist_der(VectorXs::Zero(model->get_state()->get_ndx())),
            J(Matrix6xs::Zero(6, model->get_state()->get_nv())),
            dv_dx(Matrix6xs::Zero(6, model->get_state()->get_ndx())) {
        // Check that proper shared data has been passed
        DataCollectorMultibodyTpl<Scalar> *d = dynamic_cast<DataCollectorMultibodyTpl<Scalar> *>(shared);
        if (d == NULL) {
            throw_pretty("Invalid argument: the shared data should be derived from DataCollectorMultibodyTpl");
        }

        // Avoids data casting at runtime
        pinocchio = d->pinocchio;
    }

    pinocchio::GeometryData geometry;
    pinocchio::DataTpl<Scalar>* pinocchio;
    Scalar dist_sqrt;
    Vector6s v;
    VectorXs dist_der;
    Matrix6xs J;
    Matrix6xs dv_dx;
    using Base::r;
    using Base::Ru;
    using Base::Rx;
    using Base::shared;
};

}   // namespace crocoddyl

#include "crocoddyl/multibody/residuals/obstacle-avoidance.hxx"

#endif   // CROCODDYL_MULTIBODY_RESIDUALS_OBSTACLE_AVOIDANCE_HPP_
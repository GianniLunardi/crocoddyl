// Author: Gianni Lunardi, University of Trento 2022

#include "python/crocoddyl/multibody/multibody.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance-sqr.hpp"

namespace crocoddyl {
namespace python {

void exposeResidualObstacleSqr() {
    bp::register_ptr_to_python<boost::shared_ptr<ResidualModelObstacleAvoidanceSqr> >();

    bp::class_<ResidualModelObstacleAvoidanceSqr, bp::bases<ResidualModelAbstract> >(
        "ResidualModelObstacleSqr",
        bp::init<
            boost::shared_ptr<StateMultibody>,
            std::size_t,
            boost::shared_ptr<pinocchio::GeometryModel>,
            pinocchio::PairIndex,
            pinocchio::FrameIndex,
            pinocchio::ReferenceFrame,
            double
        >( bp::args("self", "state", "nu", "geom_model", "pair_id", "frame_id", "type", "beta") )
    )
        .def<void (ResidualModelObstacleAvoidanceSqr::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calc", &ResidualModelObstacleAvoidanceSqr::calc, bp::args("self", "data", "x", "u") )
        .def<void (ResidualModelObstacleAvoidanceSqr::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calcDiff", &ResidualModelObstacleAvoidanceSqr::calcDiff, bp::args("self", "data", "x", "u") )
        .def("createData", &ResidualModelObstacleAvoidanceSqr::createData, bp::with_custodian_and_ward_postcall<0, 2>(),
            bp::args("self", "collector") )
        .add_property("pair_id", &ResidualModelObstacleAvoidanceSqr::get_pair_id)
        .add_property("frame_id", &ResidualModelObstacleAvoidanceSqr::get_frame_id)
        .add_property("type", &ResidualModelObstacleAvoidanceSqr::get_type, &ResidualModelObstacleAvoidanceSqr::set_type);

    bp::register_ptr_to_python<boost::shared_ptr<ResidualDataObstacleAvoidanceSqr> >();

    bp::class_<ResidualDataObstacleAvoidanceSqr, bp::bases<ResidualDataAbstract> >(
        "ResidualDataObstacleSqr",
        bp::init<ResidualModelObstacleAvoidanceSqr*, DataCollectorAbstract*>(
            bp::args("self", "model", "collector")
            )[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 3> >()]
    )
        .add_property("pinocchio",
                      bp::make_getter(&ResidualDataObstacleAvoidanceSqr::pinocchio, bp::return_internal_reference<>() ))
        .add_property("geometry",
                      bp::make_getter(&ResidualDataObstacleAvoidanceSqr::geometry, bp::return_internal_reference<>() ));
}

}   // namespace python
}   // namespace crocoddyl
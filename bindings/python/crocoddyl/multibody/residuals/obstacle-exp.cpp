// Author: Gianni Lunardi, University of Trento 2022

#include "python/crocoddyl/multibody/multibody.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance-exp.hpp"

namespace crocoddyl {
namespace python {

void exposeResidualObstacleExp() {
    bp::register_ptr_to_python<boost::shared_ptr<ResidualModelObstacleAvoidanceExp> >();

    bp::class_<ResidualModelObstacleAvoidanceExp, bp::bases<ResidualModelAbstract> >(
        "ResidualModelObstacleExp",
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
        .def<void (ResidualModelObstacleAvoidanceExp::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calc", &ResidualModelObstacleAvoidanceExp::calc, bp::args("self", "data", "x", "u") )
        .def<void (ResidualModelObstacleAvoidanceExp::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &,
                                                         const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calcDiff", &ResidualModelObstacleAvoidanceExp::calcDiff, bp::args("self", "data", "x", "u") )
        .def("createData", &ResidualModelObstacleAvoidanceExp::createData, bp::with_custodian_and_ward_postcall<0, 2>(),
            bp::args("self", "collector") )
        .add_property("pair_id", &ResidualModelObstacleAvoidanceExp::get_pair_id)
        .add_property("frame_id", &ResidualModelObstacleAvoidanceExp::get_frame_id)
        .add_property("type", &ResidualModelObstacleAvoidanceExp::get_type, &ResidualModelObstacleAvoidanceExp::set_type);

    bp::register_ptr_to_python<boost::shared_ptr<ResidualDataObstacleAvoidanceExp> >();

    bp::class_<ResidualDataObstacleAvoidanceExp, bp::bases<ResidualDataAbstract> >(
        "ResidualDataObstacleExp",
        bp::init<ResidualModelObstacleAvoidanceExp*, DataCollectorAbstract*>(
            bp::args("self", "model", "collector")
        )[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 3> >()]
    )
        .add_property("pinocchio",
                      bp::make_getter(&ResidualDataObstacleAvoidanceExp::pinocchio, bp::return_internal_reference<>() ))
        .add_property("geometry",
                      bp::make_getter(&ResidualDataObstacleAvoidanceExp::geometry, bp::return_internal_reference<>() ));
}

}   // namespace python
}   // namespace crocoddyl
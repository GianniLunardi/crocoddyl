// Author: Gianni Lunardi, University of Trento 2022

#include "python/crocoddyl/multibody/multibody.hpp"
#include "crocoddyl/multibody/residuals/obstacle-avoidance.hpp"

namespace crocoddyl {
namespace python {

void exposeResidualObstacleAvoidance() {
    bp::register_ptr_to_python<boost::shared_ptr<ResidualModelObstacleAvoidance> >();

    bp::class_<ResidualModelObstacleAvoidance, bp::bases<ResidualModelAbstract> >(
        "ResidualModelObstacleAvoidance",
        bp::init<boost::shared_ptr<StateMultibody>, std::size_t, boost::shared_ptr<pinocchio::GeometryModel>,
            pinocchio::PairIndex, pinocchio::FrameIndex, pinocchio::ReferenceFrame, double>(
                bp::args("self", "state", "nu", "geom_model", "pair_id", "frame_id", "type", "beta") )
    )
        .def<void (ResidualModelObstacleAvoidance::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                                      const Eigen::Ref<const Eigen::VectorXd>&,
                                                      const Eigen::Ref<const Eigen::VectorXd>&)>(
            "calc", &ResidualModelObstacleAvoidance::calc, bp::args("self", "data", "x", "u") )
        .def<void (ResidualModelObstacleAvoidance::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                                      const Eigen::Ref<const Eigen::VectorXd>&,
                                                      const Eigen::Ref<const Eigen::VectorXd>&)>(
            "calcDiff", &ResidualModelObstacleAvoidance::calcDiff, bp::args("self", "data", "x", "u") )
        .def("createData", &ResidualModelObstacleAvoidance::createData, bp::with_custodian_and_ward_postcall<0, 2>(),
            bp::args("self", "collector") )
        .add_property("pair_id",&ResidualModelObstacleAvoidance::get_pair_id)
        .add_property("frame_id", &ResidualModelObstacleAvoidance::get_frame_id)
        .add_property("type", &ResidualModelObstacleAvoidance::get_type, &ResidualModelObstacleAvoidance::set_type);

    bp::register_ptr_to_python<boost::shared_ptr<ResidualDataObstacleAvoidance> >();

    bp::class_<ResidualDataObstacleAvoidance, bp::bases<ResidualDataAbstract> >(
        "ResidualDataObstacleAvoidance",
        bp::init<ResidualModelObstacleAvoidance*, DataCollectorAbstract*>(
            bp::args("self", "model", "collector")
                )[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 3> >()]
    )
        .add_property("pinocchio",
                      bp::make_getter(&ResidualDataObstacleAvoidance::pinocchio, bp::return_internal_reference<>() ))
        .add_property("geometry",
                       bp::make_getter(&ResidualDataObstacleAvoidance::geometry, bp::return_internal_reference<>() ));
}


}   // namespace crocoddyl
}   // namespace python

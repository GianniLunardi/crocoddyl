// Author: Gianni Lunardi, University of Trento 2022

#include "python/crocoddyl/multibody/multibody.hpp"
#include "crocoddyl/multibody/residuals/fly-high-sqr.hpp"

namespace crocoddyl {
namespace python {

void exposeResidualFlyHighSqr() {
    bp::register_ptr_to_python<boost::shared_ptr<ResidualModelFlyHighSqr> >();

    bp::class_<ResidualModelFlyHighSqr, bp::bases<ResidualModelAbstract> >(
        "ResidualModelFlyHighSqr",
        bp::init<
            boost::shared_ptr<StateMultibody>,
            std::size_t,
            pinocchio::FrameIndex,
            double
        >( bp::args("self", "state", "nu", "frame_id", "beta") )
    )
        .def<void (ResidualModelFlyHighSqr::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calc", &ResidualModelFlyHighSqr::calc, bp::args("self", "data", "x", "u") )
        .def<void (ResidualModelFlyHighSqr::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calcDiff", &ResidualModelFlyHighSqr::calcDiff, bp::args("self", "data", "x", "u") )
        .def("createData", &ResidualModelFlyHighSqr::createData, bp::with_custodian_and_ward_postcall<0, 2>(),
            bp::args("self", "collector") )
        .add_property("frame_id", &ResidualModelFlyHighSqr::get_frame_id)
        .add_property("type", &ResidualModelFlyHighSqr::get_beta);

    bp::register_ptr_to_python<boost::shared_ptr<ResidualDataFlyHighSqr> >();

    bp::class_<ResidualDataFlyHighSqr, bp::bases<ResidualDataAbstract> >(
        "ResidualDataFlyHighSqr",
        bp::init<ResidualModelFlyHighSqr*, DataCollectorAbstract*>(
            bp::args("self", "model", "collector")
        )[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 3> >()]
    )
        .add_property("pinocchio",
            bp::make_getter(&ResidualDataFlyHighSqr::pinocchio, bp::return_internal_reference<>() ));
}

}   // namespace python
}   // namespace crocoddyl

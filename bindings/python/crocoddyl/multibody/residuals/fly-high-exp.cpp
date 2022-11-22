// Author: Gianni Lunardi, University of Trento 2022

#include "python/crocoddyl/multibody/multibody.hpp"
#include "crocoddyl/multibody/residuals/fly-high-exp.hpp"

namespace crocoddyl {
namespace python {

void exposeResidualFlyHighExp() {
    bp::register_ptr_to_python<boost::shared_ptr<ResidualModelFlyHighExp> >();

    bp::class_<ResidualModelFlyHighExp, bp::bases<ResidualModelAbstract> >(
        "ResidualModelFlyHighExp",
        bp::init<
            boost::shared_ptr<StateMultibody>,
            std::size_t,
            pinocchio::FrameIndex,
            double
        >( bp::args("self", "state", "nu", "frame_id", "beta") )
    )
        .def<void (ResidualModelFlyHighExp::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calc", &ResidualModelFlyHighExp::calc, bp::args("self", "data", "x", "u") )
        .def<void (ResidualModelFlyHighExp::*)(const boost::shared_ptr<ResidualDataAbstract> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &,
                                               const Eigen::Ref<const Eigen::VectorXd> &)>(
            "calcDiff", &ResidualModelFlyHighExp::calcDiff, bp::args("self", "data", "x", "u") )
        .def("createData", &ResidualModelFlyHighExp::createData, bp::with_custodian_and_ward_postcall<0, 2>(),
            bp::args("self", "collector") )
        .add_property("frame_id", &ResidualModelFlyHighExp::get_frame_id)
        .add_property("type", &ResidualModelFlyHighExp::get_gamma);

    bp::register_ptr_to_python<boost::shared_ptr<ResidualDataFlyHighExp> >();

    bp::class_<ResidualDataFlyHighExp, bp::bases<ResidualDataAbstract> >(
        "ResidualDataFlyHighExp",
        bp::init<ResidualModelFlyHighExp*, DataCollectorAbstract*>(
            bp::args("self", "model", "collector")
        )[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 3> >()]
    )
        .add_property("pinocchio",
            bp::make_getter(&ResidualDataFlyHighExp::pinocchio, bp::return_internal_reference<>() ));
}

}   // namespace python
}   // namespace crocoddyl

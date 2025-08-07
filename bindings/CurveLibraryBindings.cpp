#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <memory>

#include "utils/Angle.h"
#include "utils/Position.h"
#include "curves/ClothoidCurve.h"
#include "curves/CurveList.h"
#include "curves/HermiteSplineCurve.h"
#include "utils/G2Solve3Arc.h"

namespace py = pybind11;
class PyBaseCurve : public BaseCurve {
public:
    using BaseCurve::BaseCurve;  // Inherit constructors

    // Override virtuals so Python subclasses can provide implementations
    Position getPosition(double value, double h = 0.01) override {
        PYBIND11_OVERRIDE_PURE(
            Position,   // Return type
            BaseCurve,  // Parent class
            getPosition,// Function name in C++
            value, h    // Arguments
        );
    }

    Position getDerivative(double value) override {
        PYBIND11_OVERRIDE_PURE(
            Position,
            BaseCurve,
            getDerivative,
            value
        );
    }

    // Optional override for isBackward (not pure, but virtual)
    bool isBackward() override {
        PYBIND11_OVERRIDE(
            bool,
            BaseCurve,
            isBackward
        );
    }

    Position getLastPosition() override {
        PYBIND11_OVERRIDE(
            Position,
            BaseCurve,
            getLastPosition
        );
    }
    std::vector<double> getParameters() override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<double>,  // Return type
            BaseCurve,            // Parent class
            getParameters         // Function name
        );
    }
};
PYBIND11_MODULE(curve_library, m) {
    m.doc() = "Python bindings for CurveLibrary";

    // ------------------------
    // Angle class
    // ------------------------
    py::class_<Angle>(m, "Angle")
        // Constructors / Factories
        .def_static("from_degrees", &Angle::fromDegrees)
        .def_static("from_radians", &Angle::fromRadians)

        // Conversion methods
        .def("to_degrees", &Angle::toDegrees)
        .def("to_radians", &Angle::toRadians)
        .def("warp_angle", &Angle::warpAngle, py::return_value_policy::reference)

        // Arithmetic operators
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(py::self / double())
        .def(double() * py::self)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__repr__", [](const Angle &a) {
            return "<Angle " + std::to_string(a.toDegrees()) + " deg>";
        });

    // ------------------------
    // Position class
    // ------------------------
    py::class_<Position>(m, "Position")
        .def(py::init<double, double, Angle, double>(),
             py::arg("x")=0.0, py::arg("y")=0.0,
             py::arg("angle")=AngleConstants::ZERO, py::arg("curvature")=0.0)

        // Getters (as read-only Python properties)
        .def_property_readonly("x", &Position::getX)
        .def_property_readonly("y", &Position::getY)
        .def_property_readonly("angle", &Position::getAngle)
        .def_property_readonly("curvature", &Position::getCurvature)

        // Methods
        .def("get_xy", &Position::getXY)
        .def("get_xyrad", &Position::getXYRad)
        .def("get_xyradcurv", &Position::getXYRadCurv)
        .def("get_distance", &Position::getDistance)
        .def("get_vector_angle", &Position::getVectorAngle)

        // Operators
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(py::self / double())
        .def(double() * py::self)

        .def("__repr__", [](const Position &p) {
            return "<Position x=" + std::to_string(p.getX()) +
                   ", y=" + std::to_string(p.getY()) +
                   ", angle=" + std::to_string(p.getAngle().toDegrees()) +
                   ", curv=" + std::to_string(p.getCurvature()) + ">";
        });


    // Expose BaseCurve with trampoline
    py::class_<BaseCurve, PyBaseCurve /* trampoline */, std::shared_ptr<BaseCurve>>(m, "BaseCurve")
        .def(py::init<double, double>(), py::arg("minValue"), py::arg("maxValue"))
        .def("getPosition", &BaseCurve::getPosition)
        .def("getDerivative", &BaseCurve::getDerivative)
        .def("getMinValue", &BaseCurve::getMinValue)
        .def("getMaxValue", &BaseCurve::getMaxValue)
        .def("getLength", py::overload_cast<const double>(&BaseCurve::getLength), py::arg("h")=0.01)
        .def("getLengthBetween", py::overload_cast<double,double,double>(&BaseCurve::getLength),
             py::arg("ti"), py::arg("t_end"), py::arg("h")=0.01)
        .def("getValueForLength", &BaseCurve::getValueForLength,
             py::arg("ti"), py::arg("length"), py::arg("h")=0.01)
        .def("findNearest", py::overload_cast<Position,double>(&BaseCurve::findNearest),
             py::arg("pos"), py::arg("h"))
        .def("findNearestDetailed", py::overload_cast<Position,double,double,int,double,double>(&BaseCurve::findNearest),
             py::arg("pos"), py::arg("guess"), py::arg("h"), py::arg("maxIter"), py::arg("t_min"), py::arg("t_max"))
        .def("isBackward", &BaseCurve::isBackward)
        .def("getLastPosition", &BaseCurve::getLastPosition)
        .def("getParameters", &BaseCurve::getParameters)
        .def("getFullCurve", &BaseCurve::getFullCurve)
        .def("getCurveType", &BaseCurve::getCurveType);
    // ------------------------
    // ClothoidCurve
    // ------------------------
    py::class_<ClothoidCurve, BaseCurve, std::shared_ptr<ClothoidCurve>>(m, "ClothoidCurve")
        .def(py::init<Position,double,double,double>(),
             py::arg("start"), py::arg("start_curv"),
             py::arg("curvature_rate"), py::arg("length"))

        // Static factory methods
        .def_static("get_from_angle",
            &ClothoidCurve::getClothoidCurve,
            py::arg("start"), py::arg("end_angle"), py::arg("initial_curv"), py::arg("length"))

        .def_static("get_from_delta",
            &ClothoidCurve::getClothoidCurveDelta,
            py::arg("start"), py::arg("delta_angle"), py::arg("initial_curv"), py::arg("length"))

        .def_static("get_from_pos_and_theta_kappa", &ClothoidCurve::getClothoidCurveG0)

        .def_static("get_from_positions",
            [](Position start, Position end, double Dmax, double dmin) {
                auto result = G2Solve3Arc();
                if (result.build(start, end, Dmax, dmin) != -1) {
                    return result.getCurveList();
                }
                return std::make_shared<CurveList>();
            }, py::arg("start"), py::arg("end"), py::arg("Dmax")=0.0, py::arg("dmin")=0.0)

        // Methods
        .def("get_position", &ClothoidCurve::getPosition)
        .def("get_derivative", &ClothoidCurve::getDerivative)
        .def("get_length", &ClothoidCurve::getLength)
        .def("get_value_for_length", &ClothoidCurve::getValueForLength)

        // Accessors
        .def("start1", &ClothoidCurve::start1)
        .def("start_curvature", &ClothoidCurve::start_curvature)
        .def("curvature_rate", &ClothoidCurve::curvature_rate);

    // ------------------------
    // G2Solve3Arc
    // ------------------------
    py::class_<G2Solve3Arc>(m, "G2Solve3Arc")
        .def(py::init<>())

        .def("build", py::overload_cast<Position, Position, double, double>(&G2Solve3Arc::build))
        .def("build_detailed",
             py::overload_cast<double,double,double,double,double,double,double,double,double,double>(
                 &G2Solve3Arc::build),
             py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("kappa0"),
             py::arg("x1"), py::arg("y1"), py::arg("theta1"), py::arg("kappa1"),
             py::arg("Dmax")=0, py::arg("dmax")=0)
        .def("build_fixed_length", &G2Solve3Arc::build_fixed_length)

        .def("get_segment0_curve", &G2Solve3Arc::getSegment0Curve)
        .def("get_segment1_curve", &G2Solve3Arc::getSegment1Curve)
        .def("get_segment_middle_curve", &G2Solve3Arc::getSegmentMiddleCurve)

        .def("reverse", &G2Solve3Arc::reverse);
    // ------------------------
    // CurveList
    // ------------------------
    py::class_<CurveList, std::shared_ptr<CurveList>, BaseCurve>(m, "CurveList")
        .def(py::init<>())
        .def("add_curve", &CurveList::addCurve,
             py::arg("curve"),
             "Add a single curve to the list")
        .def("add_curve_list", &CurveList::addCurveList,
             py::arg("list"),
             "Merge another CurveList into this one")
        .def("get_curve_count", &CurveList::getCurveCount,
             "Return the number of curves in the list")
        .def("get_position", &CurveList::getPosition,
             py::arg("s"), py::arg("h")=0.01,
             "Evaluate the curve list at an arc-length position")
        .def("get_derivative", &CurveList::getDerivative,
             py::arg("s"),
             "Get derivative at an arc-length position")
        .def("get_length", &CurveList::getLength,
             py::arg("ti"), py::arg("t_end"), py::arg("h")=0.01,
             "Total length between ti and t_end")
        .def("get_value_for_length", &CurveList::getValueForLength,
             py::arg("ti"), py::arg("length"), py::arg("h")=0.01,
             "Get parameter value for a given arc-length increment")
        .def("__len__", &CurveList::getCurveCount)
        .def("__repr__", [](const CurveList &cl) {
            return "<CurveList with " + std::to_string(cl.getCurveCount()) + " curves>";
        });

    py::class_<HermiteSplineCurve<5>, std::shared_ptr<HermiteSplineCurve<5>>, BaseCurve>(m, "HermiteSplineCurve")
        .def(py::init<std::array<Position, 6>>())
        .def_static("build", &HermiteSplineCurve<5>::getSplineFromStartEnd, py::arg("start"), py::arg("end"), py::arg("L"))
        .def("eval", &HermiteSplineCurve<5>::eval<3>, py::arg("s"), "Get the spline evaluation");
    py::enum_<CurveFactory::CurveTypes>(m, "CurveType")
        .value("BASE", CurveFactory::CurveTypes::BASE)
        .value("ARC", CurveFactory::CurveTypes::ARC)
        .value("BEZIER", CurveFactory::CurveTypes::BEZIER)
        .value("CLOTHOID", CurveFactory::CurveTypes::CLOTHOID)
        .value("LIST", CurveFactory::CurveTypes::LIST)
        .value("SPLINE", CurveFactory::CurveTypes::SPLINE)
        .export_values();
    m.def("get_base_curve", &CurveFactory::getBaseCurve,
      py::arg("parameters"),
      "Construct a curve from serialized parameters");

}

/**
 * @file sbopy.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "infra/infra_context.hpp"
#include "infra/log.hpp"
#include "infra/status.hpp"
#include "math/precision.hpp"
#include "math/special_functions.hpp"

namespace py = pybind11;
using namespace sbody;

constexpr LogLevel SBOPY_LOG_LEVEL = LogLevel::OFF;
constexpr FILE *SBOPY_STD_OUT = nullptr;
constexpr FILE *SBOPY_STD_ERR = nullptr;

void bindSStarModule(py::module_ &m);
void bindHotspotModule(py::module_ &m);
void bindShadowModule(py::module_ &m);

std::vector<real_t> convertFromNumpyDoubleArray(py::array_t<double> array) {
    std::vector<real_t> result(array.size());
    for (ssize_t i = 0; i < array.size(); ++i)
        result[i] = array.at(i);
    return result;
}

std::vector<uint8_t> convertFromNumpyBoolArray(py::array_t<bool> array) {
    std::vector<uint8_t> result(array.size());
    for (ssize_t i = 0; i < array.size(); ++i)
        result[i] = array.at(i) ? 1 : 0;
    return result;
}

class EmceeTestClass {
  protected:
    double mean_x_;
    double mean_y_;
    double std_x_;
    double std_y_;

  public:
    EmceeTestClass(py::array_t<double> mean, py::array_t<double> std)
        : mean_x_(mean.at(0)), mean_y_(mean.at(1)), std_x_(std.at(0)), std_y_(std.at(1)) {};
    double LogProb(py::array_t<double> theta) const {
        return -0.5 * (power2((theta.at(0) - mean_x_) / std_x_) + power2((theta.at(1) - mean_y_) / std_y_));
    }
};

class NumpyTestClass {
  protected:
    std::vector<real_t> data_;

  public:
    NumpyTestClass(py::array_t<double> data, std::array<std::vector<double>, 2> vec_data) {
        data_ = convertFromNumpyDoubleArray(data);
        for (auto &sub_vec : vec_data) {
            for (auto i : sub_vec)
                std::cout << i << " ";
            std::cout << std::endl;
        }
    };
    void print() {
        for (auto i : data_)
            std::cout << i << " ";
        std::cout << std::endl;
    }
};

PYBIND11_MODULE(sbopy, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: cmake_example
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

    m.attr("__version__") = SBODY_VERSION;

    py::enum_<StatusCode>(m, "StatusCode")
        .value("SUCCESS", StatusCode::SUCCESS)
        .value("CONTINUE", StatusCode::CONTINUE)
        .value("FAILURE", StatusCode::FAILURE)
        .value("INVALID_ARGUMENT", StatusCode::INVALID_ARGUMENT)
        .value("DOMAIN_ERROR", StatusCode::DOMAIN_ERROR)
        .value("DIVISION_BY_ZERO", StatusCode::DIVISION_BY_ZERO)
        .value("MAX_ITERATIONS_EXCEEDED", StatusCode::MAX_ITERATIONS_EXCEEDED)
        .value("SINGULAR_MATRIX", StatusCode::SINGULAR_MATRIX)
        .value("NUMERIC_ERROR", StatusCode::NUMERIC_ERROR)
        .value("TRUST_REGION_TOO_SMALL", StatusCode::TRUST_REGION_TOO_SMALL)
        .value("SCALING_REQUIRED", StatusCode::SCALING_REQUIRED)
        .value("STEP_SIZE_TOO_SMALL", StatusCode::STEP_SIZE_TOO_SMALL)
        .value("NO_INTERSECTION_WITH_TARGET_SURFACE", StatusCode::NO_INTERSECTION_WITH_TARGET_SURFACE)
        .value("FALL_INTO_BLACK_HOLE", StatusCode::FALL_INTO_BLACK_HOLE)
        .value("INPUT_OUTPUT_ERROR", StatusCode::INPUT_OUTPUT_ERROR)
        .value("NOT_IMPLEMENTED", StatusCode::NOT_IMPLEMENTED)
        .value("COMPUTATION_TIMEOUT", StatusCode::COMPUTATION_TIMEOUT);

    py::enum_<LogLevel>(m, "LogLevel")
        .value("TRACE", LogLevel::TRACE)
        .value("DEBUG", LogLevel::DEBUG)
        .value("INFO", LogLevel::INFO)
        .value("WARNING", LogLevel::WARNING)
        .value("ERROR", LogLevel::ERROR)
        .value("OFF", LogLevel::OFF);

    py::class_<InfraContext>(m, "InfraContext").def(py::init<LogLevel>(), py::arg("log_level") = SBOPY_LOG_LEVEL);

    py::class_<EmceeTestClass>(m, "EmceeTestClass")
        .def(py::init<py::array_t<double>, py::array_t<double>>())
        .def("log_prob", &EmceeTestClass::LogProb);

    py::class_<NumpyTestClass>(m, "NumpyTestClass")
        .def(py::init<py::array_t<double>, std::array<std::vector<double>, 2>>())
        .def("print", &NumpyTestClass::print);

    bindSStarModule(m);
    bindHotspotModule(m);
    bindShadowModule(m);

    m.def("elliptic_integral", &ellipticIntegral<double>, R"pbdoc(
		Elliptic Integral
	)pbdoc");

    m.def("elliptic_integral_2_complex", &ellipticIntegral2Complex<double>, R"pbdoc(
		Elliptic Integral with 2 complex roots
	)pbdoc");

    m.def("elliptic_integral_4_complex", &ellipticIntegral4Complex<double>, R"pbdoc(
		Elliptic Integral with 4 complex roots
	)pbdoc");
}

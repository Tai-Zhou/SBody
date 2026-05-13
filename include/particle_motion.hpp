/**
 * @file particle_motion.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-16
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <chrono>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "config/config.hpp"
#include "infra/infra_context.hpp"
#include "infra/log.hpp"
#include "metric/metrics.hpp"
#include "object/object.hpp"
#include "object/particle.hpp"
#include "units.hpp"
#include "view/camera.hpp"
#include "view/view.hpp"

namespace sbody {
template <typename T>
int particleMotion(RuntimeConfig<T> &config) {
    Unit unit(config.mass_);
    auto wall_time_start = std::chrono::steady_clock::now();
    int t_step_number = config.t_step_number_;
    int metric = 1;
    size_t Hamiltonian = 0;
    size_t PN = 1;
    size_t ray = 1;
    int rest_mass = 1;
    bool display_progress = true;
    std::string store_format = "NumPy";
    T inc = M_PI * 0. / 180., eps = 0.;
    // T a = 8.3, e = 0., inclination = M_PI * 162. / 180., periapsis = M_PI
    // * 198.9 / 180. + M_PI_2, ascending_node = M_PI * 25.1 / 180., true_anomaly
    // = M_PI_2; //
    // phi=[2.73633242 3.92974873 3.32166381 3.2093593 3.67372211 5.18824159
    // | 3.19861806 2.63708292 3.05259405]
    T a = 4.999960135538237, e = 0.884649, inclination = M_PI * 134.567 / 180., periapsis = M_PI * 66.263 / 180.,
      ascending_node = M_PI * 228.171 / 180., true_anomaly = M_PI;
    T input_parameter[10];
    std::unique_ptr<View<T>> viewPtr;
    std::unique_ptr<std::thread> shadowPtr;
    std::unique_ptr<Camera<T>> camera_ptr;
    std::unique_ptr<std::thread> lens_ptr;
    T tFinal = config.t_final_ * unit.s;
    tFinal = 18 * unit.yr;
    T t = 8.3564958220072 * unit.yr, tStep = 0., tRec = tFinal / t_step_number;
    if (metric == 0)
        ray = 0;
    std::unique_ptr<Metric<T>> main_metric = createMetric(config.metric_, config.metric_params_);
    // metric::setMetric(metric, PN, mass, spin, charge, NUT);
    std::string strFormat =
        fmt::format(" ({:f},{:f})[{:f},{:f}]", config.metric_params_.spin_, config.metric_params_.NUT_, inc, eps);
    InfraContext infra_context;
    if (display_progress) {
        infra_context.progress_bar_ = new ProgressBar(metric, std::string("?") + strFormat);
    }
    if (ray & 5) {
        // viewPtr = make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN),
        // 8180. * unit.pc, inc, fmt::format("view
        // ({:.1f},{:.1f},{:.1f})[{:f},{:f}]", spin, charge, NUT, inc, eps));
        // viewPtr = make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN),
        // 8180. * unit.pc, inc, fmt::format("view
        // conical({:.2f},{:.6f},{:.6f},{:.6f},{:.6f})[{:f},{:f}]", 10.6, M_PI *
        // 0.75, M_PI * 7. / 9., 0.01, 0.45 / 10.6, inc, eps)); viewPtr =
        // make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN), 8180. *
        // unit.pc, inc, fmt::format("view
        // conical({:f},{:f},{:f},{:f},{:f})[{:f},{:f}]", a, inclination, periapsis,
        // e, ascending_node, inc, eps));
        viewPtr = std::make_unique<View<T>>(*main_metric, 8180. * unit.pc, inc, eps);
        if (ray & 4)
            shadowPtr = std::make_unique<std::thread>(&View<T>::Shadow, viewPtr.get(), "shadow", std::ref(infra_context));
    }
    if (ray & 10) {
        camera_ptr = std::make_unique<Camera<T>>(*main_metric, 1000, 5e-2, 1.e3, inc, eps);
        if (ray & 8)
            lens_ptr = std::make_unique<std::thread>(&Camera<T>::Lens, camera_ptr.get(), std::ref(infra_context));
    }
    // Integrator integ(metric::function, metric::jacobian, metric != 0);
    TimeSystem star_time = TimeSystem::T;
    Particle star_0(*main_metric, star_time, DynamicalSystem::LAGRANGIAN, config.absolute_accuracy_, config.relative_accuracy_,
                    false);
    star_0.initializeSchwarzschildKeplerianApocenter(a * unit.mpc, e, inclination, periapsis, ascending_node, inc, eps);
    // star_0.InitializeGeodesic(a, inclination, periapsis, ascending_node,
    // -0.16707659553531468, 0.3822615764261866, inc, eps);
    // star_0.InitializeSchwarzschildKeplerianApocenterHarmonic(a * unit.mpc, e,
    // inclination, periapsis, ascending_node, inc, eps);
    // star_0.InitializeConical(10.6, M_PI * 0.75, M_PI * 7. / 9., 0.01, 0.45
    // / 10.6); Integrator &&integrator = main_metric->GetIntegrator(metric != 0);
    // NumPy rec(main_metric->Name() + strFormat, {12});
    // NumPy rec(fmt::format("Hotspot a={:.1f} e={:.2f} i={:.6f} o={:.2f}", a, e,
    // inclination, periapsis), {12});
    File rec(FileFormat::NUMPY, fmt::format("mcmc test", a, inclination, periapsis, e, ascending_node), {17});
    boost::numeric::ublas::bounded_vector<T, 8> position;
    boost::numeric::ublas::bounded_vector<T, 5> view_info;
    int status = 0;
    T h = 1., stepPercent = 100. / t_step_number;
    h = -1;
    auto object_pointer_collector = std::vector<Object<T> *>{&star_0};
    if (status = star_0.integratorApply(t, 0., h); status != 0)
        infra_context.logger_.warning("star_0.IntegratorApply() = {}", status);
    h = 1;
    for (int i = 1; i <= t_step_number; ++i) {
        tStep += tRec;
        status = star_0.integratorApply(t, tStep, h);
        if (status > 0)
            infra_context.logger_.error("main status = {}", status);
        star_0.position(position);
        if (Hamiltonian)
            main_metric->hamiltonianToLagrangian(position);
        if (ray & 1)
            viewPtr->Trace(position, star_time, view_info, config.absolute_accuracy_, config.relative_accuracy_, true, false,
                           infra_context);
        if (ray & 2)
            camera_ptr->Trace(object_pointer_collector, infra_context);
        sphericalToCartesian(position);
        rec.save(position);
        rec.save(t / unit.s);
        rec.save(star_0.energy());
        rec.save(star_0.angularMomentum());
        rec.save(star_0.carterConstant());
        view_info[3] /= unit.s;
        rec.save(view_info);
        if (display_progress)
            infra_context.progress_bar_->setProgress(0, i * stepPercent);
        if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - wall_time_start).count() >=
            config.t_calculation_)
            break;
    }
    if (display_progress) {
        infra_context.progress_bar_->setComplete(0, std::string("!") + strFormat);
    }
    if (ray & 2)
        camera_ptr->Save("camera");
    if (ray & 4)
        shadowPtr->join();
    if (ray & 8)
        lens_ptr->join();
    return StatusCode::SUCCESS;
}
} // namespace sbody

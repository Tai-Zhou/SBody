/**
 * @file benchmark.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "benchmark.hpp"

#include <chrono>
#include <memory>

#include "infra/infra_context.hpp"
#include "infra/status.hpp"
#include "math/precision.hpp"
#include "metric/kerr.hpp"
#include "object/particle.hpp"
#include "units.hpp"
#include "view/view.hpp"

using namespace std;

namespace sbody {
int serialBenchmark() {
    auto t_start = chrono::steady_clock::now();
    Unit unit(4.261e6);
    uint T_STEP_NUMBER = 10000;
    real_t t = 0., t1 = 0., t_step = 3600. * unit.s / T_STEP_NUMBER;
    auto main_metric = make_unique<Kerr<real_t>>(0.3);
    // unique_ptr<Metric> main_metric = make_unique<Schwarzschild>();
    auto view_ptr =
        make_unique<View<real_t>>(*main_metric, 8180. * unit.pc, boost::math::constants::half_pi<real_t>() - 0.01, 0.);
    Particle<real_t> star_0(*main_metric, TimeSystem::T, DynamicalSystem::LAGRANGIAN, 1e-12, 1e-12, false);
    star_0.initializeKeplerian(4., 0.1, 0.2, 0., 0., 0., M_PI_4);
    real_t temp[4];
    boost::numeric::ublas::bounded_vector<real_t, 8> position;
    boost::numeric::ublas::bounded_vector<real_t, 5> view_info;
    int status = 0;
    real_t h = 1., stepPercent = 100. / T_STEP_NUMBER;
    indicators::show_console_cursor(false);
    ProgressBar bars(2, "Benchmarking...");
    InfraContext infra_context;
    for (uint i = 1; i <= T_STEP_NUMBER; ++i) {
        t1 += t_step;
        status = star_0.integratorApply(t, t1, h);
        if (status > 0) {
            infra_context.logger_.error("Benchmark error, main status = {}", status);
            return StatusCode::FAILURE;
        }
        star_0.position(temp);
        temp[0] = t / unit.s;
        temp[1] = star_0.energy();
        temp[2] = star_0.angularMomentum();
        temp[3] = star_0.carterConstant();
        view_ptr->Trace(position, TimeSystem::T, view_info, 1e-12, 1e-12, false, true, infra_context);
        sphericalToCartesian(position);
        bars[0].set_progress(i * stepPercent);
        bars[0].set_option(indicators::option::PrefixText{fmt::format(
            "Anticipated score is: {}? ",
            i * 200000000000ul /
                (T_STEP_NUMBER * chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - t_start).count()))});
    }
    bars.setComplete(
        0, fmt::format("Final score is: {}! ",
                       200000000000ul /
                           chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - t_start).count()));
    return StatusCode::SUCCESS;
}
int parallelBenchmark(int thread_num) {
    cout << thread_num << endl;
    return StatusCode::NOT_IMPLEMENTED;
}
} // namespace sbody

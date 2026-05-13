/**
 * @file main.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <csignal>

#include <CLI/CLI.hpp>
#include <fmt/core.h>

#include "benchmark.hpp"
#include "config/config.hpp"
#include "math/precision.hpp"
#include "particle_motion.hpp"

using namespace sbody;

void InterruptHandler(int signum) {
    fmt::print("\r\033[KSBody exiting...\n");
    indicators::show_console_cursor(true);
    exit(signum);
}

void registerParticleMotion(CLI::App &app, RuntimeConfig<real_t> &config) {
    app.add_option("-m,--mass", config.mass_, "mass of the source (M_sun)")->default_str(fmt::format("{}", config.mass_));
    app.add_option("-n,--metric", config.metric_, "Post-Newtonian/Schwarzschild/Kerr/KerrTaubNUT")
        ->default_str(fmt::format("{}", config.metric_name_));
    app.add_option("-s,--spin", config.metric_params_.spin_, "spin of the source")
        ->default_str(fmt::format("{}", config.metric_params_.spin_));
    app.add_option("-l,--nut", config.metric_params_.NUT_, "NUT charge of the source")
        ->default_str(fmt::format("{}", config.metric_params_.NUT_));
    app.add_option("-R,--sa", config.kepler_[0], "semimajor axis of the star")
        ->default_str(fmt::format("{}", config.kepler_[0]));
    app.add_option("-E,--ec", config.kepler_[1], "eccentricity of the star")->default_str(fmt::format("{}", config.kepler_[1]));
    app.add_option("-I,--in", config.kepler_[2], "inclination of the star (deg)")
        ->default_str(fmt::format("{}", config.kepler_[2]));
    app.add_option("-o,--pe", config.kepler_[3], "position angle of periapsis of the star (deg)")
        ->default_str(fmt::format("{}", config.kepler_[3]));
    app.add_option("-O,--an", config.kepler_[4], "ascending node of the star (deg)")
        ->default_str(fmt::format("{}", config.kepler_[4]));
    app.add_option("-t,--time", config.t_final_, "time limit (s)")->default_str(fmt::format("{}", config.t_final_));
    app.add_option("-k,--rec", config.t_step_number_, "record number")->default_str(fmt::format("{}", config.t_step_number_));
    app.add_option("-c,--tcal", config.t_calculation_, "real time limit of calculation (s)")
        ->default_str(fmt::format("{}", config.t_calculation_));
    app.add_option("-a,--abs", config.absolute_accuracy_, "absolute accuracy")
        ->default_str(fmt::format("{}", config.absolute_accuracy_));
    app.add_option("-r,--rel", config.relative_accuracy_, "relative accuracy")
        ->default_str(fmt::format("{}", config.relative_accuracy_));
    app.add_option("-i,--inc", config.inclination_, "inclination of the BH")
        ->default_str(fmt::format("{}", config.inclination_));
    app.add_option("-e,--eps", config.position_angle_, "epsilon of the BH")
        ->default_str(fmt::format("{}", config.position_angle_));
    app.add_option("-f,--format", config.output_dir_, "storage folder")->default_str(fmt::format("{}", config.output_dir_));
    app.add_option("-f,--format", config.output_format_, "storage format")
        ->default_str(fmt::format("{}", config.output_format_));
    app.add_flag("-b,--pbar", config.progress_bar_, "display progress bar");
}

// void registerShadow(CLI::App &app, RuntimeConfig &config) { // TODO:
// }

void registerBenchmarks(CLI::App &app) {
    auto cli_benchmark = app.add_subcommand("benchmark", "run benchmark");
    cli_benchmark->require_subcommand(1);
    auto cli_benchmark_serial = cli_benchmark->add_subcommand("serial");
    auto cli_benchmark_parallel = cli_benchmark->add_subcommand("parallel");
    size_t thread_num = 0;
    cli_benchmark_parallel->add_option("-t,--threads", thread_num)->default_str("auto");
    cli_benchmark_serial->callback(serialBenchmark);
    cli_benchmark_parallel->callback([thread_num]() { parallelBenchmark(thread_num); });
}

int main(int argc, char *argv[]) {
    CLI::App cli_sbody{"A relativistic ray tracing program."};
    cli_sbody.set_version_flag("-V,--version", fmt::format("SBody {} ({})\n\nC++: {}\nCompiler: {} {}\nCMake: {}\nBoost: {}",
                                                           SBODY_VERSION, PRECISION_NAME, SBODY_CPP_STANDARD, SBODY_COMPILER_ID,
                                                           SBODY_COMPILER_VERSION, SBODY_CMAKE_VERSION, SBODY_BOOST_VERSION));
    RuntimeConfig<real_t> runtime_cfg;
    registerParticleMotion(cli_sbody, runtime_cfg);
    registerBenchmarks(cli_sbody);
    cli_sbody.footer("SUPPORT:\n  https://github.com/tai-zhou/SBody");
    cli_sbody.get_formatter()->enable_footer_formatting(false);
    // CLI11_PARSE(cli_sbody, argc, argv);
    try {
        cli_sbody.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return cli_sbody.exit(e);
    }
    if (argc == 1) {
        std::cout << cli_sbody.help();
        return 0;
    }
    if (runtime_cfg.validate() == false) {
        fmt::print(stderr, "Invalid configuration!\n");
        return 1;
    }
    signal(SIGINT, InterruptHandler);
    return particleMotion(runtime_cfg);
}

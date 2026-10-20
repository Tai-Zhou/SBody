/**
 * @file camera.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include "infra/log.hpp"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "infra/infra_context.hpp"
#include "infra/status.hpp"
#include "object/object.hpp"
#include "view/view.hpp"

namespace sbody {
template <typename T>
class Camera : public View<T> {
  protected:
    const std::size_t pixel_;
    const T half_angle_;
    std::vector<boost::numeric::ublas::bounded_vector<T, 8>> initial_positions_;
    std::vector<T> initial_times_;
    std::vector<std::vector<T>> screen_;

  public:
    /**
     * @brief Construct a new camera object
     *
     * @param pixel
     * @param half_angle
     * @param r
     * @param theta
     * @param file_name
     */
    Camera(const Metric<T> &metric, std::size_t pixel, T half_angle, T r, T theta, T iota)
        : View<T>(metric, r, theta, iota), pixel_(pixel), half_angle_(half_angle) {
        screen_ = std::vector<std::vector<T>>(pixel, std::vector<T>(pixel));
        initial_positions_ = std::vector<boost::numeric::ublas::bounded_vector<T, 8>>(pixel * pixel);
        const T pixel_size = 2. * half_angle * r / pixel, t1 = r + 100.;
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_fehlberg78<boost::numeric::ublas::bounded_vector<T, 8>>>
            integration_stepper;
        for (std::size_t i = 0; i < pixel; ++i)
            for (std::size_t j = 0; j < pixel; ++j)
                this->InitializePhoton(initial_positions_[i * pixel + j], initial_times_[i * pixel + j],
                                       pixel_size * (i - 0.5 * pixel + 0.5), pixel_size * (j - 0.5 * pixel + 0.5));
#pragma omp parallel for
        for (int p = pixel * pixel - 1; p >= 0; --p) {
            this->metric_.normalizeNullGeodesic(initial_positions_[p], 1.);
            this->metric_.lagrangianToHamiltonian(initial_positions_[p]);
            boost::numeric::odeint::integrate_adaptive(integration_stepper, this->integration_system_, initial_positions_[p],
                                                       initial_times_[p], t1, -1.0);
        }
    }

    /**
     * @brief
     *
     * @return int
     */
    int Trace(std::vector<Object<T> *> &object_list, const InfraContext &infra_context) {
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_fehlberg78<boost::numeric::ublas::bounded_vector<T, 8>>>
            integration_stepper;
        const T t1 = -1000.;
#pragma omp parallel for
        for (int p = pixel_ * pixel_ - 1; p >= 0; --p) {
            int i = p / pixel_;
            int j = p - i * pixel_;
            boost::numeric::ublas::bounded_vector<T, 8> photon, last;
            T t = initial_times_[p], h = 1.0;
            int status = 0;
            std::copy(initial_positions_[p].begin(), initial_positions_[p].end(), photon.begin());
            while (status <= 0 && t > t1) {
                std::copy(photon.begin(), photon.end(), last.begin());
                integration_stepper.try_step(this->integration_system_, photon, t, h);
                for (auto object_pointer : object_list)
                    if (object_pointer->hit(photon, last))
                        screen_[i][j] = object_pointer->dopplerFactor(photon, TimeSystem::T); // FIXME: if multi objects
                if (screen_[i][j] > boost::math::tools::root_epsilon<T>())
                    break;
            }
            if (status > 0)
                infra_context.logger_.warning("Camera::Trace() status = {}", status);
        }
        return StatusCode::SUCCESS;
    }

    /**
     * @brief
     *
     * @return int
     */
    int Lens(InfraContext &infra_context) {
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_fehlberg78<boost::numeric::ublas::bounded_vector<T, 8>>>
            integration_stepper;
        const T t1 = -1000. * this->r_, pixelPerAngle = 0.5 * pixel_ / half_angle_;
        File rec(FileFormat::NUMPY, "lens", {2});
        int progress_bar_index = -1;
        if (infra_context.progress_bar_ != nullptr) {
            indicators::BlockProgressBar bar{
                indicators::option::ForegroundColor{indicators::Color(4)},
                indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                indicators::option::MaxProgress(pixel_ * pixel_),
                indicators::option::PrefixText("? Lens"),
                indicators::option::ShowElapsedTime{true},
                indicators::option::ShowRemainingTime{true}};
            progress_bar_index = infra_context.progress_bar_->push_back(bar);
        }
        for (std::size_t i = 0; i < pixel_; ++i)
            for (std::size_t j = 0; j < pixel_; ++j) {
                std::size_t idx = i * pixel_ + j;
                boost::numeric::ublas::bounded_vector<T, 8> photon;
                T t = initial_times_[idx], h = -1;
                int status = 0;
                std::copy(initial_positions_[idx].begin(), initial_positions_[idx].end(), photon.begin());
                boost::numeric::odeint::integrate_adaptive(integration_stepper, this->integration_system_, photon, t, t1, -1.0);
                if (status > 0)
                    infra_context.logger_.warning("Camera::Lens() status = {}", status);
                if (photon[1] <= 3.)
                    rec.save({NAN, NAN});
                else {
                    this->metric_.hamiltonianToLagrangian(photon);
                    if (photon[2] < 0)
                        photon[2] += boost::math::constants::pi<T>();
                    sphericalToCartesian(photon);
                    rec.save({photon[6] * pixelPerAngle,
                              (photon[7] * this->sin_theta_ - photon[5] * this->cos_theta_) * pixelPerAngle});
                }
                if (progress_bar_index >= 0)
                    infra_context.progress_bar_->tick(progress_bar_index);
            }
        if (progress_bar_index >= 0)
            infra_context.progress_bar_->setComplete(progress_bar_index, "! lens");
        return StatusCode::SUCCESS;
    }

    /**
     * @brief
     *
     * @return int
     */
    int Save(std::string file_name) {
        File record(FileFormat::NUMPY, file_name, {2});
        for (const std::vector<T> &line : screen_)
            record.save(line);
        return StatusCode::SUCCESS;
    }
};
} // namespace sbody

/**
 * @file progress.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "infra/progress.hpp"

#include <string>

namespace sbody {
ProgressBar::ProgressBar(int color, std::string prefix_text) {
    indicators::show_console_cursor(false);
    indicators::BlockProgressBar bar{
        indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
        indicators::option::ForegroundColor{indicators::Color(color)},
        indicators::option::PrefixText{prefix_text},
        indicators::option::ShowElapsedTime{true},
        indicators::option::ShowRemainingTime{true},
    };
    this->push_back(bar);
}
ProgressBar ::~ProgressBar() { indicators::show_console_cursor(true); }
void ProgressBar::tick(int index) { this->operator[](index).tick(); }
void ProgressBar::setProgress(int index, double progress) { this->operator[](index).set_progress(progress); }
/**
 * @brief Set a specific progress bar as completed.
 *
 * @param index index of the progress bar in the container
 * @param prefix prefix shown in the completed progress bar
 */
void ProgressBar::setComplete(int index, std::string prefix) {
    this->operator[](index).set_progress(100.);
    this->operator[](index).set_option(indicators::option::PrefixText{prefix});
    this->operator[](index).mark_as_completed();
}

} // namespace sbody

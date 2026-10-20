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

#pragma once

#include <cstdint>
#include <string>

#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/dynamic_progress.hpp>

namespace sbody {
/**
 * @brief
 *
 */
class ProgressBar : public indicators::DynamicProgress<indicators::BlockProgressBar> {
  public:
    ProgressBar(int color, std::string prefix_text);
    ~ProgressBar();
    void tick(int index);
    /**
     * @brief Set the Progress object
     *
     * @param index
     * @param progress
     */
    void setProgress(int index, double progress);
    /**
     * @brief Set a specific progress bar as completed.
     *
     * @param index index of the progress bar in the container
     * @param prefix prefix shown in the completed progress bar
     */
    void setComplete(int index, std::string prefix);
};

} // namespace sbody

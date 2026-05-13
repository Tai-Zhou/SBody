/**
 * @file infra_context.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-04-08
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <filesystem>

#include "file/file.hpp"
#include "infra/log.hpp"

namespace sbody {
class Logger;
class ProgressBar;

class InfraContext {
  public:
    Logger logger_;
    std::filesystem::path output_dir_;
    FileFormat output_format_;
    ProgressBar *progress_bar_;
    InfraContext(LogLevel log_level = LogLevel::DEBUG, FILE *std_out = stdout, FILE *std_err = stderr,
                 std::filesystem::path output_dir = "data", FileFormat output_format = FileFormat::NUMPY,
                 ProgressBar *progress_bar = nullptr)
        : logger_(log_level, std_out, std_err), output_dir_(output_dir), output_format_(output_format),
          progress_bar_(progress_bar) {}
};
} // namespace sbody

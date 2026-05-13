/**
 * @file log.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <cstdio>
#include <string>

#include <fmt/core.h>
#include <fmt/format.h>

// #include "Utility.hpp"

namespace sbody {
enum class LogLevel {
    TRACE = -2,
    DEBUG = -1,
    INFO = 0,
    WARNING = 1,
    ERROR = 2,
    OFF = 3,
};

class Logger {
  protected:
    LogLevel log_level_;
    FILE *std_out_;
    FILE *std_err_;

  public:
    Logger(LogLevel log_level, FILE *std_out, FILE *std_err) : log_level_(log_level), std_out_(std_out), std_err_(std_err) {}
    template <typename... T>
    void print(std::string format, T &&...args) const {
        fmt::print(std_out_, format, std::forward<T>(args)...);
    }
    template <typename... T>
    void println(std::string format, T &&...args) const {
        fmt::println(std_out_, format, std::forward<T>(args)...);
    }
    template <typename... T>
    void trace(std::string format, T &&...args) const {
        if (log_level_ > LogLevel::TRACE || std_out_ == nullptr)
            return;
        fmt::println(std_err_, "\033[100m[TRACE]\033[0m " + format, std::forward<T>(args)...);
    }
    template <typename... T>
    void debug(std::string format, T &&...args) const {
        if (log_level_ > LogLevel::DEBUG || std_out_ == nullptr)
            return;
        fmt::println(std_err_, "\033[102m[DEBUG]\033[0m " + format, std::forward<T>(args)...);
    }
    template <typename... T>
    void info(std::string format, T &&...args) const {
        if (log_level_ > LogLevel::INFO || std_out_ == nullptr)
            return;
        fmt::println(std_err_, "\033[104m[INFO]\033[0m " + format, std::forward<T>(args)...);
    }
    template <typename... T>
    void warning(std::string format, T &&...args) const {
        if (log_level_ > LogLevel::WARNING || std_err_ == nullptr)
            return;
        fmt::println(std_err_, "\033[103m[WARN]\033[0m " + format, std::forward<T>(args)...);
    }
    template <typename... T>
    void error(std::string format, T &&...args) const {
        if (log_level_ > LogLevel::ERROR || std_err_ == nullptr)
            return;
        fmt::println(std_err_, "\033[101m[ERROR]\033[0m " + format, std::forward<T>(args)...);
    }
};
} // namespace sbody

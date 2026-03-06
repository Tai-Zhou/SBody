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

#ifndef SBODY_LOG_HPP
#define SBODY_LOG_HPP

#include <string>

#include <fmt/core.h>
#include <fmt/format.h>

// #include "Utility.hpp"

namespace sbody {
	template <typename... T>
	void PrintlnBold(std::string format, T &&...args) {
		fmt::println(format, "\033[1m", "\033[0m", std::forward<T>(args)...);
	}
	template <typename... T>
	void PrintlnWarning(std::string format, T &&...args) {
		fmt::println(stdout, "\033[103m[WRN]\033[0m " + format, std::forward<T>(args)...);
	}
	template <typename... T>
	void PrintlnError(std::string format, T &&...args) {
		fmt::println(stdout, "\033[101m[ERR]\033[0m " + format, std::forward<T>(args)...);
	}
} // namespace sbody

#endif

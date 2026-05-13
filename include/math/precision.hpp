/**
 * @file precision.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-16
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

namespace sbody {
#if SBODY_PRECISION == 0
	using real_t = double;
	constexpr const char *PRECISION_NAME = "double";
#elif SBODY_PRECISION == 1
	using real_t = long double;
	constexpr const char *PRECISION_NAME = "long double";
#endif
} // namespace sbody

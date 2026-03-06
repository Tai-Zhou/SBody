/**
 * @file error.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */
#ifndef SBODY_ERROR_HPP
#define SBODY_ERROR_HPP

#include <string>

namespace sbody {
	enum StatusCode {
		// General status
		CONTINUE = -1,
		SUCCESS = 0,
		FAILURE = 1,

		// Numerical
		INVALID_ARGUMENT,
		DOMAIN_ERROR,
		DIVISION_BY_ZERO,
		MAX_ITERATIONS_EXCEEDED,
		SINGULAR_MATRIX,
		NUMERIC_ERROR,
		TRUST_REGION_TOO_SMALL,
		SCALING_REQUIRED,
		STEP_SIZE_TOO_SMALL,

		// Physical
		NO_INTERSECTION_WITH_TARGET_SURFACE,
		FALL_INTO_BLACK_HOLE,

		// Other
		INPUT_OUTPUT_ERROR,
		NOT_IMPLEMENTED,
		COMPUTATION_TIMEOUT,
	};

	std::string to_string(StatusCode code);
} // namespace sbody
#endif

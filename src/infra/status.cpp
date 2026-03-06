/**
 * @file status.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-05
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "infra/status.hpp"

namespace sbody {
	std::string to_string(StatusCode code) {
		switch (code) {
		case CONTINUE:
			return "CONTINUE";
		case SUCCESS:
			return "SUCCESS";
		case FAILURE:
			return "FAILURE";
		case INVALID_ARGUMENT:
			return "INVALID_ARGUMENT";
		case DOMAIN_ERROR:
			return "DOMAIN_ERROR";
		case DIVISION_BY_ZERO:
			return "DIVISION_BY_ZERO";
		case MAX_ITERATIONS_EXCEEDED:
			return "MAX_ITERATIONS_EXCEEDED";
		case SINGULAR_MATRIX:
			return "SINGULAR_MATRIX";
		case NUMERIC_ERROR:
			return "NUMERIC_ERROR";
		case TRUST_REGION_TOO_SMALL:
			return "TRUST_REGION_TOO_SMALL";
		case SCALING_REQUIRED:
			return "SCALING_REQUIRED";
		case STEP_SIZE_TOO_SMALL:
			return "STEP_SIZE_TOO_SMALL";
		case NO_INTERSECTION_WITH_TARGET_SURFACE:
			return "NO_INTERSECTION_WITH_TARGET_SURFACE";
		case FALL_INTO_BLACK_HOLE:
			return "FALL_INTO_BLACK_HOLE";
		case INPUT_OUTPUT_ERROR:
			return "INPUT_OUTPUT_ERROR";
		case NOT_IMPLEMENTED:
			return "NOT_IMPLEMENTED";
		case COMPUTATION_TIMEOUT:
			return "COMPUTATION_TIMEOUT";
		default:
			return "UNKNOWN";
		}
	}
} // namespace sbody

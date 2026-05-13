/**
 * @file benchmark.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

namespace sbody {
	int serialBenchmark();
	int parallelBenchmark(int thread_num);
} // namespace sbody

/**
 * @file numpy.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_NUMPY_HPP
#define SBODY_NUMPY_HPP

#include <array>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "file/writer.hpp"

namespace sbody {
	/**
	 * @brief
	 *
	 */
	class NumPy : public Writer {
	  private:
		/// Dimension of the NumPy file.
		const std::vector<std::size_t> dimension_;
		/// Size of a row
		int row_size_;
		/// Total size of all saved bytes
		int file_size_;

	  public:
		/**
		 * @brief Construct a new NumPy object
		 *
		 * @param file_path
		 * @param dimension
		 */
		NumPy(std::filesystem::path file_path, std::vector<std::size_t> dimension);
		~NumPy();
		void save(const double data[], uint length);
		void save(double datum);
	};
} // namespace sbody
#endif

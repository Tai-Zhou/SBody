/**
 * @file numpy.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "file/numpy.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "file/writer.hpp"
#include "infra/log.hpp"

namespace sbody {
	NumPy::NumPy(std::filesystem::path file_path, std::vector<std::size_t> dimension) : Writer(file_path.replace_extension(".npy"), std::ios::binary | std::ios::out), dimension_(dimension), file_size_(0) {
		row_size_ = 1;
		for (int i : dimension)
			if (i < 1)
				PrintlnError("NumPy file dimension warning (smaller than 1)");
			else
				row_size_ *= i;
		write("\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (                                                                  \n", 128);
	}
	NumPy::~NumPy() {
		int row_number = ceil(file_size_ / row_size_);
		for (int i = row_number * row_size_; i > file_size_; --i)
			write("\0\0\0\0\0\0\0\0", 8);
		seekP(61);
		std::string row_number_string = std::to_string(row_number);
		write(row_number_string.c_str(), row_number_string.length());
		for (int i : dimension_)
			if (i < 1)
				write(", 1", 3);
			else {
				row_number_string = ", " + std::to_string(i);
				write(row_number_string.c_str(), row_number_string.length());
			}
		write("), }", 4);
	}
	void NumPy::save(const double data[], uint length) {
		file_size_ += length;
		write(reinterpret_cast<const char *>(data), 8 * length);
	}
	void NumPy::save(double datum) {
		file_size_ += 1;
		write(reinterpret_cast<const char *>(&datum), 8);
	}
} // namespace sbody

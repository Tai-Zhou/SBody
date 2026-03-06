/**
 * @file csv.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "file/csv.hpp"

#include <fstream>

namespace sbody {
	CSV::CSV(std::filesystem::path file_path, char separator) : Writer(file_path.replace_extension(".csv"), std::ios::out), separator_(separator) {}
	void CSV::save(const double data[], uint length) {
		for (int i = 0; i < length; ++i) {
			std::string datum_s = std::to_string(data[i]);
			write(datum_s.c_str(), datum_s.length());
			put(separator_);
		}
		seekP(-1, std::ios::cur);
		put('\n');
	}
	void CSV::save(const double datum) {
		std::string datum_s = std::to_string(datum);
		write(datum_s.c_str(), datum_s.length());
		put(separator_);
	}
} // namespace sbody

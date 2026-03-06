/**
 * @file writer.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "file/writer.hpp"

#include <filesystem>
#include <fstream>

#include "infra/status.hpp"

namespace sbody {
	Writer::Writer(const std::filesystem::path &file_path, std::ios::openmode mode) {
		if (file_path.stem().empty()) {
			null_output_ = true;
			return;
		}
		if (auto dir = file_path.parent_path(); !std::filesystem::exists(dir))
			std::filesystem::create_directories(dir);
		file_buffer_.open(file_path, mode);
		if (!file_buffer_.is_open()) {
			throw std::ios_base::failure("Writer: failed to open file.");
		}
		file_buffer_.exceptions(std::ios::failbit | std::ios::badbit);
	}
	Writer::~Writer() {
		file_buffer_.close();
	}
	void Writer::seekP(std::streamoff offset, std::ios_base::seekdir dir) {
		if (!null_output_)
			file_buffer_.seekp(offset, dir);
	}
	void Writer::put(char c) {
		if (null_output_)
			return;
		file_buffer_.put(c);
	}
	void Writer::write(const char *s, std::size_t n) {
		if (null_output_)
			return;
		file_buffer_.write(s, n);
	}
} // namespace sbody

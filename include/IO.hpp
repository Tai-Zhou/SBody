/**
 * @file IO.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_IO_HPP
#define SBODY_IO_HPP

#include <array>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <fmt/core.h>
#include <fmt/format.h>
#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/dynamic_progress.hpp>

#include "Utility.hpp"

namespace SBody {
	template <typename Type>
	auto format_as(Type e) {
		return fmt::underlying(e);
	}
	template <typename... Type>
	void PrintlnBold(std::string format, Type &&...args) {
		fmt::println(format, "\033[1m", "\033[0m", std::forward<Type>(args)...);
	}
	template <typename... Type>
	void PrintlnWarning(std::string format, Type &&...args) {
		fmt::println(stdout, "\033[103m[WRN]\033[0m " + format, std::forward<Type>(args)...);
	}
	template <typename... Type>
	void PrintlnError(std::string format, Type &&...args) {
		fmt::println(stdout, "\033[101m[ERR]\033[0m " + format, std::forward<Type>(args)...);
	}
	/**
	 * @brief
	 *
	 */
	class ProgressBar : public indicators::DynamicProgress<indicators::BlockProgressBar> {
	  public:
		ProgressBar(int color, std::string prefix_text) {
			indicators::show_console_cursor(false);
			indicators::BlockProgressBar bar{
				indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
				indicators::option::ForegroundColor{indicators::Color(color)},
				indicators::option::PrefixText{prefix_text},
				indicators::option::ShowElapsedTime{true},
				indicators::option::ShowRemainingTime{true},
			};
			this->push_back(bar);
		}
		~ProgressBar() {
			indicators::show_console_cursor(true);
		}
		/**
		 * @brief Set a specific progress bar as completed.
		 *
		 * @param index index of the progress bar in the container
		 * @param prefix prefix shown in the completed progress bar
		 */
		void SetComplete(int index, std::string prefix) {
			this->operator[](index).set_progress(100.);
			this->operator[](index).set_option(indicators::option::PrefixText{prefix});
			this->operator[](index).mark_as_completed();
		}
	};
	/**
	 * @brief
	 *
	 */
	class File {
	  protected:
		bool null_output_;
		/// File buffer to the file
		std::ofstream file_buffer_;
		/**
		 * @brief Construct a new file object
		 *
		 * @param file_path
		 * @param mode
		 */
		File(const std::filesystem::path &file_path, std::ios::openmode mode) {
			if (file_path.stem().empty()) {
				null_output_ = true;
				return;
			}
			if (auto dir = file_path.parent_path(); !std::filesystem::exists(dir))
				std::filesystem::create_directories(dir);
			file_buffer_.open(file_path, mode);
			if (!file_buffer_.is_open()) {
				PrintlnError("IO file not open");
				throw StatusCode::INPUT_OUTPUT_ERROR;
			}
			file_buffer_.exceptions(std::ios::failbit | std::ios::badbit);
		}
		void SeekP(std::streamoff offset, std::ios_base::seekdir dir = std::ios::beg) {
			if (!null_output_)
				file_buffer_.seekp(offset, dir);
		}
		int Put(char c) {
			if (null_output_)
				return StatusCode::SUCCESS;
			try {
				file_buffer_.put(c);
			} catch (const std::ios_base::failure &e) {
				throw e;
			}
			return StatusCode::SUCCESS;
		}
		int Write(const char *s, std::size_t n) {
			if (null_output_)
				return StatusCode::SUCCESS;
			try {
				file_buffer_.write(s, n);
			} catch (const std::ios_base::failure &e) {
				throw e;
			}
			return StatusCode::SUCCESS;
		}

	  public:
		/// Destroy the file object
		virtual ~File() {
			file_buffer_.close();
		}
		/**
		 * @brief Save the data.
		 *
		 * @param data
		 * @return int
		 */
		virtual int Save(const std::vector<double> &data) = 0;
		/**
		 * @brief Save the data.
		 *
		 * @param data Data to be saved
		 * @param length Length of the data
		 * @return int
		 */
		virtual int Save(const double data[], uint length) = 0;
	};
	/**
	 * @brief
	 *
	 */
	// template <typename Type>
	class NumPy : public File {
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
		NumPy(std::filesystem::path file_path, std::vector<std::size_t> dimension) : File(file_path.replace_extension(".npy"), std::ios::binary | std::ios::out), dimension_(dimension), file_size_(0) {
			row_size_ = 1;
			for (int i : dimension)
				if (i < 1)
					PrintlnError("NumPy file dimension warning (smaller than 1)");
				else
					row_size_ *= i;
			Write("\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (                                                                  \n", 128);
		}
		~NumPy() {
			int row_number = ceil(file_size_ / row_size_);
			for (int i = row_number * row_size_; i > file_size_; --i)
				Write("\0\0\0\0\0\0\0\0", 8);
			SeekP(61);
			std::string row_number_string = std::to_string(row_number);
			Write(row_number_string.c_str(), row_number_string.length());
			for (int i : dimension_)
				if (i < 1)
					Write(", 1", 3);
				else {
					row_number_string = ", " + std::to_string(i);
					Write(row_number_string.c_str(), row_number_string.length());
				}
			Write("), }", 4);
		}
		template <std::size_t N>
		int Save(const std::array<double, N> &data) {
			file_size_ += N;
			return Write(reinterpret_cast<const char *>(data.data()), 8 * N);
		}
		template <std::size_t N>
		int Save(const boost::numeric::ublas::bounded_vector<double, N> &data) {
			file_size_ += N;
			return Write(reinterpret_cast<const char *>(data.data().begin()), 8 * N);
		}
		int Save(const std::vector<double> &data) {
			file_size_ += data.size();
			return Write(reinterpret_cast<const char *>(data.data()), 8 * data.size());
		}
		int Save(const double *data, uint length) {
			file_size_ += length;
			return Write(reinterpret_cast<const char *>(data), 8 * length);
		}
		int Save(double data) {
			file_size_ += 1;
			return Write(reinterpret_cast<const char *>(&data), 8);
		}
	};
	/**
	 * @brief
	 *
	 */
	class CSV : public File {
	  private:
		/// Separator of the file.
		char separator_;

	  public:
		/**
		 * @brief Construct a new CSV object.
		 *
		 * @param file_path
		 * @param separator Separator of the file. If `separator == '\t'`, the file would be a Tab-Separated Values (TSV).
		 */
		CSV(std::filesystem::path file_path, char separator = ',') : File(file_path.replace_extension(".npy"), std::ios::out), separator_(separator) {}
		int Save(const std::vector<double> &data) {
			for (const double &element : data) {
				std::string elements = std::to_string(element);
				Write(elements.c_str(), elements.length());
				Put(separator_);
			}
			SeekP(-1, std::ios::cur);
			return Put('\n');
		}
	};
} // namespace SBody

#endif

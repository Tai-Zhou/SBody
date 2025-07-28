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

#ifndef SBODY_IO_H
#define SBODY_IO_H

#include <array>
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
		ProgressBar() {
			indicators::show_console_cursor(false);
			indicators::BlockProgressBar bar{
				indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
				indicators::option::ForegroundColor{indicators::Color(4)},
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
		/// File buffer to the file
		std::filebuf file_buffer_;
		/**
		 * @brief Construct a new file object
		 *
		 * @param file_name
		 * @param mode
		 */
		File(std::string file_name, std::ios::openmode mode) {
			if (file_name.empty())
				file_buffer_.open("/dev/null", mode);
			else
				file_buffer_.open("../data/" + file_name, mode);
			if (!file_buffer_.is_open()) {
				PrintlnError("IO file not open");
				exit(1);
			}
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
		const std::vector<int> dimension_;
		/// Size of a row
		int row_size_;
		/// Total size of all saved bytes
		int file_size_;

	  public:
		/**
		 * @brief Construct a new NumPy object
		 *
		 * @param file_name
		 * @param dimension
		 */
		NumPy(std::string file_name, std::vector<int> dimension) : File(file_name + ".npy", std::ios::binary | std::ios::out), dimension_(dimension), file_size_(0) {
			row_size_ = 1;
			for (int i : dimension)
				if (i < 1)
					PrintlnError("NumPy file dimension warning (smaller than 1)");
				else
					row_size_ *= i;
			file_buffer_.sputn("\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (                                                                  \n", 128);
		}
		~NumPy() {
			int row_number = ceil(file_size_ / row_size_);
			for (int i = row_number * row_size_; i > file_size_; --i)
				file_buffer_.sputn("\0\0\0\0\0\0\0\0", 8);
			file_buffer_.pubseekpos(61);
			std::string row_number_string = std::to_string(row_number);
			file_buffer_.sputn(row_number_string.c_str(), row_number_string.length());
			for (int i : dimension_)
				if (i < 1)
					file_buffer_.sputn(", 1", 3);
				else {
					row_number_string = ", " + std::to_string(i);
					file_buffer_.sputn(row_number_string.c_str(), row_number_string.length());
				}
			file_buffer_.sputn("), }", 4);
		}
		template <std::size_t N>
		int Save(const std::array<double, N> &data) {
			file_size_ += N;
			file_buffer_.sputn(reinterpret_cast<const char *>(data.data()), 8 * N);
			return Status::SUCCESS;
		}
		int Save(const std::vector<double> &data) {
			file_size_ += data.size();
			file_buffer_.sputn(reinterpret_cast<const char *>(data.data()), 8 * data.size());
			return Status::SUCCESS;
		}
		int Save(const double *data, uint length) {
			file_size_ += length;
			file_buffer_.sputn(reinterpret_cast<const char *>(data), 8 * length);
			return Status::SUCCESS;
		}
		int Save(double data) {
			file_size_ += 1;
			file_buffer_.sputn(reinterpret_cast<const char *>(&data), 8);
			return Status::SUCCESS;
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
		 * @param file_name
		 * @param separator Separator of the file. If `separator == '\t'`, the file would be a Tab-Separated Values (TSV).
		 */
		CSV(std::string file_name, char separator = ',') : File(file_name + ".csv", std::ios::out), separator_(separator) {}
		int Save(const std::vector<double> &data) {
			for (const double &element : data) {
				std::string elements = std::to_string(element);
				file_buffer_.sputn(elements.c_str(), elements.length());
				file_buffer_.sputc(separator_);
			}
			file_buffer_.pubseekoff(-1, std::ios::cur);
			file_buffer_.sputc('\n');
			return Status::SUCCESS;
		}
	};
} // namespace SBody

#endif

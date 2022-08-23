/**
 * @file IO.h
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

#include <fstream>
#include <string>
#include <vector>

#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/dynamic_progress.hpp>

namespace SBody {
	/**
	 * @brief
	 *
	 */
	class ProgressBar {
	  public:
		/// Display of progress bars.
		static bool display_;
		/// Progress bar container object
		static indicators::DynamicProgress<indicators::BlockProgressBar> bars_;
		/**
		 * @brief Set a specific progress bar as completed.
		 *
		 * @param index index of the progress bar in the container
		 * @param prefix prefix shown in the completed progress bar
		 */
		static void SetComplete(int index, std::string prefix);
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
		File(std::string file_name, std::ios::openmode mode);

	  public:
		/// Destroy the file object
		virtual ~File();
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
		virtual int Save(const double *data, const int length) = 0;
	};
	/**
	 * @brief
	 *
	 */
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
		NumPy(std::string file_name, std::vector<int> dimension);
		/// Destroy the NumPy object
		~NumPy();
		int Save(const std::vector<double> &data);
		int Save(const double *data, const int length);
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
		CSV(std::string file_name, char separator = ',');
		int Save(const std::vector<double> &data);
	};
#ifdef WITH_CFITSIO
	/**
	 * @brief
	 *
	 */
	class FITS : public file {
	  public:
		/**
		 * @brief Construct a new FITS object
		 *
		 * @param file_name
		 */
		FITS(std::string file_name);
		int Save(const std::vector<double> &data);
	};
#endif
} // namespace SBody

#endif

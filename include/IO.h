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
	class file {
	  protected:
		/// File buffer to the file
		std::filebuf fileBuf;
		/**
		 * @brief Construct a new file object
		 *
		 * @param fileName
		 * @param mode
		 */
		file(std::string fileName, std::ios::openmode mode);

	  public:
		/// Destroy the file object
		virtual ~file();
		/**
		 * @brief Save the coming data.
		 *
		 * @param data
		 * @return int
		 */
		virtual int save(const std::vector<double> &data) = 0;
	};
	/**
	 * @brief
	 *
	 */
	class NumPy : public file {
	  private:
		/// Dimension of the NumPy file.
		const std::vector<int> dim;
		/// Size of a row
		int dimSize;
		/// Total size of all saved bytes
		int fileSize;

	  public:
		/**
		 * @brief Construct a new NumPy object
		 *
		 * @param fileName
		 * @param dimension
		 */
		NumPy(std::string fileName, std::vector<int> dimension);
		/// Destroy the NumPy object
		~NumPy();
		/**
		 * @brief Save the coming data.
		 *
		 * @param data
		 * @return int
		 */
		int save(const std::vector<double> &data);
	};
	/**
	 * @brief
	 *
	 */
	class CSV : public file {
	  private:
		/// Separator of the file.
		char sep;

	  public:
		/**
		 * @brief Construct a new CSV object.
		 *
		 * @param fileName
		 * @param separator Separator of the file. If `sep == '\t'`, the file would be a Tab-Separated Values (TSV).
		 */
		CSV(std::string fileName, char separator = ',');
		/**
		 * @brief Save the coming data.
		 *
		 * @param data
		 * @return int
		 */
		int save(const std::vector<double> &data);
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
		 * @param fileName
		 */
		FITS(std::string fileName);
		/**
		 * @brief Save the coming data.
		 *
		 * @param data
		 * @return int
		 */
		int save(const std::vector<double> &data);
	};
#endif
} // namespace SBody

#endif

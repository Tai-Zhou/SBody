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
	namespace IO {
		/// global switch to control the display of progress bars. `0` for not display, `1` for display.
		extern int displayProgressBar;
		/// Progress bar container object
		extern indicators::DynamicProgress<indicators::BlockProgressBar> progressBars;
		/**
		 * @brief Set a specific progress bar as completed.
		 *
		 * @param index index of the progress bar in the container
		 * @param prefix prefix shown in the progress bar
		 */
		void progressBarComplete(int index, std::string prefix);
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
		class CSV : public file {
		  private:
			/// Separator of the file.
			char sep;

		  public:
			/**
			 * @brief Construct a new CSV object
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
		int save();
		int load();
	} // namespace IO
} // namespace SBody

#endif

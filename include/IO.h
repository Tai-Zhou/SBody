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
		extern int displayProgressBar;
		extern indicators::DynamicProgress<indicators::BlockProgressBar> progressBars;
		void progressBarComplete(int index, std::string prefix);
		class file {
		  protected:
			std::filebuf fileBuf;
			file(std::string fileName, std::ios::openmode mode);

		  public:
			virtual ~file();
			virtual int save(const std::vector<double> &data) = 0;
		};
		class NumPy : public file {
		  private:
			const std::vector<int> dim;
			int dimSize;
			int fileSize;

		  public:
			NumPy(std::string fileName, std::vector<int> dimension);
			~NumPy();
			int save(const std::vector<double> &data);
		};
		class CSV : public file {
		  private:
			char sep;

		  public:
			CSV(std::string fileName, char separator = ',');
			int save(const std::vector<double> &data);
		};
#ifdef WITH_CFITSIO
		class FITS : public file {
		  public:
			FITS(std::string fileName);
			int save(const std::vector<double> &data);
		};
#endif
		int save();
		int load();
	} // namespace IO
} // namespace SBody

#endif

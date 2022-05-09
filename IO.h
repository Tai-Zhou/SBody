#ifndef _IO_H
#define _IO_H

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
			const int column;
			int row;

		  public:
			NumPy(std::string fileName, int columnNumber);
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
		class FITS : public file {
		  public:
			FITS(std::string fileName);
			int save(const std::vector<double> &data);
		};
		int save();
		int load();
	} // namespace IO
} // namespace SBody

#endif

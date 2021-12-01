#ifndef _IO_H
#define _IO_H

#include <fstream>
#include <string>
#include <vector>

#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>

namespace SBody {
	namespace IO {
		extern indicators::BlockProgressBar progressBar;
		class file {
		  protected:
			const std::string fileName;
			file(std::string fileName);
		};
		// T should be int or double
		template <typename T>
		class NumPy : public file {
		  public:
			NumPy(std::string fileName);
			int save(const std::vector<std::vector<T>> &data);
		};
		template <typename T>
		class CSV : public file {
		  public:
			CSV(std::string fileName);
			int save(const std::vector<std::vector<T>> &data);
		};
		template <typename T>
		class FITS : public file {
		  public:
			FITS(std::string fileName);
			int save(const std::vector<std::vector<T>> &data);
		};
		int save();
		int load();
	} // namespace IO
} // namespace SBody

#endif

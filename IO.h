#ifndef _IO_H
#define _IO_H

#include <fstream>
#include <string>
#include <vector>

namespace SBody {
	namespace IO {
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
		int save();
		int load();
	} // namespace IO
} // namespace SBody

#endif

#ifndef _IO_H
#define _IO_H

#include <fstream>
#include <string>
#include <vector>

namespace SBody {
	namespace IO {
		template <typename T>
		int NumPySave(const std::vector<std::vector<T>> &data, std::string fileName);
		template <typename T>
		int CSVSave(const std::vector<std::vector<T>> &data, std::string fileName);
		int save();
		int load();
	} // namespace IO
} // namespace SBody

#endif

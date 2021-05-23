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
		template <typename T>
		class NumPy : public file {
		  public:
			NumPy(std::string fileName);
			int save(const std::vector<std::vector<T>> &data);
		};
		template <typename T>
		int CSVSave(const std::vector<std::vector<T>> &data, std::string fileName) {
			std::ofstream CSVFile("out/" + fileName + ".csv", std::ios::out);
			for (const std::vector<T> &line : data) {
				for (const T element : line)
					CSVFile << element << ',';
				CSVFile.seekp(-1, std::ios::cur);
				CSVFile << std::endl;
			}
			CSVFile.close();
			return 0;
		}
		int save();
		int load();
	} // namespace IO
} // namespace SBody

#endif

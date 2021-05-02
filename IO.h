#ifndef _IO_H
#define _IO_H

#include <fstream>
#include <string>
#include <vector>

namespace IO {
	template <typename T>
	int NumPySave(const std::vector<std::vector<T>> &data, std::string fileName) {
		char NumPyHead[10] = {'\x93', 'N', 'U', 'M', 'P', 'Y', '\x01', '\x00', '\x76', '\x00'};
		std::ofstream NumPyFile(fileName + ".npy", std::ios::binary | std::ios::out);
		NumPyFile.write(NumPyHead, 10);
		NumPyFile << "{'descr': '<" << (sizeof(T) == 4 ? "i4" : "f8")
				  << "', 'fortran_order': False, 'shape': (                                                                  \n";
		for (const std::vector<T> &line : data)
			NumPyFile.write((char *)(&line[0]), sizeof(T) * line.size());
		NumPyFile.seekp(61);
		NumPyFile << data.size() << ", " << data[0].size() << "), }";
		NumPyFile.close();
		return 0;
	}
	template <typename T>
	int CSVSave(const std::vector<std::vector<T>> &data, std::string fileName) {
		std::ofstream CSVFile(fileName + ".csv", std::ios::out);
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
#endif

#include "IO.h"

#include <array>
#include <fstream>

#include <fitsio.h>

using namespace std;

namespace SBody {
	namespace IO {
		template <typename T>
		int NumPySave(const vector<vector<T>> &data, string fileName) {
			char NumPyHead[10] = {'\x93', 'N', 'U', 'M', 'P', 'Y', '\x01', '\x00', '\x76', '\x00'};
			ofstream NumPyFile(fileName + ".npy", ios::binary | ios::out);
			NumPyFile.write(NumPyHead, 10);
			NumPyFile << "{'descr': '<" << (sizeof(T) == 4 ? "i4" : "f8")
					  << "', 'fortran_order': False, 'shape': (                                                                  \n";
			for (const vector<T> &line : data)
				NumPyFile.write((char *)(&line[0]), sizeof(T) * line.size());
			NumPyFile.seekp(61);
			NumPyFile << data.size() << ", " << data[0].size() << "), }";
			NumPyFile.close();
			return 0;
		}
		template int NumPySave(const vector<vector<double>> &data, string fileName);
		template int NumPySave(const vector<vector<int>> &data, string fileName);
		template <typename T>
		int NumPySave(const vector<T> &data, string fileName) {
			char NumPyHead[10] = {'\x93', 'N', 'U', 'M', 'P', 'Y', '\x01', '\x00', '\x76', '\x00'};
			ofstream NumPyFile(fileName + ".npy", ios::binary | ios::out);
			NumPyFile.write(NumPyHead, 10);
			NumPyFile << "{'descr': '<" << (sizeof(T) == 4 ? "i4" : "f8")
					  << "', 'fortran_order': False, 'shape': (                                                                  \n";
			for (const T &line : data)
				NumPyFile.write((char *)(line.data()), sizeof(double) * line.size());
			NumPyFile.seekp(61);
			NumPyFile << data.size() << ", " << data[0].size() << "), }";
			NumPyFile.close();
			return 0;
		}
		template int NumPySave(const vector<array<double, 3>> &data, string fileName);
		template <typename T>
		int CSVSave(const vector<vector<T>> &data, string fileName) {
			ofstream CSVFile(fileName + ".csv", ios::out);
			for (const vector<T> &line : data) {
				for (const T element : line)
					CSVFile << element << ',';
				CSVFile.seekp(-1, ios::cur);
				CSVFile << endl;
			}
			CSVFile.close();
			return 0;
		}
		template int CSVSave(const vector<vector<double>> &data, string fileName);
		template int CSVSave(const vector<vector<int>> &data, string fileName);
		int save() {
			return 0;
		}
		int load() {
			return 0;
		}
	} // namespace IO
} // namespace SBody

#include "IO.h"

#include <fitsio.h>

using namespace std;

namespace SBody {
	namespace IO {
		file::file(string fileName) : fileName(fileName) {}
		template <typename T>
		NumPy<T>::NumPy(string fileName) : file(fileName) {}
		template <typename T>
		int NumPy<T>::save(const vector<vector<T>> &data) {
			char NumPyHead[10] = {'\x93', 'N', 'U', 'M', 'P', 'Y', '\x01', '\x00', '\x76', '\x00'};
			ofstream NumPyFile("out/" + fileName + ".npy", ios::binary | ios::out);
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
		template class NumPy<int>;
		template class NumPy<double>;
		int save() {
			return 0;
		}
		int load() {
			return 0;
		}
	} // namespace IO
} // namespace SBody

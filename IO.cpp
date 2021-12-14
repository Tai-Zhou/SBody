#include "IO.h"

#include <fitsio.h>

using namespace std;

namespace SBody {
	namespace IO {
		indicators::BlockProgressBar progressBar{
			indicators::option::BarWidth{80},
			indicators::option::Start{"["},
			indicators::option::End{"]"},
			indicators::option::PrefixText{"Calculating..."},
			indicators::option::ForegroundColor{indicators::Color::cyan},
			indicators::option::ShowElapsedTime{true},
			indicators::option::ShowRemainingTime{true},
			indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
		file::file(string fileName) : fileName(fileName) {}
		template <typename T>
		NumPy<T>::NumPy(string fileName) : file(fileName) {}
		template <typename T>
		int NumPy<T>::save(const vector<vector<T>> &data) {
			ofstream NumPyFile("out/" + fileName + ".npy", ios::binary | ios::out);
			NumPyFile.write("\x93NUMPY\x01\x00\x76\x00", 10);
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
		template <typename T>
		CSV<T>::CSV(string fileName) : file(fileName) {}
		template <typename T>
		int CSV<T>::save(const vector<vector<T>> &data) {
			ofstream CSVFile("out/" + fileName + ".csv", ios::out);
			for (const vector<T> &line : data) {
				for (const T element : line)
					CSVFile << element << ',';
				CSVFile.seekp(-1, ios::cur);
				CSVFile << endl;
			}
			CSVFile.close();
			return 0;
		}
		template class CSV<int>;
		template class CSV<double>;
		template <typename T>
		FITS<T>::FITS(string fileName) : file(fileName) {}
		template <typename T>
		int FITS<T>::save(const vector<vector<T>> &data) {
			fitsfile *fptr;
			char card[FLEN_CARD];
			int status = 0, nkeys;

			/*fits_open_file(&fptr, argv[1], READONLY, &status);
			fits_get_hdrspace(fptr, &nkeys, NULL, &status);

			for (int i = 1; i <= nkeys; i++) {
				fits_read_record(fptr, i, card, &status);
				printf("%s\n", card);
			}
			fits_close_file(fptr, &status);*/
			if (status)
				fits_report_error(stderr, status);
			return (status);
		}
		template class FITS<int>;
		template class FITS<double>;
		int save() {
			return 0;
		}
		int load() {
			return 0;
		}
	} // namespace IO
} // namespace SBody

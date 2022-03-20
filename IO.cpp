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
		file::file(string fileName, ios::openmode mode) {
			fileBuf.open("out/" + fileName, mode);
		}
		file::~file() {
			fileBuf.close();
		}
		template <typename T>
		NumPy<T>::NumPy(string fileName, int columnNumber) : file(fileName + ".npy", ios::binary | ios::out), column(columnNumber), row(0) {
			fileBuf.sputn("\x93NUMPY\x01\x00\x76\x00", 10);
			fileBuf.sputn("{'descr': '<", 12);
			fileBuf.sputn(sizeof(T) == 4 ? "i4" : "f8", 2);
			fileBuf.sputn("', 'fortran_order': False, 'shape': (                                                                  \n", 104);
		}
		template <typename T>
		NumPy<T>::~NumPy() {
			fileBuf.pubseekpos(61);
			string rows = to_string(row), columns = to_string(column);
			fileBuf.sputn((rows + ", " + columns + "), }").c_str(), rows.length() + columns.length() + 6);
		}
		template <typename T>
		int NumPy<T>::save(const vector<T> &data) {
			++row;
			fileBuf.sputn(reinterpret_cast<const char *>(data.data()), sizeof(T) * column);
			return 0;
		}
		template class NumPy<int>;
		template class NumPy<double>;
		template <typename T>
		CSV<T>::CSV(string fileName, char separator) : file(fileName + ".csv", ios::out), sep(separator) {}
		template <typename T>
		int CSV<T>::save(const vector<T> &data) {
			for (const T &element : data) {
				string elements = to_string(element);
				fileBuf.sputn(elements.c_str(), elements.length());
				fileBuf.sputc(sep);
			}
			fileBuf.pubseekoff(-1, ios::cur);
			fileBuf.sputc('\n');
			return 0;
		}
		template class CSV<int>;
		template class CSV<double>;
		template <typename T>
		FITS<T>::FITS(string fileName) : file(fileName + ".fits", ios::binary | ios::out) {}
		template <typename T>
		int FITS<T>::save(const vector<T> &data) {
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

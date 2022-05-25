#include "IO.h"

#include <fitsio.h>

using namespace std;
using namespace indicators;

namespace SBody {
	namespace IO {
		int displayProgressBar = 0;
		BlockProgressBar mainProgressBar{
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{vector<FontStyle>{FontStyle::bold}}};
		DynamicProgress<BlockProgressBar> progressBars(mainProgressBar);
		void progressBarComplete(int index, string prefix) {
			progressBars[index].set_option(option::PrefixText{prefix});
			progressBars[index].mark_as_completed();
		}
		file::file(string fileName, ios::openmode mode) {
			fileBuf.open("../data" + fileName, mode);
			if (!fileBuf.is_open()) {
				cerr << "[!] IO file not open" << endl;
				exit(1);
			}
		}
		file::~file() {
			fileBuf.close();
		}
		NumPy::NumPy(string fileName, int columnNumber) : file(fileName + ".npy", ios::binary | ios::out), column(columnNumber), row(0) {
			fileBuf.sputn("\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (                                                                  \n", 128);
		}
		NumPy::~NumPy() {
			fileBuf.pubseekpos(61);
			string rows = to_string(row), columns = to_string(column);
			fileBuf.sputn((rows + ", " + columns + "), }").c_str(), rows.length() + columns.length() + 6);
		}
		int NumPy::save(const vector<double> &data) {
			++row;
			fileBuf.sputn(reinterpret_cast<const char *>(data.data()), 8 * column);
			return 0;
		}
		CSV::CSV(string fileName, char separator) : file(fileName + ".csv", ios::out), sep(separator) {}
		int CSV::save(const vector<double> &data) {
			for (const double &element : data) {
				string elements = to_string(element);
				fileBuf.sputn(elements.c_str(), elements.length());
				fileBuf.sputc(sep);
			}
			fileBuf.pubseekoff(-1, ios::cur);
			fileBuf.sputc('\n');
			return 0;
		}
		FITS::FITS(string fileName) : file(fileName + ".fits", ios::binary | ios::out) {}
		int FITS::save(const vector<double> &data) {
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
		int save() {
			return 0;
		}
		int load() {
			return 0;
		}
	} // namespace IO
} // namespace SBody
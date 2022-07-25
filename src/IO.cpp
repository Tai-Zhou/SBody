/**
 * @file IO.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "IO.h"

#ifdef WITH_CFITSIO
#include <fitsio.h>
#endif

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
		void progressBarComplete(int index, std::string prefix) {
			progressBars[index].set_progress(100.);
			progressBars[index].set_option(option::PrefixText{prefix});
			progressBars[index].mark_as_completed();
		}
		file::file(string fileName, ios::openmode mode) {
			if (fileName.empty())
				fileBuf.open("/dev/null", mode);
			else
				fileBuf.open("../data/" + fileName, mode);
			if (!fileBuf.is_open()) {
				cerr << "[!] IO file not open" << endl;
				exit(1);
			}
		}
		file::~file() {
			fileBuf.close();
		}
		NumPy::NumPy(string fileName, vector<int> dimension) : file(fileName + ".npy", ios::binary | ios::out), dim(dimension), fileSize(0) {
			dimSize = 1;
			for (int i : dim)
				if (i < 1)
					cerr << "[?] NumPy file dimension warning (smaller than 1)" << endl;
				else
					dimSize *= i;
			fileBuf.sputn("\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (                                                                  \n", 128);
		}
		NumPy::~NumPy() {
			int dimNum = ceil(fileSize / dimSize);
			for (int i = dimNum * dimSize; i > fileSize; --i)
				fileBuf.sputn("\0\0\0\0\0\0\0\0", 8);
			fileBuf.pubseekpos(61);
			string dimNumS = to_string(dimNum);
			fileBuf.sputn(dimNumS.c_str(), dimNumS.length());
			for (int i : dim)
				if (i < 1)
					fileBuf.sputn(", 1", 3);
				else {
					dimNumS = ", " + to_string(i);
					fileBuf.sputn(dimNumS.c_str(), dimNumS.length());
				}
			fileBuf.sputn("), }", 4);
		}
		int NumPy::save(const vector<double> &data) {
			fileSize += data.size();
			fileBuf.sputn(reinterpret_cast<const char *>(data.data()), 8 * data.size());
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
#ifdef WITH_CFITSIO
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
#endif
		int save() {
			return 0;
		}
		int load() {
			return 0;
		}
	} // namespace IO
} // namespace SBody

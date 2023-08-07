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

#include <fmt/core.h>

#ifdef WITH_CFITSIO
#include <fitsio.h>
#endif

using namespace std;

namespace SBody {
	bool ProgressBar::display_ = false;
	indicators::BlockProgressBar main_bar{
		indicators::option::ShowElapsedTime{true},
		indicators::option::ShowRemainingTime{true},
		indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
	indicators::DynamicProgress<indicators::BlockProgressBar> ProgressBar::bars_(main_bar);
	void ProgressBar::SetComplete(int index, std::string prefix) {
		bars_[index].set_progress(100.);
		bars_[index].set_option(indicators::option::PrefixText{prefix});
		bars_[index].mark_as_completed();
	}
	File::File(string file_name, ios::openmode mode) {
		if (file_name.empty())
			file_buffer_.open("/dev/null", mode);
		else
			file_buffer_.open("../data/" + file_name, mode);
		if (!file_buffer_.is_open()) {
			PrintlnError("IO file not open");
			exit(1);
		}
	}
	File::~File() {
		file_buffer_.close();
	}
	NumPy::NumPy(string file_name, vector<int> dimension) : File(file_name + ".npy", ios::binary | ios::out), dimension_(dimension), file_size_(0) {
		row_size_ = 1;
		for (int i : dimension)
			if (i < 1)
				PrintlnError("NumPy file dimension warning (smaller than 1)");
			else
				row_size_ *= i;
		file_buffer_.sputn("\x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False, 'shape': (                                                                  \n", 128);
	}
	NumPy::~NumPy() {
		int row_number = ceil(file_size_ / row_size_);
		for (int i = row_number * row_size_; i > file_size_; --i)
			file_buffer_.sputn("\0\0\0\0\0\0\0\0", 8);
		file_buffer_.pubseekpos(61);
		string row_number_string = to_string(row_number);
		file_buffer_.sputn(row_number_string.c_str(), row_number_string.length());
		for (int i : dimension_)
			if (i < 1)
				file_buffer_.sputn(", 1", 3);
			else {
				row_number_string = ", " + to_string(i);
				file_buffer_.sputn(row_number_string.c_str(), row_number_string.length());
			}
		file_buffer_.sputn("), }", 4);
	}
	int NumPy::Save(const vector<double> &data) {
		file_size_ += data.size();
		file_buffer_.sputn(reinterpret_cast<const char *>(data.data()), 8 * data.size());
		return 0;
	}
	int NumPy::Save(const double data[], uint length) {
		file_size_ += length;
		file_buffer_.sputn(reinterpret_cast<const char *>(data), 8 * length);
		return 0;
	}
	CSV::CSV(string file_name, char separator) : File(file_name + ".csv", ios::out), separator_(separator) {}
	int CSV::Save(const vector<double> &data) {
		for (const double &element : data) {
			string elements = to_string(element);
			file_buffer_.sputn(elements.c_str(), elements.length());
			file_buffer_.sputc(separator_);
		}
		file_buffer_.pubseekoff(-1, ios::cur);
		file_buffer_.sputc('\n');
		return 0;
	}
#ifdef WITH_CFITSIO
	FITS::FITS(string file_name) : file(file_name + ".fits", ios::binary | ios::out) {}
	int FITS::Save(const vector<double> &data) {
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
} // namespace SBody

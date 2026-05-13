/**
 * @file file.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "file/file.hpp"

#include "file/csv.hpp"
#include "file/numpy.hpp"

namespace sbody {

FileFormat getFileFormat(const std::string &name) {
    std::string lower_name;
    std::transform(name.begin(), name.end(), std::back_inserter(lower_name), ::tolower);
    if (lower_name == "numpy")
        return FileFormat::NUMPY;
    else if (lower_name == "csv")
        return FileFormat::CSV;
    throw std::invalid_argument("Invalid file format: " + name);
}

File::File(FileFormat format, std::filesystem::path file_path, std::vector<std::size_t> feature_shape) : format_(format) {
    if (format_ == FileFormat::NUMPY)
        writer_ = std::make_unique<NumPy>(file_path, feature_shape);
    else if (format_ == FileFormat::CSV)
        writer_ = std::make_unique<CSV>(file_path);
}
void File::save(const std::vector<double> &data) { writer_->save(data.data(), data.size()); }
void File::save(const double data[], uint length) { writer_->save(data, length); }
void File::save(double datum) { writer_->save(datum); }
} // namespace sbody

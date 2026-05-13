/**
 * @file file.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include "file/writer.hpp"

namespace sbody {
enum class FileFormat {
    NUMPY,
    CSV,
};

/**
 * @brief Get the File Format enum
 *
 * @param name
 * @return FileFormat
 */
FileFormat getFileFormat(const std::string &name);

class File {
  protected:
    FileFormat format_;
    std::unique_ptr<Writer> writer_;

  public:
    File(FileFormat format, std::filesystem::path file_path, std::vector<std::size_t> feature_shape);
    template <std::size_t N>
    void save(const std::array<double, N> &data) {
        writer_->save(data.data(), N);
    }
    template <std::size_t N>
    void save(const boost::numeric::ublas::bounded_vector<double, N> &data) {
        static_assert(N > 0, "N must be greater than 0");
        writer_->save(&data[0], N);
    }
    void save(const std::vector<double> &data);
    void save(const double data[], uint length);
    void save(double datum);
};
} // namespace sbody

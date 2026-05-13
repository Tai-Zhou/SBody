/**
 * @file writer.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <filesystem>
#include <fstream>

namespace sbody {
/**
 * @brief
 *
 */
class Writer {
  protected:
    bool null_output_;
    /// File buffer to the file
    std::ofstream file_buffer_;
    /**
     * @brief Construct a new file object
     *
     * @param file_path
     * @param mode
     */
    Writer(const std::filesystem::path &file_path, std::ios::openmode mode);
    void seekP(std::streamoff offset, std::ios_base::seekdir dir = std::ios::beg);
    void put(char c);
    void write(const char *s, std::size_t n);

  public:
    /// Destroy the file object
    virtual ~Writer();
    /**
     * @brief Save the data.
     *
     * @param data Data to be saved
     * @param length Length of the data
     */
    virtual void save(const double data[], uint length) = 0;
    /**
     * @brief Save the data.
     *
     * @param datum Data to be saved
     */
    virtual void save(const double datum) = 0;
};
} // namespace sbody

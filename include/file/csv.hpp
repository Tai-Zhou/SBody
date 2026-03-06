/**
 * @file csv.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_CSV_HPP
#define SBODY_CSV_HPP

#include "file/writer.hpp"

namespace sbody {
	/**
	 * @brief
	 *
	 */
	class CSV : public Writer {
	  private:
		/// Separator of the file.
		char separator_;

	  public:
		/**
		 * @brief Construct a new CSV object.
		 *
		 * @param file_path
		 * @param separator Separator of the file. If `separator == '\t'`, the file would be a Tab-Separated Values (TSV).
		 */
		CSV(std::filesystem::path file_path, char separator = ',');
		void save(const double data[], uint length) override;
		void save(const double datum) override;
	};
} // namespace sbody
#endif

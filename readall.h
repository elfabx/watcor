#ifndef READALL_H
#define READALL_H
#include <vector>
#include <string>

/** \defgroup readall File reading utility
 * @{
 */

//! Read all lines of a (text) file into a vector.
/**
 * Each line of the file becomes a string in the vector. The lines are added to
 * the vector, which is not cleared. Consecutive calls to different files
 * can thus be used to concatenate the files to form a single input vector.
 *
 * \param name file name to read
 * \param lines handle to the vector to which lines are added
 * \return the number of lines successfully read
 * \throw runtime_error if file cannot be opened/read
 */
long readall(std::string name, std::vector<std::string> &lines);

/**@}*/

#endif

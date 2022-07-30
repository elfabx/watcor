#ifndef GRO_H
#define GRO_H
#include "model.h"
#include <vector>
#include <string>
#include <iostream>

//! Modify water molecules in a gro file to match model wm (see model.h)
/**
  *  \param os the output stream to write results to
  *  \param lines vector made up of the lines of the input gro file
  *  \param wm the water model to be used in the output
  *  \return the number of molecules changed
  *  \sa model.h
  *  \throws gro_error indicates error in parsing the input file
*/
int process_gro(std::ostream &os, const std::vector<std::string> &lines,
                const model &wm);

//! Exception class to reflect error in parsing gro file
class gro_error: public std::runtime_error {
public:
    //! Constructor of gro_error class
    /**
     * \param msg description of the error
     * \param line the line in which the error occurred
     */
    gro_error(const std::string & msg = "", const std::string & line = "") :
        std::runtime_error(msg+"; current line:\n"+line) {}
};

#endif

#include "readall.h"
#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>

long readall(std::string name, std::vector<std::string> &lines) {

    std::ifstream inp;

    inp.open(name);
    if (!inp.good()) {
        std::string msg{ "cannot open '" + name + "': " };
        msg += strerror(errno);
        throw(std::runtime_error(msg));
    }
    
    long n{ 0 }; // nr. of lines read
    std::string l; // line buffer
    while (std::getline(inp,l)) {
        lines.push_back(l);
        ++n;
    }
    
    if (!inp.eof()) {
        std::string msg { "error while reading '" + name + "': " };
        msg += strerror(errno);
        throw(std::runtime_error(msg));
    }
    
    return n;

}

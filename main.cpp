#include "model.h"
#include "readall.h"
#include "gro.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cerrno>
#include <cstring>
#include <vector>

//! Print unix style usage information
/**
 * \param a  name of the current executable
*/
void print_help(const std::string a) {
    std::cout << "Usage: " << a << " [-m model] infile [outfile]\n";
    std::cout << "Convert MD coordinate file for use with a different ";
    std::cout << "water model.\n\n";
    std::cout << "Supported models:\n";
    std::vector<std::string> m = model::catalog();
    for (auto i = m.begin(); i != m.end(); ++i) {
        std::cout << "  " << *i;
        if (i == m.begin()) { std::cout << " (default)"; }
        std::cout << std::endl;
    }
}

// don't use enum class, so that it works as int
enum return_t {
    RET_OK = 0,
    RET_COMMAND_ERROR,
    RET_FILE_IO_ERROR,
    RET_FILE_FORMAT_ERROR
};

int main(int argc, char **argv) {

    // give help if requested
    
    if (argc == 1) { print_help(argv[0]); return RET_OK; }
    std::string arg{ argv[1] };
    if (arg == "-h" || arg == "--help") { print_help(argv[0]); return RET_OK; }
    
    // select water model
    
    int n{ 1 }; // index of current command line argument
    model m; // selected model
    if (arg == "-m") {
        if (argc < 4) { print_help(argv[0]); return RET_COMMAND_ERROR; }
        std::string nm{ argv[2] };
        if (!m.initialise(nm)) {
            print_help(argv[0]);
            return RET_COMMAND_ERROR;
        }
        n = 3;
    } else {
        m.initialise(0); // default model is the first
    }

    // open input file
    
    std::vector<std::string> lines{};
    long rd; // lines read
    try {
        rd = readall(argv[n],lines);
    }
    catch (const std::runtime_error & e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        return RET_FILE_IO_ERROR;
    }
    
    //basic format check
    
    if (rd < 1) {
        std::cerr << argv[0] << ": cannot process input: '";
        std::cerr << argv[n] << "' is empty" << std::endl;
        return RET_FILE_FORMAT_ERROR;
    }
    
    // open output file or use stdout if none given
    
    ++n;
    std::ofstream of;
    std::ostream *out{ &std::cout };
    if (n < argc) {
        of.open(argv[n]);
        if (!of.good()) {
            std::cerr << argv[0] << ": cannot open '" << argv[n];
            std::cerr << "': " << strerror(errno) << std::endl;
            return RET_FILE_IO_ERROR;
        }
        out = &of;
    }
    
    // produce output
    
    int wf; // water mols. found & modified
    try {
        wf = process_gro(*out,lines,m);
    }
    catch(const gro_error & e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        std::cerr << "  in '" << argv[n-1] << "'" << std::endl;
        return RET_FILE_FORMAT_ERROR;
    }
    
    std::clog << "Processed " << wf << " water molecules.\n";
    
    // close & check output for errors
    
    out->flush();  // error state only correct after flush (also done by endl)
    if (of.is_open()) { of.close(); }
    // close flushes, but explicit flush still needed for stdout
    if (!*out) {
        std::cerr << argv[0] << ": error writing results" << std::endl;
        return RET_FILE_IO_ERROR;
    }
    
    return RET_OK;

}

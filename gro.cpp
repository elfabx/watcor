#include "gro.h"
#include <cstdio>


// extract atom name from .gro atom line; remove spaces
std::string atom_name(const std::string l) {
    if (l.length() < 15) {
        std::string msg{ "file format error (atom name)" };
        throw(gro_error(msg,l));
    }
    std::string nm{ l.substr(10,5) };
    nm.erase(0,nm.find_first_not_of(" "));
    nm.erase(nm.find_last_not_of(" ")+1);
    return nm;
}

// return only first 2 chars of atom name as uppercase
// (done to simplify comparisons)
std::string standardise(const std::string n) {
    std::string v{ "" };
    for (auto c: n) {
        if (c >= 'a' && c <= 'z') {
            c += 'A' - 'a';            
        }
        v += c;
    }
    if (v.size() < 2) { return v; }
    return v.substr(0,2);
}

// return coordinates in Angstroms from .gro atom line
void coordinates(const std::string l, double &x, double &y, double &z) {
    std::string msg{ "file format error (coordinates)" };
    if (l.length() < 44) {
        throw(gro_error(msg,l));
    }
    try {
        x = std::stod(l.substr(20,8))*10.0;
        y = std::stod(l.substr(28,8))*10.0;
        z = std::stod(l.substr(36,8))*10.0;
    }
    catch(const std::invalid_argument &e) {
        throw(gro_error(msg,l));
    }
}

// replace the atom counter at the beginning of a .gro atom line by c
// and remove velocities
std::string update_line(std::string l, size_t c) {
    
    if (l.length() < 44) {
        std::string msg { "file format error (atom)" };
        throw(gro_error(msg, l));
    };
    
    char buf[8];
    std::snprintf(buf, 6, "%5ld", c);
    
    return l.substr(0,15) + buf + l.substr(20,24);
    
}

// replace coords in addition to update_line(l,c)
std::string update_line(std::string l, size_t c, double x, double y, double z)
{
    std::string t = update_line(l,c); // fix the counter
    
    char buf[32]; // print new coordinates
    std::snprintf(buf, 25, "%8.3f%8.3f%8.3f", x, y, z);
    
    return t.substr(0,20) + buf;
}

// replace atom name in addition to above
std::string update_line(std::string l, std::string atnam, size_t c,
                        double x, double y, double z) {
    
    std::string t { update_line(l,c,x,y,z) };
    
    char buf[6];
    std::snprintf(buf, 6, "%5s", atnam.data());    
    t.replace(10, 5, buf);
    
    return t;
}

int process_gro(std::ostream &os, const std::vector<std::string> &lines,
                const model &wm) {
    
    // how many atoms will we need for each water molecule
    int model_size{ wm.size() };
    
    size_t n{ lines.size() };

    if (n < 5) {
        std::string msg{ "file too short to contain a water molecule" };
        std::string l{ std::to_string(n) };
        throw(gro_error(msg,l));
    }
    
    size_t na{ 0 }; // number of atoms
    long nas{ 0l };  // signed version for reading

    try {
        nas = std::stol(lines[1]);
    }    
    catch (const std::invalid_argument & e) {
        nas = -1;
    }
    // unreadable or negative nr. of atoms
    if (nas < 0) {
        std::string msg{ "file format error (atom count)" };
        throw(gro_error(msg,lines[1]));
    }

    na = static_cast<size_t>(nas); // safe now: nas >= 0
    
    if (n < na + 2) {
        std::string msg{ "file too short for " };
        msg += std::to_string(na) + " atoms";
        throw(gro_error(msg,lines[1]));
    }

    // count water molecules to figure out how many atoms we will have
    // identify consecutive atoms with names starting OW, HW, HW
    // optionally followed by one or more of MW|LP|EP (last is Amber name)
    size_t cur{ 2 }; // current line: first atom
    std::string an{}; // current atom name
    int nw{ 0 }; // number of water molecules
    int nwa{ 0 }; // number of atoms in water molecules

    // iterate over atom lines (excl. last 2, but should be followed by 2 HW)
    while (cur < na) {
        an = standardise(atom_name(lines[cur]));
 //       std::cout << an << ":" << lines[cur] << "\n";
        if (an == "OW" && standardise(atom_name(lines[cur+1])) == "HW"
                    && standardise(atom_name(lines[cur+2])) == "HW") {
            ++nw;
            nwa += 2; cur += 2; //add 2: the while loop will run at least once
            an = "MW";
            while ((cur < na + 2) && (an=="MW" || an=="LP" || an=="EP")) {
                ++cur; ++nwa;
                an = standardise(atom_name(lines[cur]));
            }
        } else {
            ++cur;
        }
    } 
    // NOTE: we may have 2 trailing atom lines with non-water atoms
    // plus the box, but we can ignore them here (not in printing)
    

    // Now we can write the modified file
    
    int modified{ 0 }; // number of water molecules processed
    
    os << lines[0] << '\n'; // title line written unchanged
    os << na - nwa + nw*model_size << '\n'; // new number of atoms
    
    cur = 2; // rewind for printing of atoms
    size_t counter{ 1 }; // for atom numbering in file
    double x0, x1, x2, y0, y1, y2, z0, z1, z2; // coords of 3 water atoms
    std::vector<double> extras{};  // coords of extra sites (xyz order)
    while (cur < na) {
        an = standardise(atom_name(lines[cur]));
        if (an == "OW" && standardise(atom_name(lines[cur+1])) == "HW"
                    && standardise(atom_name(lines[cur+2])) == "HW") {

            // extract coords of OW, HW1, HW2
            coordinates(lines[cur],x0,y0,z0);
            coordinates(lines[cur+1],x1,y1,z1);
            coordinates(lines[cur+2],x2,y2,z2);
        
            // idealise coordinates & store extra sites in extras
            extras = wm.transform(x0,y0,z0,x1,y1,z1,x2,y2,z2);
        
            // convert from Angstrom to nm for gro format
            x0 /= 10.0; y0 /= 10.0; z0 /= 10.0;
            x1 /= 10.0; y1 /= 10.0; z1 /= 10.0;
            x2 /= 10.0; y2 /= 10.0; z2 /= 10.0;
        
            // write updated water atoms
            os << update_line(lines[cur], counter, x0, y0, z0) << '\n';
            os << update_line(lines[cur+1], counter+1, x1, y1, z1) << '\n';
            os << update_line(lines[cur+2], counter+2, x2, y2, z2) << '\n';
            
            // write possible extra sites
            // M site (4-site models)
            if (model_size == 4) {
                os << update_line(lines[cur+2], "MW", counter+3, extras[0]/10.0,
                                  extras[1]/10.0, extras[2]/10.0);
                os << '\n';
            }

            // LP sites (5-site models)
            if (model_size == 5) {
                os << update_line(lines[cur+2], "LP1", counter+3,
                        extras[0]/10.0, extras[1]/10.0, extras[2]/10.0);
                os << '\n';
                os << update_line(lines[cur+2], "LP2", counter+4,
                        extras[3]/10.0, extras[4]/10.0, extras[5]/10.0);
                os << '\n';

            }

            // skip extra sites of original model if present
            cur += 2; //add 2: the while loop will run at least once
            an = "MW";
            while ((cur < na + 2) && (an=="MW" || an=="LP" || an=="EP")) {
                ++cur;
                an = standardise(atom_name(lines[cur]));
            }
            counter += model_size;
        } else {
            // replace atom counter, remove velocities
            os << update_line(lines[cur],counter) << '\n';
            ++cur;
            ++counter;
        }
    } 
    
    // copy the rest of the file to output
    while (cur < lines.size()) {
        if (cur < na + 2) {       // still atoms
            os <<  update_line(lines[cur],counter) << '\n';
            ++cur; ++counter;
        } else {
            os << lines[cur] << '\n';
            ++cur;
        }
    }

    
    
    modified = nw;
    
    
    return modified;
}

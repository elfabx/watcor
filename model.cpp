#include "model.h"
#include <cmath>
#include <stdexcept>
// for debug: #include <iostream>

// needed to assign memory in linker
constexpr model_param const model::models[];

// degree to radian conversion
const double degree = 0.017453292519943295769139;


std::vector<std::string> model::catalog() {
    std::vector<std::string> m;
    m.clear();
    for (int i = 0; i < n_models; ++i) {
        m.push_back(models[i].name);
    }
    return m;
}


bool model::initialise(int id) {
    if (id >= n_models || id < 0) {
        is_initialised = false;
    } else {
        parameters = models[id];
        is_initialised = true;
    }
    return is_initialised;
}


bool model::initialise(const std::string &name) {
    is_initialised = false;
    bool m{ false };
    for (int i = 0; i < n_models; ++i) {
        if (name == models[i].name) { m = initialise(i); }
    }
    return m;
}

void model::check() const {
    if (!is_initialised) {
        throw(std::logic_error("access to uninitialised model"));
    }
}

int model::size() const {
    check();
    int s{ 3 };
    if (std::fabs(parameters.rOM) > 1.0e-4) { ++s; }
    if (std::fabs(parameters.rOL) > 1.0e-4) { s += 2; }
    return s;
}

static double vec_length(double x, double y, double z)
{
    return std::sqrt(x*x + y*y + z*z);
}

std::vector<double> model::transform(double &xO, double &yO, double &zO,
                                  double &x1, double &y1, double &z1,
                                  double &x2, double &y2, double &z2) const {
    check();
    // O-H vectors and their lengths (v1, v2)
    double vx1{ x1 - xO };
    double vy1{ y1 - yO };
    double vz1{ z1 - zO };
    double lv1 { vec_length(vx1,vy1,vz1) };
    double vx2{ x2 - xO }; 
    double vy2{ y2 - yO }; 
    double vz2{ z2 - zO };
    double lv2{ vec_length(vx2, vy2, vz2) };
    
    // if either O-H is too short, bail out (cut-off 0.0001 A)
    if (lv1 < 1.0e-4 || lv2 < 1.0e-4) {
        throw(std::runtime_error("bad input water structure"));
    }
 
    // normalise O-H vectors (v1, v2)
    vx1 /= lv1; vy1 /= lv1; vz1 /= lv1;
    vx2 /= lv2; vy2 /= lv2; vz2 /= lv2;
    
    // bisector direction: sum of unit vectors along O-H bonds
    double ax{ vx1 + vx2 };
    double ay{ vy1 + vy2 };
    double az{ vz1 + vz2 };
    double la{ vec_length(ax, ay, az) }; // length of bisector
    
    // if OH vectors are collinear this is either 0 or 2
    if (la < 1.0e-4 || la > 1.9999) {
        throw(std::runtime_error("bad input water structure"));
    }
    // normalise bisector (this is unit vector a)
    ax /= la; ay /= la; az /= la;
    
    // diff v1-v2 is roughly the H...H direction (perpendicular to a)
    // (this is now unit vector b)
    double bx{ vx1 - vx2 };
    double by{ vy1 - vy2 };
    double bz{ vz1 - vz2 };
    double lb{ vec_length(bx, by, bz) };
    bx /= lb; by /= lb; bz /= lb; // O-Hs not collinear, so length > 0
    
    // new O-H vectors: [ a*cos(angle/2) +/- b*sin(angle/2) ] * length
    double cosa{ std::cos(parameters.angle*degree/2.0)*parameters.rOH };
    double acx{ ax * cosa }; // components along bisector (vec. a)
    double acy{ ay * cosa };
    double acz{ az * cosa };
    double sinb{ std::sin(parameters.angle*degree/2.0)*parameters.rOH };
    double bsx{ bx * sinb }; // components along vector b 
    double bsy{ by * sinb }; 
    double bsz{ bz * sinb };
    
    // add these components to position of O atom
    x1 = xO + acx + bsx;
    x2 = xO + acx - bsx;
    y1 = yO + acy + bsy;
    y2 = yO + acy - bsy;
    z1 = zO + acz + bsz;
    z2 = zO + acz - bsz;
    
    // build up extra sites
    std::vector<double> sites{};
    
    // M site is along the bisector of H-O-H, ie. 'a' vector
    if (std::fabs(parameters.rOM) > 1.0e-4) {
        double xm{ xO + ax * parameters.rOM };
        double ym{ yO + ay * parameters.rOM };
        double zm{ zO + az * parameters.rOM };
        sites.push_back(xm);
        sites.push_back(ym);
        sites.push_back(zm);
    }

    // calculate lone pair sites if rOL is not 
    if (std::fabs(parameters.rOL) > 1.0e-4) {
        // c vector is perp. to water plane (cross product of a & b)
        double cx{ ay * bz - az * by };
        double cy{ az * bx - ax * bz };
        double cz{ ax * by - ay * bx };

        // components of O-Lp vector
        // along a: - cos(lpangle/2)*rOL
        double cosl{ std::cos(parameters.lpangle*degree/2.0)*parameters.rOL };
        // along c: +/- sin(lpangle/2)*rOL
        double sinl{ std::sin(parameters.lpangle*degree/2.0)*parameters.rOL };
        // add xyz components of these vectors to O coordinates -> Lp coords
        double xl1{ xO + cx * sinl - ax * cosl };
        double xl2{ xO - cx * sinl - ax * cosl };
        double yl1{ yO + cy * sinl - ay * cosl };
        double yl2{ yO - cy * sinl - ay * cosl };
        double zl1{ zO + cz * sinl - az * cosl };
        double zl2{ zO - cz * sinl - az * cosl };

        sites.push_back(xl1);
        sites.push_back(yl1);
        sites.push_back(zl1);
        sites.push_back(xl2);
        sites.push_back(yl2);
        sites.push_back(zl2);

    }

    return sites;
    
}

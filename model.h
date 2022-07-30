#ifndef MODEL_H
#define MODEL_H
#include <string>
#include <vector>


//! Geometric parameters of a water model
/**
 * \note rOM and rOL are zero in models without M or LP sites
 */
struct model_param {
    const char *name; //!< model name
    double rOH;       //!< O-H distance
    double angle;     //!< H-O-H angle
    double rOM;       //!< O-M distance (along bisector of HOH angle)
    double rOL;       //!< O-Lone Pair distance
    double lpangle;   //!< Lp-O-Lp angle (in plane perp. to HOH plane)
};

//! Class to set up and perform geometric caclulations using a water model
class model {
public:
    //! Constructor; takes no parameters
    model() : is_initialised{ false } {}

    //! return list of known model names
    static std::vector<std::string> catalog();

    //! set up model params based on model name
    /**
     * Models must be succesfully initialised before any other method
     * can be called, except catalog().
     *
     * \param name descriptive name of the model
     * \return true if set-up is successful
     */
    bool initialise(const std::string &name);

    //! set up model params bases on sequence number of model
    /**
     * Models must be succesfully initialised before any other method
     * can be called, except catalog().
     *
     * \param id sequence number of model (see catalog)
     * @return true if set-up is successful
     */
    bool initialise(int id);

    int size() const;  //!< gives the number of sites in current model
    
    //! change water coordinates to idealised model geometry
    /**
     * O remains in place, H atoms are placed at equal angle from the original
     * bisector of the H-O-H angle
     *
     * \param[in,out] xO,yO,zO coordinates of the O atom
     * \param[in,out] x1,y1,z1,x2,y2,z2 coordinates of the two H atoms
     *
     * \returns a vector of either x, y, z coordinates for a single M site
     *     or x1, y1, z1, x2, y2, z2 for the two Lp sites (both are optional)
     *
     * \note The algorithm in gro.cpp does not work with both M and Lp sites
     *     althogh this routine would support it (I know of no such model)
        */
    std::vector<double> transform(double &xO, double &yO, double &zO,  // Ow
                           double &x1, double &y1, double &z1,         // Hw1
                           double &x2, double &y2, double &z2) const;  // Hw2
protected:
    model_param parameters; //!< a copy of the current model parameters
    bool is_initialised;    //!< flag to show the model is initialised
    void check() const;    //!< throw a logic_error if model is not initialised

    static constexpr int n_models{ 12 }; //!< number of known models

    //! contains all parameters for the known models
    static constexpr model_param const models[] = {
        {"tip3p",    0.9572, 104.52, 0.00,   0.0, 0.0 },
        {"tip3p-fb", 1.0118, 108.15, 0.00,   0.0, 0.0 },
        {"spc/e",    1.0000, 109.47, 0.00,   0.0, 0.0 },
        {"spc/fw",   1.0120, 113.24, 0.00,   0.0, 0.0 },
        {"spc/eb",   1.0100, 109.47, 0.00,   0.0, 0.0 },
        {"opc3",    0.97888, 109.47, 0.00,   0.0, 0.0 },
        {"opc",     0.87243, 103.60, 0.1594, 0.0, 0.0 },
        {"tip4p",    0.9572, 104.52, 0.15,   0.0, 0.0 },
        {"tip4p-ew", 0.9572, 104.52, 0.125,  0.0, 0.0 },
        {"tip4p-fb", 0.9572, 104.52, 0.10527,0.0, 0.0 },
        {"tip5p",    0.9572, 104.52, 0.00,   0.7, 109.47},
        {"tip5p-e",  0.9572, 104.52, 0.00,   0.7, 109.47}
    };
};


#endif

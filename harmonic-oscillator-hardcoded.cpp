#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>

//#############################################################################################
//################################### Define Constants ########################################
//#############################################################################################

const double m1 = 1.67e-27;  // Mass of the first atom (kg)
const double m2 = 1.67e-27;  // Mass of the second atom (kg)
const double k = 500;        // Spring constant (N/m)
const double dt = 1e-16;     // Time step (s)
const double t_end = 1e-12;  // End time (s)

//#############################################################################################
//#################################### Define Classes #########################################
//#############################################################################################


//##############################################################################################
//#################################### Define Functions ########################################
//##############################################################################################


class Atom {
private:
    std::string atomName;
    std::string atomType;
    double x, y, z;

public:
    Atom(const std::string& name, const std::string& type) : atomName(name), atomType(type) {}

    // Getter methods
    std::string getAtomName() const { return atomName; }
    std::string getAtomType() const { return atomType; }
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }

    // Setter methods
    void setAtomName(const std::string& name) { atomName = name; }
    void setAtomType(const std::string& type) { atomType = type; }
    void setCoordinates(double x_coord, double y_coord, double z_coord) {
        x = x_coord;
        y = y_coord;
        z = z_coord;
    }
};

std::vector<Atom> readpdb(char *filename)
{	

	std::vector<Atom> atoms;
	std::ifstream pdb(filename);
	std::string line;


	if (!file.is_open()) 
    {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return ""; // Return an empty string if file opening fails
    }

    while (std::getline(file, line))
    {
    	if (line.substr(0, 6) == "ATOM")
    	{
    		std::istringstream iss(line.substr(30, 55));
    		double x, y, z;
    		iss >> x >> y >> z;
    		std::string name = line.substr(12, 4);
    		std::string type = line.substr(76, 2);
    		Atom atom(name, type);
    		atom.setCoordinates(x, y, z);
    		atoms.push_back(atom);
    	}
    } 

	return atoms;
}


int main(int argc, char const *argv[])
{

	char *pdb;
    char *nstep;
    char *out;

	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-pdb") == 0)
		{
			pdb = argv[i + 1];
		}

        if (strcmp(argv[i], "-nstep") == 0)
        {
            forcefield = argv[i + 1];
        }

        if (strcmp(argv[i], "-out") == 0)
        {
            forcefield = argv[i + 1];
        }

	}

	return 0;
}
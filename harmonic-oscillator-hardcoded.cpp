#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstring>

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


//##############################################################################################
//#################################### Define Functions ########################################
//##############################################################################################

std::vector<Atom> readpdb(std::string filename)
{	

	std::vector<Atom> atoms;
	std::ifstream file(filename);
	std::string line;


	if (!file.is_open()) 
    {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return std::vector<Atom> ();
    }

    while (std::getline(file, line))
    {	

    	if (line.substr(0, 4) == "ATOM")
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

void writePDBModelToTrajectory(const std::string& filename, const std::vector<Atom>& atoms, int modelNumber) 
{
    std::ofstream outFile(filename, std::ios_base::app); // Open the file in append mode
    
    if (!outFile.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write MODEL record
    outFile << "MODEL     " << std::setw(4) << std::setfill(' ') << modelNumber << std::endl;

    // Write atom records
    for (const auto& atom : atoms) 
    {
        outFile << std::fixed << std::setprecision(3); // Set precision for coordinates
        outFile << "ATOM  " << std::setw(5) << std::setfill(' ') << std::right << 1 << " "; // Atom serial number
        outFile << std::left << std::setw(4) << atom.getAtomName(); // Atom name
        outFile << "   " << atom.getAtomType() << " "; // Atom type
        outFile << " " << " " << std::setw(4) << std::right << 1 << " "; // Residue sequence number
        outFile << "    " << std::setw(8) << std::right << atom.getX(); // X coordinate
        outFile << std::setw(8) << std::right << atom.getY(); // Y coordinate
        outFile << std::setw(8) << std::right << atom.getZ(); // Z coordinate
        outFile << "  1.00  0.00           " << std::endl; // Default occupancy and temperature factor
    }

    // Write ENDMDL record
    outFile << "ENDMDL" << std::endl;

    outFile.close();
}


std::array<double, 3> calculate_force(const std::array<double, 3>& pos1, const std::array<double, 3>& pos2)
{
    std::array<double, 3> force;

    for (int i = 0; i < 3; ++i) 
    {
        force[i] = -k * (pos1[i] - pos2[i]);
    }
    return force;
}

// Runge-Kutta 4th-order integration step
void runge_kutta_step(std::array<double, 3>& pos1, std::array<double, 3>& vel1, std::array<double, 3>& pos2, std::array<double, 3>& vel2) 
{    
    std::array<double, 3> k1_v1, k1_v2, k1_x1, k1_x2;
    std::array<double, 3> k2_v1, k2_v2, k2_x1, k2_x2;
    std::array<double, 3> k3_v1, k3_v2, k3_x1, k3_x2;
    std::array<double, 3> k4_v1, k4_v2, k4_x1, k4_x2;

    // Compute k1
    auto force = calculate_force(pos1, pos2);
    for (int i = 0; i < 3; ++i) 
    {
        k1_v1[i] = dt * force[i] / m1;
        k1_v2[i] = -dt * force[i] / m2;
        k1_x1[i] = dt * vel1[i];
        k1_x2[i] = dt * vel2[i];
    }

    // Compute k2
    std::array<double, 3> pos1_temp = pos1, pos2_temp = pos2;
    std::array<double, 3> vel1_temp = vel1, vel2_temp = vel2;
    for (int i = 0; i < 3; ++i) 
    {
        pos1_temp[i] += k1_x1[i] / 2;
        pos2_temp[i] += k1_x2[i] / 2;
        vel1_temp[i] += k1_v1[i] / 2;
        vel2_temp[i] += k1_v2[i] / 2;
    }
    force = calculate_force(pos1_temp, pos2_temp);
    for (int i = 0; i < 3; ++i) 
    {
        k2_v1[i] = dt * force[i] / m1;
        k2_v2[i] = -dt * force[i] / m2;
        k2_x1[i] = dt * vel1_temp[i];
        k2_x2[i] = dt * vel2_temp[i];
    }

    // Compute k3
    pos1_temp = pos1;
    pos2_temp = pos2;
    vel1_temp = vel1;
    vel2_temp = vel2;
    for (int i = 0; i < 3; ++i) 
    {
        pos1_temp[i] += k2_x1[i] / 2;
        pos2_temp[i] += k2_x2[i] / 2;
        vel1_temp[i] += k2_v1[i] / 2;
        vel2_temp[i] += k2_v2[i] / 2;
    }
    force = calculate_force(pos1_temp, pos2_temp);
    for (int i = 0; i < 3; ++i) 
    {
        k3_v1[i] = dt * force[i] / m1;
        k3_v2[i] = -dt * force[i] / m2;
        k3_x1[i] = dt * vel1_temp[i];
        k3_x2[i] = dt * vel2_temp[i];
    }

    // Compute k4
    pos1_temp = pos1;
    pos2_temp = pos2;
    vel1_temp = vel1;
    vel2_temp = vel2;
    for (int i = 0; i < 3; ++i) {
        pos1_temp[i] += k3_x1[i];
        pos2_temp[i] += k3_x2[i];
        vel1_temp[i] += k3_v1[i];
        vel2_temp[i] += k3_v2[i];
    }
    force = calculate_force(pos1_temp, pos2_temp);
    for (int i = 0; ++i < 3; ++i) 
    {
        k4_v1[i] = dt * force[i] / m1;
        k4_v2[i] = -dt * force[i] / m2;
        k4_x1[i] = dt * vel1_temp[i];
        k4_x2[i] = dt * vel2_temp[i];
    }

    // Update positions and velocities
    for (int i = 0; i < 3; ++i) 
    {
        pos1[i] += (k1_x1[i] + 2 * k2_x1[i] + 2 * k3_x1[i] + k4_x1[i]) / 6;
        pos2[i] += (k1_x2[i] + 2 * k2_x2[i] + 2 * k3_x2[i] + k4_x2[i]) / 6;
        vel1[i] += (k1_v1[i] + 2 * k2_v1[i] + 2 * k3_v1[i] + k4_v1[i]) / 6;
        vel2[i] += (k1_v2[i] + 2 * k2_v2[i] + 2 * k3_v2[i] + k4_v2[i]) / 6;
    }
}

//##############################################################################################
//#################################### Main Function ###########################################
//##############################################################################################


int main(int argc, char *argv[])
{

	char *pdbName;
    char *outFileName;

	for (int i = 0; i < argc; i++)
	{
        if (std::strcmp(argv[i], "-pdb") == 0 && i + 1 < argc)
		{
			pdbName = argv[i + 1];
		}

        if (std::strcmp(argv[i], "-out") == 0 && i + 1 < argc)
        {
            outFileName = argv[i + 1];
        }
	}

   
    std::string pdb = std::string(pdbName);
    std::string outFile = std::string(outFileName);

     // Check if required arguments are provided
    if (pdb.empty() || outFile.empty()) 
    {
        std::cerr << "Usage: " << argv[0] << " -pdb <pdb_file> -out <output_file>" << std::endl;
        return 1;
    }

	std::vector <Atom> atoms = readpdb(pdb);

	if (atoms.size() < 1) 
	{
    	std::cerr << "Error: Insufficient atoms in the vector." << std::endl;
    	// Handle the error appropriately, e.g., return from the function or exit the program
	}

	std::array<double, 3> x1 = {atoms[0].getX(), atoms[0].getY(), atoms[0].getZ()}; // Initial position of atom 1 (m)
    std::array<double, 3> v1 = {0.0, 0.0, 0.0}; // Initial velocity of atom 1 (m/s)
    std::array<double, 3> x2 = {atoms[1].getX(), atoms[1].getY(), atoms[1].getZ()}; // Initial position of atom 2 (m)
    std::array<double, 3> v2 = {0.0, 0.0, 0.0};

    double time = 0.0;
    int modelNumber = 0;

    while (time < t_end)
    {	

    	writePDBModelToTrajectory(outFile, atoms, modelNumber);

    	runge_kutta_step(x1, v1, x2, v2);
    	
    	atoms[0].setCoordinates(x1[0], x1[1], x1[2]);
    	atoms[1].setCoordinates(x2[0], x2[1], x2[2]);
    
    	modelNumber++;
    	time += dt;
    }

	return 0;
}
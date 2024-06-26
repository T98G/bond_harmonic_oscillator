#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <regex>
#include <string>
#include <cctype>


//#############################################################################################
//################################### Define Constants ########################################
//#############################################################################################

//const double m1 = 1.67e-27;  // Mass of the first atom (kg)
//const double m2 = 1.67e-27;  // Mass of the second atom (kg)
//const double k = 500;        // Spring constant (N/m)
const double dt = 1e-16;     // Time step (s)
const double t_end = 1e-13;  // End time (s)

//#############################################################################################
//#################################### Define Class  ##########################################
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

//Define functions to search the forcefield

std::string extractString(const std::string& input)
{
    size_t start = input.find('"');
    if (start == std::string::npos) // If no opening quotation mark is found
        return "";
    size_t end = input.find('"', start + 1);
    if (end == std::string::npos) // If no closing quotation mark is found
        return "";
    return input.substr(start + 1, end - start - 1);
}

std::string readff(const std::string& filename, std::vector<std::string>& visitedFiles) 
{
    // Check if the file has already been visited to prevent infinite loops
    if (std::find(visitedFiles.begin(), visitedFiles.end(), filename) != visitedFiles.end()) 
    {
        // File has already been visited, return an empty string
        return "";
    }
    // Add the current file to the visited files list
    visitedFiles.push_back(filename);

    // Open the file
    std::ifstream file(filename);

    // Check if the file is opened successfully
    if (!file.is_open()) 
    {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return ""; // Return an empty string if file opening fails
    }

    std::regex regexPattern("#include");
    std::string line = "";
    std::string nextFile = "";

    // Find the last slash in the path to extract the directory
    size_t lastSlash = filename.find_last_of('/');
    std::string directory = filename.substr(0, lastSlash + 1);

    // Resulting string to store the concatenation of all included files
    std::string result = "";

    // Read and check each line until the end of the file
    while (std::getline(file, line)) 
    {
        if (std::regex_search(line, regexPattern)) 
        {
            nextFile = extractString(line); 

            std::cout << "Including file: " << nextFile << std::endl;

            // Construct the full path for the next file
            std::string nextFilePath = directory + nextFile;

            // Recursively read the next file and append its content to the result
            result += readff(nextFilePath, visitedFiles);
        }
        else 
        {
            result += line + "\n"; // Append the line to the result
        }
    }

    // Close the file
    file.close();

    return result;
}

std::string extractAtomType(const std::string& ffContents, const std::string& atomName) 
{
    std::istringstream iss(ffContents);
    std::string line;
    std::string atomType;

    // Construct the regex pattern for finding the atom type
    std::regex atomPattern("^\\s*\\d+\\s+(\\w+)\\s+\\d+\\s+\\w+\\s+" + atomName + "\\s+");

    // Search for the atom type based on the pattern
    while (std::getline(iss, line)) 
    {
        std::smatch match;
        if (std::regex_search(line, match, atomPattern)) 
        {
            std::cout << "Match found: " << match[0] << std::endl; // Print the entire match
            atomType = match[1]; // Extract the captured group containing the atom type
            break;  // Stop searching after finding the atom type
        }
    }

    return atomType;
}

double findAtomsinFF(const std::string& ffContents, const std::string& atomType) {
    std::istringstream iss(ffContents);
    std::string line;

    // Adjusted regex to match atom type up to 4 characters followed by a floating-point number
    std::regex regexAtom(atomType);

    double mass = 0;

    // Read and check each line until the end of the file
    while (std::getline(iss, line)) {
        std::smatch match;

        if (std::regex_search(line, match, regexAtom)) {

            std::cout << line << std::endl;

            // Find the position of the atom type in the line
            size_t pos = line.find(atomType);
            if (pos != std::string::npos) {
                // Extract the number after the atom type
                std::string numberStr = line.substr(pos + atomType.length());
                try {
                    mass = std::stod(numberStr);
                    std::cout << mass << std::endl;
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Conversion error: " << e.what() << std::endl;
                }
            }
            break;  // Assuming we only need the first match
        }
    }

    return mass;
}

double findBondinFF(const std::string& ffContents, const std::string& atomType1, const std::string& atomType2) 
{
    std::istringstream iss(ffContents);
    std::string line;

    std::string pattern1 = atomType1 + "\\s+" + atomType2 + "\\s+\\d+\\s+\\d+\\.\\d+\\s+(\\d+\\.\\d+)";
    std::string pattern2 = atomType2 + "\\s+" + atomType1 + "\\s+\\d+\\s+\\d+\\.\\d+\\s+(\\d+\\.\\d+)";
    std::regex regexPattern1(pattern1);
    std::regex regexPattern2(pattern2);

    bool inBondTypesSection = false;

    // Read and check each line until the end of the file
    while (std::getline(iss, line)) 
    {
        std::smatch match;

        if (line.find("[ bondtypes ]") != std::string::npos) 
        {
            inBondTypesSection = true;
            continue;
        }
        if (inBondTypesSection && line.find("[") != std::string::npos) 
        {
            break; // Exit loop when the next section is encountered
        }
        if (inBondTypesSection) 
        {
            if (std::regex_search(line, match, regexPattern1) || std::regex_search(line, match, regexPattern2)) 
            {
                return std::stod(match[1].str()); // Return the captured spring constant
            }
        }
    }

    return 0.0; // Return 0.0 if the pattern is not found
}


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

std::array<double, 3> calculate_force(const std::array<double, 3>& pos1, const std::array<double, 3>& pos2, double k)
{
    std::array<double, 3> force;

    for (int i = 0; i < 3; ++i) 
    {
        force[i] = -k * (pos1[i] - pos2[i]);
    }
    return force;
}

// Runge-Kutta 4th-order integration step
void runge_kutta_step(std::array<double, 3>& pos1, std::array<double, 3>& vel1, std::array<double, 3>& pos2, std::array<double, 3>& vel2, double m1, double m2, double k) 
{    
    std::array<double, 3> k1_v1, k1_v2, k1_x1, k1_x2;
    std::array<double, 3> k2_v1, k2_v2, k2_x1, k2_x2;
    std::array<double, 3> k3_v1, k3_v2, k3_x1, k3_x2;
    std::array<double, 3> k4_v1, k4_v2, k4_x1, k4_x2;

    // Compute k1
    auto force = calculate_force(pos1, pos2, k);
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

    force = calculate_force(pos1_temp, pos2_temp, k);
    
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
    
    force = calculate_force(pos1_temp, pos2_temp, k);
    
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
    
    for (int i = 0; i < 3; ++i)
    {
        pos1_temp[i] += k3_x1[i];
        pos2_temp[i] += k3_x2[i];
        vel1_temp[i] += k3_v1[i];
        vel2_temp[i] += k3_v2[i];
    }
   
    force = calculate_force(pos1_temp, pos2_temp, k);
    
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
    // Declare pointers for file names and initialize them to nullptr
    char *pdbName = nullptr;
    char *outFileName = nullptr;
    char *top = nullptr;

    // Parse the command-line arguments
    for (int i = 0; i < argc; i++)
    {
        if (std::strcmp(argv[i], "-pdb") == 0 && i + 1 < argc)
        {
            pdbName = argv[i + 1];
        }
        else if (std::strcmp(argv[i], "-out") == 0 && i + 1 < argc)
        {
            outFileName = argv[i + 1];
        }
        else if (std::strcmp(argv[i], "-top") == 0 && i + 1 < argc)
        {
            top = argv[i + 1];
        }
    }

    // Check if required arguments are provided
    if (pdbName == nullptr || outFileName == nullptr || top == nullptr)
    {
        std::cerr << "Usage: " << argv[0] << " -pdb <pdb_file> -top <topology file> -out <output_file>" << std::endl;
        return 1;
    }

    std::string topol = std::string(top);
    std::vector<std::string> visitedFiles;
    std::string ffContents = readff(topol, visitedFiles);

    std::string pdb = std::string(pdbName);
    std::string outFile = std::string(outFileName);


    std::vector <Atom> atoms = readpdb(pdb);

    // Atom names to extract the atom types
    std::string atomName1 = atoms[0].getAtomName();
    std::string atomName2 = atoms[1].getAtomName();

    // Extract atom types
    std::string atomType1 = extractAtomType(ffContents, atomName1);
    std::string atomType2 = extractAtomType(ffContents, atomName2);

    double k = findBondinFF(ffContents, atomType1, atomType2);

    double m1 = findAtomsinFF(ffContents, atomType1) * 1e-27;
    double m2 = findAtomsinFF(ffContents, atomType2) * 1e-27;

    //std::cout << m1 << std::endl;
    //std::cout << m2 << std::endl;

    if (atoms.size() < 2)
    {
        std::cerr << "Error: Insufficient atoms in the vector." << std::endl;
        return 1;
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

        runge_kutta_step(x1, v1, x2, v2, m1, m2, k);

        atoms[0].setCoordinates(x1[0], x1[1], x1[2]);
        atoms[1].setCoordinates(x2[0], x2[1], x2[2]);

        modelNumber++;
        time += dt;
    }

    return 0;
}
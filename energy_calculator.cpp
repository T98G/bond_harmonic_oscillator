#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <regex>

class Atom
{
private:
	std::string atomName;
	std::string atomType;

	//Add coordinates for atoms

public:
	Atom(const std::string& name, const std::string& type) : atomName(name), atomType(type) {}

	// Getter methods
	std::string getAtomName() const { return atomName; }
	std::string getAtomType() const { return atomType; }

	// Setter methods
    void setAtomName(const std::string& name) { atomName = name; }
    void setAtomType(const std::string& type) { atomType = type; }

};

class Bond
{
private:

	Atom atom1;
	Atom atom2;
	double eq_distance;

public:
	Bond(const Atom& atm1, const Atom& atm2, double eq_dist) : atom1(atm1), atom2(atm2), eq_distance(eq_dist) {}

};

std::string findBondinFF(const std::string& filename, const std::string& variable1, const std::string& variable2) {
    // Open the file
    std::ifstream file(filename);

    // Check if the file is opened successfully

    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return ""; // Return an empty string if file opening fails
    }

    // Construct the regular expression pattern with variable spacing
    std::string pattern1 = variable1 + "\\s+" + variable2; // \\s+ matches one or more whitespace characters
    std::string pattern2 = variable2 + "\\s+" + variable1;
    std::regex regexPattern1(pattern1);
    std::regex regexPattern2(pattern2);

    std::string line;
    bool inBondTypesSection = false;

    // Read and check each line until the end of the file
    while (std::getline(file, line)) {
        // Check if the line matches the regular expression pattern

    	// Check if the line contains the start of the bond types section
        if (line.find("[ bondtypes ]") != std::string::npos) {
            inBondTypesSection = true;
            continue; // Skip this line and move to the next
        }

        // Check if the line contains the start of the next section
        if (inBondTypesSection && line.find("[") != std::string::npos) {
            // We've reached the end of the bond types section
            break; // Exit the loop
        }

        if (std::regex_search(line, regexPattern1) || std::regex_search(line, regexPattern1)) {
            // Close the file
            file.close();
            return line; // Return the matched line
        }
    }

    // Close the file
    file.close();
    return ""; // Return an empty string if the pattern is not found
}



int main(int argc, char *argv[])
{	

	char *pdb;
	char *topol;
	char *forcefield;
	char *mol;
	char *atom1;
	char *atom2;

	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-pdb") == 0)
		{
			pdb = argv[i + 1];
		}

		if (strcmp(argv[i], "-top") == 0)
		{
			topol = argv[i + 1];
		}

		if (strcmp(argv[i], "-ff") == 0)
		{
			forcefield = argv[i + 1];
		}

		if (strcmp(argv[i], "-mol") == 0)
		{
			mol = argv[i + 1];
		}

		if (strcmp(argv[i], "-atom1") == 0)
		{
			atom1 = argv[i + 1];
		}

		if (strcmp(argv[i], "-atom2") == 0)
		{
			atom2 = argv[i + 1];
		}

	}

	std::string atomType1 = "CA";
	std::string atomType2 = "OH";

	std::string line = findBondinFF(forcefield, atomType1, atomType2);

	std::cout << line << std::endl; 


	return 0;
}
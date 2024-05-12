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
	Bond(const Atom& atm1, const Atom& atm2, double eq_dist) : atom1(atm1), atom2(atm2), eq_distance(eq_len) {}

};

std::string readff(char fileName, std::string atomType1, std::string atomType2)
{
	std::ifstream file(fileName);

	if (!file.is_open()) 
	{
    	std::cerr << "Failed to open Force Field file." << std::endl;
        return std::string ""; // return error code
    }

    std::string line;
    std::regex regexPattern(atomType1 + "\\s" + atomType2);
    std::regex regexPattern(atomType2 + "\\s" + atomType1);

    // Read and check each line until the end of the file
    
    while (std::getline(file, line))
    {
    	if (std::regex_search(line, regexPattern1) || std::regex_search(line, regexPattern2))
    	{
    		return line;
    	}
    }

    return std::string "";
}

std::string readtop(char fileName, char *atomName1, char *atomName2)
{
	
}



int main(int argc, char *argv[])
{	

	char *pdb;
	char *topol;
	char *forcefield;
	char *mol;
	char *atom1;
	char *atom2;

	for (int i = 0; i++; i < argc)
	{
		if (argv[i] == "-pdb")
		{
			pdb = argv[i + 1];
		}

		if (argv[i] == "-top")
		{
			topol = argv[i + 1];
		}

		if (argv[i] == "-ff")
		{
			forcefield = argv[i + 1];
		}

		if (argv[i] == "-mol")
		{
			mol = argv[i + 1];
		}

		if (argv[i] == "-atom1")
		{
			atom1 = argv[i + 1];
		}

		if (argv[i] == "-atom2")
		{
			atom2 = argv[i + 1];
		}

	}


	return 0;
}
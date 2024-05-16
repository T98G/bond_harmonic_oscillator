#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <regex>
#include <vector>


//##########################################
//########## Define Classes ################
//##########################################


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

//##########################################
//########## Define Functions ##############
//##########################################

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


std::string findBondinFF(const std::string ffContents, const std::string& atom1, const std::string& atom2)
{


    // Read the finalResult line by line
    std::istringstream iss(ffContents);
    std::string line;

    std::string atomType1;
    std::string atomType2;

    while(std::getline(iss, line))
    {

        if (line.find("[ atoms ]") != std::string::npos) 
        {
            inBondTypesSection = true;
            continue; // Skip this line and move to the next
        }

        if if (inBondTypesSection && line.find("[") != std::string::npos)
        {
            break;
        }

        if (line.find(atom1) != std::string::npos)
        {
            atom1 = line.substr(13, 4);
        }

        if (line.find(atom2) != std::string::npos)
        {
            atom2 = line.substr(13, 4);
        }
    }


    std::string pattern1 = atomType1 + "\\s+" + atomType2; // \\s+ matches one or more whitespace characters
    std::string pattern2 = atomType2 + "\\s+" + atomType1;
    std::regex regexPattern1(pattern1);
    std::regex regexPattern2(pattern2);

    std::string line;
    bool inBondTypesSection = false;

    // Read and check each line until the end of the file
    while (std::getline(iss, line)) 
    {
        // Check if the line contains the start of the bond types section
        if (line.find("[ bondtypes ]") != std::string::npos) 
        {
            inBondTypesSection = true;
            continue; // Skip this line and move to the next
        }

        // Check if the line contains the start of the next section
        if (inBondTypesSection && line.find("[") != std::string::npos) 
        {
            // We've reached the end of the bond types section
            break; // Exit the loop
        }

        if (std::regex_search(line, regexPattern1) || std::regex_search(line, regexPattern1)) 
        {
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



	return 0;
}
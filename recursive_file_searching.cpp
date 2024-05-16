#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <regex>

std::string extractString(const std::string& input) {
    size_t start = input.find('"');
    if (start == std::string::npos) // If no opening quotation mark is found
        return "";
    size_t end = input.find('"', start + 1);
    if (end == std::string::npos) // If no closing quotation mark is found
        return "";
    return input.substr(start + 1, end - start - 1);
}

std::string readff(const std::string& filename, std::vector<std::string>& visitedFiles) {
    // Check if the file has already been visited to prevent infinite loops
    if (std::find(visitedFiles.begin(), visitedFiles.end(), filename) != visitedFiles.end()) {
        // File has already been visited, return an empty string
        return "";
    }
    // Add the current file to the visited files list
    visitedFiles.push_back(filename);

    // Open the file
    std::ifstream file(filename);

    // Check if the file is opened successfully
    if (!file.is_open()) {
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
    while (std::getline(file, line)) {
        if (std::regex_search(line, regexPattern)) {
            nextFile = extractString(line); 

            std::cout << "Including file: " << nextFile << std::endl;

            // Construct the full path for the next file
            std::string nextFilePath = directory + nextFile;

            // Recursively read the next file and append its content to the result
            result += readff(nextFilePath, visitedFiles);
        }
        else {
            result += line + "\n"; // Append the line to the result
        }
    }

    // Close the file
    file.close();

    return result;
}

std::string findBondinFF(const std::string& ffContents, const std::string& atom1, const std::string& atom2) {
    
    std::istringstream iss(ffContents);
    std::string line;
    std::string atomType1;
    std::string atomType2;

    // Read each line of ffContents

    bool inAtomsSection = false;
    bool inBondtypesSection = false;


    while(std::getline(iss, line))
    {
        if (line.find("[ atoms ]") != std::string::npos) 
        {
            inBondtypesSection = true;
            continue;
        }

        if (inBondtypesSection && line.find("[") != std::string::npos)
            break; // Exit loop when next section is encountered
        // Extract atom types from the line
        if (line.find(atom1) != std::string::npos)
            atomType1 = line.substr(13, 4); // Example: Extracting from columns 14-17
        if (line.find(atom2) != std::string::npos)
            atomType2 = line.substr(13, 4); // Example: Extracting from columns 14-17
    
        
    }


    std::string pattern1 = atomType1 + "\\s+" + atomType2;
    std::string pattern2 = atomType2 + "\\s+" + atomType1;
    std::regex regexPattern1(pattern1);
    std::regex regexPattern2(pattern2);


    iss.clear(); 
    iss.seekg(0, std::ios::beg); 

    while (std::getline(iss, line)) 
    {
        if (line.find("[ bondtypes ]") != std::string::npos) 
        {   
            inAtomsSection = true;
            continue;  
        }

        if (inAtomsSection && line.find("[") != std::string::npos)
            break; // Exit loop when next section is encountered
        
                        
            // Check if current line matches either pattern
        if (std::regex_search(line, regexPattern1) || std::regex_search(line, regexPattern2))
            return line; // Return the matched line    
    }

    return ""; // Return an empty string if the pattern is not found
}


std::string extractAtomType(const std::string& ffContents, const std::string& atomName) {
    std::istringstream iss(ffContents);
    std::string line;
    std::string atomType;

    // Construct the regex pattern for finding the atom type
    std::regex atomPattern(atomName);

    // Search for the atom type based on the pattern
    while (std::getline(iss, line)) {

        if (std::regex_search(line, atomPattern)) {

            atomType = line; // Extract the captured group containing the atom type
            break;  // Stop searching after finding the atom type
        }
    }

    return atomType;
}



int main() {
    // Sample force field contents

    std::vector<std::string> visitedFiles;

    std::string ffContents = readff("test.top", visitedFiles);

    std::cout << ffContents << std::endl;

    // Atom names to extract the atom types
    std::string atomName1 = "O12";
    std::string atomName2 = "P";

    // Extract atom types
    std::string atomType1 = extractAtomType(ffContents, atomName1);
    std::string atomType2 = extractAtomType(ffContents, atomName2);

    // Output the extracted atom types
    std::cout << "Atom Type for " << atomName1 << ": " << atomType1 << std::endl;
    std::cout << "Atom Type for " << atomName2 << ": " << atomType2 << std::endl;

    return 0;
}





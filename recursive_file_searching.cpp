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


int main()
{	
	std::vector<std::string> visitedFiles;

	std::cout << readff("test.top", visitedFiles) << std::endl;


	return 0;
}
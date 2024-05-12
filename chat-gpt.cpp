#include <iostream>
#include <fstream>
#include <string>
#include <regex>

bool isPatternMatch(const std::string& line, const std::string& variable1, const std::string& variable2) {
    // Construct the regular expression pattern with variable spacing
    std::string pattern = variable1 + "\\s+" + variable2; // \\s+ matches one or more whitespace characters
    
    // Compile the regular expression pattern
    std::regex regexPattern(pattern);
    
    // Check if the line matches the regular expression pattern
    return std::regex_search(line, regexPattern);
}

int main() {
    // Open the file
    std::ifstream file("forcefield.ff/ffbonded.itp");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1; // return error code
    }

    std::string variable1 = "CA";
    std::string variable2 = "OH";
    std::string line;
    // Read and check each line until the end of the file
    while (std::getline(file, line)) {
        // Check if the line matches the pattern with variables
        if (isPatternMatch(line, variable1, variable2)) {
            std::cout << "Pattern found: " << line << std::endl;
            // You can perform further operations with the matching line here
        }
    }

    // Close the file
    file.close();

    return 0;
}

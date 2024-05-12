#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

class Atom
{
private:
	std::string atomName;
	std::string atomType;

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
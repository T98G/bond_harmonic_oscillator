# Harmonic Bond Oscillator

## Overview

This program simulates a harmonic bond oscillator, which models the interaction between two atoms connected by a spring. The simulation considers the forces, masses, and spring constants to compute the motion of the atoms over time.

## Features

	•	Simulates the harmonic oscillation of two atoms.
	•	Allows setting atomic properties and coordinates.
	•	Configurable time step and simulation duration.

## File Structure

	•	harmonic_bond_oscillator.cpp: The main source file containing the implementation of the harmonic bond oscillator simulation.

## Getting Started

### Prerequisites

To compile and run this program, you need:

	•	A C++ compiler (e.g., g++).
	•	Standard C++ libraries.

Compilation

To compile the program, use the following command in your terminal:

* g++ -o harmonic_bond_oscillator harmonic_bond_oscillator.cpp

## Running the Program

After compiling, run the executable:

* ./harmonic_bond_oscillator

## Program Structure

### Constants

	•	dt: Time step for the simulation (in seconds).
	•	t_end: End time for the simulation (in seconds).

Class: Atom

Represents an atom with the following properties:

	•	atomName: Name of the atom.
	•	atomType: Type of the atom.
	•	x, y, z: Coordinates of the atom in 3D space.

Methods

	•	Getters: getAtomName(), getAtomType(), getX(), getY(), getZ()
	•	Setters: setAtomName(), setAtomType(), setCoordinates(double x_coord, double y_coord, double z_coord)

## Usage

Modify the constants and parameters in the source file as needed to simulate different scenarios. The results of the simulation, including the positions and velocities of the atoms, will be output to the console or a file as implemented in the code.

### Example

Here is an example of how to define an atom and set its coordinates:

Atom hydrogen("Hydrogen", "H");
hydrogen.setCoordinates(0.0, 0.0, 0.0);

## License

This project is licensed under the MIT License.

## Acknowledgements

This program utilizes basic principles of classical mechanics to model atomic interactions and is intended for educational purposes.

# FirstChallengePacs

## Laplace equation solution in C++ with MPI and OMP

This project is a C++ application that solves the Laplace equation using both sequential and parallel methods. The application reads input data from JSON files and writes the results to VTK files.

## Features
- Supports fast solvability for the Laplace operator in the domain [0,1]x[0,1].  
- The solution can be computed sequentially, or in parallel with MPI, OMP or both, based on the user's choice.  
- Input parameters are read from a JSON file for flexibility and ease of use.  
- Provides options for Dirichlet boundary conditions, both homogeneous and non-homogenous.  
- Possibility to save the solution in a .vtk file, that can be opened from paraview.  


## Usage

1. Clone the repository:

```
git clone git@github.com:Palpal16/FirstChallengePacs.git
```

2. Modify the variable PACS_ROOT in the Makefile, to your repository with muparser library

3. Modify the file data.json with parameters of your choice. The names should be self-explanatory to the parameters meaning.  

4. Compile and run the program

```
make
./main
```

5. To run the code in parallel use:
```
mpirun -np 2 ./main
```
Where the number after -np is the number of processors used in the MPI process.


## Directory Structure and Code Explanation
[include/](/include/) and [src/](src/) Contain rispectively all the header and source files for the project.
[mesh/](mesh/): Contains the VTK files for the output of the mesh.
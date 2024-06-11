# ThirdChallengePacs

## Laplace equation solution in C++ with MPI and OMP

This project is a C++ application that solves the Laplace equation using both sequential and parallel methods. The application reads input data from JSON files and writes the results to VTK files.

## Features
- Supports fast solvability for the Laplace operator in the domain [0,1]x[0,1].  
- The solution can be computed sequentially, or in parallel with MPI, OMP or both, based on the user's choice.  
- Input parameters are read from a JSON file for flexibility and ease of use.  
- Provides options for Dirichlet boundary conditions, both homogeneous and non-homogenous.  
- Possibility to save the solution in a .vtk file, that can be opened from Paraview.  
- Performance test through a .sh file over a changing number of cores  


## Usage

1. Clone the repository:

```
git clone git@github.com:Palpal16/ThirdChallengePacs.git
```

2. Modify the variable PACS_ROOT in the Makefile, to your repository with muparser library

3. Modify the file input/data.json with parameters of your choice. The names should be self-explanatory to the parameters meaning.  

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

6. To open and view the solution on Paraview 
```
paraview output/myMesh.vtk
```

7. To analyze the performance with different number of cores
```
./scalability.sh
```


## Directory Structure and Code Explanation
- [include/](include/) and [src/](src/) contain rispectively all the header and source files for the project.  
- [output/](output/): Contains the VTK files for the output of the mesh and the performace file .txt that compares the execution time with different grid sizes and number of cores.  

- [ApproxSolution](src/approxSolution.cpp)  
Here is implemented the main class. The data is saved in a **vector of doubles**, instead the main functions are **iterate**, which does an iteration of the Jacobi method, and **solve**, that calles iterate sequentially until convergence. Inside iterate, where the main computations are computed, is applied a **omp parallel for** which is able to cut in half the execution time (for my hardware).  

- [Muparser](include/muparser_fun.hpp)  
This header is needed to create the class able to read a string and make it into a muparser element. This will be called by the read methods.  

- [Reader](src/readJson.cpp)  
In here are defined the functions that interact with the **input .json file**, there are many different functions with basic interpretation and the main function **readDataJson()** is the one used directly in the constructor of approsSolution.  

- [Solver](src/LaplaceSolver.cpp)  
Here you can find three functions for the computation of the solution, all functions compute with chrono the execution time and print it.  
**sequentialSolver()** simply constructs the approximate, calls the solve() and finally computes the error with the exact solution.  
**parallelSolver()** is a bit more complex. The idea is to construct the initial data based on the input values, then split with a scatterv (to split as evenly as possible the data) into local approximations, follows the iteration method on all processors (with care on the sharing of the local boundary conditions) and finally a gatherv.  
Moreover, since this is the most efficient method, the call for the saving of the VTK file for Paraview is here, this is done if the flag in the .json file is true.  
**parallelPerformance()** is used to analyze the compared performance between sequential and parallel execution as the grid size increases.  
The mesh size starts form the input "GridDimension" and at each iteration doubles. To avoid long and expensive tests I set the test to stop when it reaches a dimension of 300. If needed this can be changed inside the function source file.  

- [Utilities](src/utilities.cpp) and [Output](include/writeVTK.hpp)  
Here the functions are pretty obvious. Used for computing the **error** between an approxSolution element and a fucntion of doubles, **printing** a function given a grid size and writing the **output** in VTK format which is compatible with Paraview.  
Even if is quite obvious, I want to precise that the printing functions, both in utilities.cpp and in approxSolution.cpp, are set to print only the first 4 rows and columns for big matrixes, this can be changed inside the functions.  

- [Scalability](scalability.sh)
To test the **scalability** on your computer, you can run ./scalability.sh this will compile only the part of the main flagged by the environmental variable SCALABILITY_TEST which runs the performance test and saves in scalabilityResult.txt the sequential and parallel execution time of the program on the number of processors given.  
In the file myScalabilityResults.txt you can find the results obtained with my computer. I am using a docker container with 8 CPUs.  
To modify the number of processors, simply change the varibale "processors" in the .sh file.  
Remeber that after compiling the code with the scalability flag, for the other part of the code to be compiled it's necessary to run **make clan** before the **make** or equivalently **make -B**. This way the part of the code outside the ifdef is compiled.
#ifndef LAPLACESOLVER_HPP
#define LAPLACESOLVER_HPP

#include <string>

// Returns the time for the performance test
double sequentialSolver(std::string filename, std::size_t N=0);
double parallelSolver(std::string filename, std::size_t N=0);

std::vector<std::vector<double>> parallelPerformance(std::string filename, std::vector<int> &gridSizes);

#endif // LAPLACESOLVER_HPP
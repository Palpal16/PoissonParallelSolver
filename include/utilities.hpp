#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <functional>
#include <iostream>

#include "approxSolution.hpp"

double seqL2_error(std::function<double(double, double)> &u_ex, approxSolution& u_h);

void printExact(std::function<double(double, double)> &u_ex, std::size_t N);


#endif // UTILITIES_HPP
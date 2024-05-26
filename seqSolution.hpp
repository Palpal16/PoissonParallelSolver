#ifndef SEQSOLUTION_HPP
#define SEQSOLUTION_HPP

#include <array>
#include <functional>
#include <cmath>
#include <iostream>

template <std::size_t N>
class seqSolution{
	static constexpr double h = 1.0 / (N - 1);
	std::array<std::array<double,N>,N> data; //This is ordered row-wise
	std::function<double(double, double)> force; 
	//Maybe would be better to store the matrix instead of evaluating everytime??
	double boundaryValue = 0.0;
	// default set the function to return the boundary condition when it is homogeneous
	std::function<double(double, double)> boundaryFunction = [this](double x, double y) {return boundaryValue;};
	double tol = 1e-5;
	std::size_t maxIter = 5000;

public:
	// Constructor that sets central elements to 0 and boundary elements to the values given
	seqSolution(double boundary, std::function<double(double, double)> forcingTerm);
	// Constructor with non-homogenous boundary conditions
	seqSolution(std::function<double(double, double)> boundary, std::function<double(double, double)> forcingTerm);
	// Copy constructor
	seqSolution(const seqSolution& other);

	// Remember that indexes are row-wise and start from 0
	
	double operator()(const std::size_t& row, const std::size_t& col) const {return data[row][col];}
	// Non const version not needed for now
	//double& operator()(const std::size_t& row, const std::size_t& col){return data[row][col];}

	void setTolerance(double tolerance){ tol = tolerance;}
	void setMaxIter(std::size_t max){ maxIter = max; }

	std::size_t solve();
	double seqIterate();

	void print();
};

#include "seqSolutionImpl.hpp"

#endif //SEQSOLUTION_HPP
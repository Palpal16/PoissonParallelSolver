#ifndef APPROXSOLUTION_HPP
#define APPROXSOLUTION_HPP

#include<vector>
#include <functional>
#include <cmath>
#include <iostream>

#include "readJson.hpp"


class approxSolution{
	std::size_t nCols;
	double h;
	std::vector<double> data; //This is ordered row-wise
	std::vector<double> force;
	double tol;
	std::size_t maxIter;

public:
	// Constructor that reads the data from a json file
	approxSolution(std::string filename);
	//Constructor that takes only the rows chosen, from the passed data and force matirx, needed for the scatter for parallelizatoin
	approxSolution(std::vector<double> Data, std::vector<double> Force, std::size_t nCols, double tolerance, std::size_t iter): nCols(nCols), h(1.0/(nCols-1)), data(Data), force(Force), tol(tolerance), maxIter(iter){};
	// Copy constructor
	approxSolution(const approxSolution& other): nCols(other.nCols), h(other.h), data(other.data), force(other.force), tol(other.tol), maxIter(other.maxIter){};

	double operator()(const std::size_t& col, const std::size_t& row) const {return data[nCols*col+row];} // x coordinate in vertical direction
	// Non const version not needed for now
	//double& operator()(const std::size_t& col, const std::size_t& row){return data[nCols*row+col];}

	void setTolerance(double tolerance){ tol = tolerance;}
	double getTolerance() const {return tol;}
	void setMaxIter(std::size_t max){ maxIter = max;}
	std::size_t getMaxIter() const {return maxIter;}
	std::size_t getGridSize() const {return nCols;}

	void setRow(std::size_t n,const std::vector<double>& Row);
	std::vector<double> getRow(std::size_t n);

	std::vector<double> getData() const {return data;}
	std::vector<double> getForce() const {return force;}

	// Returns the number of iteration completed
	std::size_t solve();
	// Returns the error from the previous iteration
	double iterate();

	void print();
};


#endif //APPROXSOLUTION_HPP
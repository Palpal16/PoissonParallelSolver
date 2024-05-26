#include <iostream>
#include <cmath>
#include "seqSolution.hpp"

template <std::size_t N>
double seqL2_error(std::function<double(double, double)> &u_ex, seqSolution<N>& u_h);

template <std::size_t N>
void printExact(std::function<double(double, double)> &u_ex);

// Test function of the sequential solver
void sequentialTest(){
	// Set number of vertices for each dimension
	constexpr std::size_t GridDimensions = 100;

	// Set the tolerance and maximum number of iterations
	double tol = 1e-5;
	std::size_t maxIter = 5000;

	// Set the forcing term and exact solution
	std::function<double(double, double)> force = [](double x, double y) {
    return 8 * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y); };
	std::function<double(double, double)> exactSol = [](double x, double y) {
    return std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y); };

	// Set the Dirichlet boundary conditions, either homogeneous or non-homogeneous
	double boundaryValue=0.0;
	std::function<double(double, double)> boundaryFunction = [](double x, double y) {return 0.0;};

	// Generate data. Remember to change the boundary condition if you want to test non-homogeneous boundary conditions
	seqSolution<GridDimensions> test(boundaryFunction,force);

	if(tol!=1e-5)
		test.setTolerance(tol);
	if(maxIter!=5000)
		test.setMaxIter(maxIter);

	std::size_t iter = test.solve();
	std::cout << "Number of iterations: " << iter << std::endl;
	std::cout << "Approximate solution: " << std::endl;
	test.print();

	std::cout << "\nL2 error: " << seqL2_error<GridDimensions>(exactSol, test) << std::endl;

	std::cout << "\nExact solution: " << std::endl;
	printExact<GridDimensions>(exactSol);
}


template <std::size_t N>
void printExact(std::function<double(double, double)> &u_ex){
	double h=1.0/(N-1);
	std::size_t n= N<8 ? N : 8; // To avoid printing too much
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			std::cout << u_ex(h*i, h*j) << " ";
		}
		std::cout << std::endl;
	}
	if (N>8){
		for (int i = 0; i < n; i++){
			std::cout << "...   ";
		}
		std::cout << std::endl;
	}
}

template <std::size_t N>
double seqL2_error(std::function<double(double, double)> &u_ex, seqSolution<N>& u_h){
	double error = 0.0;
	double h=1.0/(N-1);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			error +=(u_ex(h*i, h*j) - u_h(i,j))*(u_ex(h*i, h*j) - u_h(i,j));
		}
	}
	return std::sqrt(h*error);
}
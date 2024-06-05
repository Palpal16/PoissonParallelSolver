#include "utilities.hpp"

double seqL2_error(std::function<double(double, double)> &u_ex, approxSolution& u_h){
    std::size_t N=u_h.getGridSize();
	double error = 0.0;
	double h=1.0/(N-1);
	for (std::size_t i = 0; i < N; i++){
		for (std::size_t j = 0; j < N; j++){
			error +=(u_ex(h*i, h*j) - u_h(i,j))*(u_ex(h*i, h*j) - u_h(i,j));
		}
	}
	return std::sqrt(h*error);
}

void printExact(std::function<double(double, double)> &u_ex, std::size_t N){
	// This function is compatible only with square grids
	double h=1.0/(N-1);
	std::size_t MaxRows = 4;  // Change here to print a different number of rows
    bool MaxPrint = N>MaxRows;
	std::size_t n= MaxPrint ? MaxRows : N; // To avoid printing too much
	for (std::size_t i = 0; i < n; i++){
		for (std::size_t j = 0; j < n; j++){
			std::cout << u_ex(h*i, h*j) << " ";
		}
        if (MaxPrint)
            std::cout << "...";
		std::cout << std::endl;
	}
	if (MaxPrint){
		for (std::size_t i = 0; i < n; i++){
			std::cout << " ⁝     ";
		}
		std::cout << std::endl;
	}
}
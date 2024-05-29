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
	double h=1.0/(N-1);
    bool Max_print = N>4;
	std::size_t n= Max_print ? 4 : N; // To avoid printing too much
	for (std::size_t i = 0; i < n; i++){
		for (std::size_t j = 0; j < n; j++){
			std::cout << u_ex(h*i, h*j) << " ";
		}
        if (Max_print)
            std::cout << "...";
		std::cout << std::endl;
	}
	if (Max_print){
		for (std::size_t i = 0; i < n; i++){
			std::cout << " ⁝     ";
		}
		std::cout << std::endl;
	}
}
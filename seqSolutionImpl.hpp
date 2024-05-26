#ifndef SEQSOLUTIONIMPL_HPP
#define SEQSOLUTIONIMPL_HPP

#include "seqSolution.hpp"

template <std::size_t N>
seqSolution<N>::seqSolution(double boundary, std::function<double(double, double)> forcingTerm): boundaryValue(boundary), force(forcingTerm){
    for(auto& row : data)
        row.fill(0.0);
    if (boundary!=0.0){
        data[0].fill(boundary);
        data[N-1].fill(boundary);
        for (int i = 1; i < N-1; i++){
            data[i][0] = boundary;
            data[i][N-1] = boundary;
        }
    }
}

template <std::size_t N>
seqSolution<N>::seqSolution(std::function<double(double, double)> boundary, std::function<double(double, double)> forcingTerm): boundaryFunction(boundary), force(forcingTerm){
    for(auto& row : data)
        row.fill(0.0);
    for (int i = 0; i < N; i++){
        data[0][i] = boundaryFunction(h*i, 0);
        data[N-1][i] = boundaryFunction(h*i, 1);
        data[i][0] = boundaryFunction(0, h*i);
        data[i][N-1] = boundaryFunction(1, h*i);
    }
}

template <std::size_t N>
seqSolution<N>::seqSolution(const seqSolution& other) {
    data = other.data;
    boundaryValue = other.boundaryValue;
    force = other.force;
}


template <std::size_t N>
std::size_t
seqSolution<N>::solve(){
    std::size_t iter = 0;
    double error = 2*tol;
    while (error > tol && iter < maxIter){
        error = seqIterate();
        iter++;
    }
    return iter;
}

template <std::size_t N>
double
seqSolution<N>::seqIterate(){
    std::array<std::array<double,N>,N> newData(data);
    // Loop over all points, except boundaries
    for (int i = 1; i < N-1; i++){
        for (int j = 1; j < N-1; j++){
            newData[i][j] = 0.25 * (data[i-1][j] + data[i+1][j] + data[i][j-1] + data[i][j+1] + h*h*force(h*(i-1), h*(j-1)));
            // In (i=0,j=0) the point is (x=0, y=0), in (i=1,j=0) the point is (x=h, y=0), etc.
        }
    }

    // Compute the error between the current and the previous iteration
    double error = 0.0;
    for (int i = 1; i < N-1; i++){
        for (int j = 1; j < N-1; j++){
            error += (data[i][j] - newData[i][j])*(data[i][j] - newData[i][j]);
        }
    }
    data = newData;
    return std::sqrt(h*error);
}

template <std::size_t N>
void
seqSolution<N>::print(){
    std::size_t n= N<8 ? N : 8; // To avoid printing too much
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            std::cout << data[i][j] << " ";
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


#endif  // SEQSOLUTIONIMPL_HPP
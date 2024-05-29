#include "approxSolution.hpp"
#include <omp.h>

approxSolution::approxSolution(std::string filename){
    std::function<double(double, double)> forcingTerm;
    bool Dir;
    double boundaryValue;
    std::function<double(double, double)> boundaryFunction;
    // Read the data from the json file
    readDataJson(filename, nCols, tol, maxIter, forcingTerm, Dir, boundaryValue, boundaryFunction);
    h = 1.0 / (nCols - 1);

    // Set all elements to 0
    data = std::vector<double>(nCols*nCols, 0.0);

    // Resize the vector to avoid reallocation
    force.resize(nCols*nCols);
    for (std::size_t i = 0; i < nCols; i++){
        for (std::size_t j = 0; j < nCols; j++){
            // Evaluate if we are at the boundary
            if(i==0 || i==nCols-1 || j==0 || j==nCols-1){
                // Set the boundary values base on the homogenoeus choice
                if(Dir && boundaryValue!=0.0)
                    data[i*nCols+j] = boundaryValue;
                else if(!Dir)
                    data[i*nCols+j] = boundaryFunction(h*j, h*i);
            }
            // Evaluate the forcing terms in all points
            // (we don't need the values of f at the boundary, but we keep them for simplicity)
            force[i*nCols+j] = forcingTerm(h*j, h*i);
            // Points are defined as follows:
            // In (i=0,j=0) the point is (x=0, y=0), in (i=1,j=0) the point is (x=0, y=h), etc.
            // This is done to get "visually" the x-axis along the horizontal axis
        }
    }
}

void approxSolution::setRow(std::size_t n,const std::vector<double>& Row){
    for (std::size_t i = 0; i < nCols; i++){
        data[n*nCols+i] = Row[i];
    }
}
std::vector<double> approxSolution::getRow(std::size_t n){
    std::vector<double> Row(nCols);
    for (std::size_t i = 0; i < nCols; i++){
        Row[i] = data[n*nCols+i];
    }
    return Row;
}

std::size_t approxSolution::solve(){
    std::size_t iter = 0;
    double error = 2*tol;
    while (error > tol && iter < maxIter){
        error = iterate();
        iter++;
    }
    return iter;
}

double approxSolution::iterate(){
    std::size_t rows = data.size()/nCols;
    std::vector<double> newData(data);
    double error = 0.0;
    
    // Loop over all points, except boundaries
    #pragma omp parallel for shared(newData) reduction(+:error)
    for (std::size_t i = 1; i < rows-1; i++){
        for (std::size_t j = 1; j < nCols-1; j++){
            // In (i=0,j=0) the point is (x=0, y=0), in (i=1,j=0) the point is (x=h, y=0), etc.
            newData[i*nCols+j] = 0.25 * (data[(i-1)*nCols+j] + data[(i+1)*nCols+j] + data[i*nCols+j-1] + data[i*nCols+j+1] + h*h*force[i*nCols+j]);
            // Compute the error between the current and the previous iteration
            error += (data[i*nCols+j] - newData[i*nCols+j])*(data[i*nCols+j] - newData[i*nCols+j]);
        }
    }
    data = newData;
    return std::sqrt(h*error);
}


void approxSolution::print(){    
    bool Max_print = nCols>4;
    std::size_t cols= Max_print ? 4 : nCols; // To avoid printing too much
    std::size_t rows = data.size()/nCols; // The matrix might not be a square during parallel execution
    bool Max_print_rows= rows>4;
    rows= (Max_print_rows) ? 4 : rows;
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++){
            std::cout << data[i*nCols+j] << " ";
        }
        if (Max_print)
            std::cout << "...";
        std::cout << std::endl;
    }
    if (Max_print_rows){
        for (std::size_t i = 0; i < cols; i++){
            std::cout << " ⁝     ";
        }
        std::cout << std::endl;
    }
}
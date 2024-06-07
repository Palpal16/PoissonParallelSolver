#include <vector>
#include <iostream>
#include <mpi.h>
#include <chrono>

#include "LaplaceSolver.hpp"
#include "approxSolution.hpp"
#include "utilities.hpp"
#include "readJson.hpp"
#include "writeVTK.hpp"


double sequentialSolver(std::string filename, std::size_t N){

    approxSolution u(filename, N);
    std::function<double(double, double)> u_ex = readExactJson(filename);

    auto start = std::chrono::high_resolution_clock::now();
    std::size_t iter = u.solve();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;

    std::cout << " ------------ Sequential approximate solution ------------ "<< std::endl;
    std::cout << "Sequential number of iterations: " << iter << std::endl;
    std::cout << "Sequential execution time: " << diff.count() << " s"<< std::endl;
    std::cout << "Sequential L2 error: " << seqL2_error(u_ex,u) << std::endl;

    if(readIfPrintJson(filename)){
        std::cout << "\nSequential approximate solution: " << std::endl;
        u.print();
    }
    std::cout << std::endl;

    return diff.count();
}


double parallelSolver(std::string filename, std::size_t N){

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int GridDimensions;
    double tolerance;
    std::size_t maxIter;

    std::vector<double> data;
    std::vector<double> force;

    if(rank==0){
        approxSolution u(filename, N);
        GridDimensions = u.getGridSize();
        tolerance = u.getTolerance();
        maxIter = u.getMaxIter();

        data = u.getData();
        force = u.getForce();
    }
    MPI_Bcast(&GridDimensions, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxIter, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    // Scatter the data and the force to all processes
    // The first and last rows are not to be computed, so we don't consider them in the split of the data
    // Moreover, we need to add a row above and below the local data to store the "local boundaries"

    auto start = std::chrono::high_resolution_clock::now();

    int global_n = GridDimensions-2;
    int local_n = (global_n % size > rank) ? global_n/size +1 : global_n/size;
    local_n = local_n + 2;

    std::vector<double> local_data(local_n*GridDimensions);
    std::vector<double> local_force(local_n*GridDimensions);

    std::vector<int> send_counts(size, 0), send_start_idx(size, 0);
    // For the send and start we consider the number of rows, then we multiply for the columns
    if(rank==0){
        for (int i = 0; i < size; i++){
            send_counts[i] = (global_n % size > i) ? global_n/size +1 : global_n/size;
            send_counts[i] = send_counts[i] + 2;
            send_counts[i] = send_counts[i]*GridDimensions;
            // First process should start from first row, the others from the last row of the previous process
            send_start_idx[i] = i==0 ? 0 : send_start_idx[i-1] + send_counts[i-1]-2*GridDimensions;
        }
    }
    MPI_Scatterv(data.data(), send_counts.data(), send_start_idx.data(), MPI_DOUBLE, local_data.data(), local_n*GridDimensions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(force.data(), send_counts.data(), send_start_idx.data(), MPI_DOUBLE, local_force.data(), local_n*GridDimensions, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    approxSolution local_u(local_data, local_force, GridDimensions, tolerance, maxIter);

    bool converged = false;
    std::size_t iter = 0;
    while(!converged && iter < maxIter){
        double local_error = local_u.iterate();
        bool local_converged = local_error < tolerance;
        MPI_Allreduce(&local_converged, &converged, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        iter++;
        
        // Send the second-last row to the next process and receive from the next process into the last row
        std::vector<double> send_buffer = local_u.getRow(local_n-2);
        std::vector<double> recv_buffer(GridDimensions);
        if (rank < size - 1) {
            MPI_Sendrecv(send_buffer.data(), GridDimensions, MPI_DOUBLE, rank + 1, 0,
                        recv_buffer.data(), GridDimensions, MPI_DOUBLE, rank + 1, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_u.setRow(local_n-1, recv_buffer);  // Set the last row
        }

        // Send the second row to the previous process and receive from the previous process into the first row
        send_buffer = local_u.getRow(1);
        if (rank > 0) {
            MPI_Sendrecv(send_buffer.data(), GridDimensions, MPI_DOUBLE, rank - 1, 1,
                        recv_buffer.data(), GridDimensions, MPI_DOUBLE, rank - 1, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_u.setRow(0, recv_buffer); // Set the first row
        }
    }

    local_data = local_u.getData();

    // Gather the local data to the root process
    // I want to receive the data into the vector data
    // The first and last rows are unchanged, so I don't consider them in the gather
    // From each process I receive local_n-2 rows, because the first and last rows are gathered from the processes above and below
    std::vector<int> recv_counts(size, 0), recv_start_idx(size, 0);
    if (rank == 0) {
        for (int i = 0; i < size; i++){
            recv_counts[i] = send_counts[i] - 2*GridDimensions;
            recv_start_idx[i] = send_start_idx[i] + GridDimensions;
        }
    }
    MPI_Gatherv(local_data.data() + GridDimensions, (local_n-2) * GridDimensions, MPI_DOUBLE,
                data.data(), recv_counts.data(), recv_start_idx.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

                
    if(rank == 0){
        approxSolution u(data, force, GridDimensions, tolerance, maxIter);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end-start;

        std::function<double(double, double)> u_ex = readExactJson(filename);
        std::cout << " ------------ Parallel approximate solution ------------ "<< std::endl;
        std::cout << "Parallel number of iterations: " << iter << std::endl;
        std::cout << "Parallel execution time: " << diff.count() << " s"<< std::endl;
        std::cout << "Parallel L2 error: " << seqL2_error(u_ex,u) << std::endl;

        if(readIfPrintJson(filename)){
            std::cout << "\nParallel approximate solution: " << std::endl;
            u.print();
        }

        if(readIfVtkJson(filename)){
            std::string vtkFile = readVtkNameJson(filename); 
            std::cout << " ------------ Writing VTK file to " << vtkFile << " ------------ "<< std::endl;
            generateVTKFile(vtkFile, data, GridDimensions);
        }
        return diff.count();
    }
    return 0.0;

}


std::vector<std::vector<double>> parallelPerformance(std::string filename, std::vector<int> &gridSizes){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int MaxGridSize = 50;
    bool finished=false;
    std::vector<std::vector<double>> Times(2, std::vector<double>());
    
    for (int i=0; ; i++){
        if(rank==0){
            if(i==0){            
                gridSizes.push_back(readDimJson(filename));
            }else if(gridSizes[i-1]*2 <= MaxGridSize){
                gridSizes.push_back(gridSizes[i-1]*2);
            }
            else{
                finished = true;
            }
            if(!finished){
                std::cout << "\n********************************************" << std::endl;
                std::cout << "Performance iteration "<< i << " with grid size " << gridSizes[i] << std::endl;
                std::cout << "********************************************" << std::endl;
                Times[0].push_back(sequentialSolver(filename, gridSizes[i]));
            }
        }
        // Broadcast the bool finished to all processes
        MPI_Bcast(&finished, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if(finished){
            break;
        }

        //Broadcast the grid size of this iteration to all processes
        int N;
        if(rank==0){
            N = gridSizes[i];
        }
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        double mom_time = parallelSolver(filename, N);
        if(rank==0){
            Times[1].push_back(mom_time);
        }
    }
    return Times;
}
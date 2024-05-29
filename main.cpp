#include <mpi.h>
#include <omp.h>

#include "LaplaceSolver.hpp"
#include "readJson.hpp"
#include "utilities.hpp"

int main(int argc, char** argv) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string filename = "data.json";

    if(rank==0 && readIfPrintJson(filename)){
        std::function<double(double, double)> u_ex = readExactJson(filename);
        std::cout << " ------------ Print of exact solution ------------ "<< std::endl;
        std::cout << "Exact solution: " << std::endl;
        printExact(u_ex, readDimJson(filename));
        std::cout << std::endl;
    }

    if (rank == 0 && readIfSeqJson(filename)) {
        sequentialSolver(filename);
    }

    parallelSolver(filename);

    MPI_Finalize();

    return 0;
}

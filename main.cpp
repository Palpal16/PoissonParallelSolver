#include <mpi.h>

#include "LaplaceSolver.hpp"
#include "readJson.hpp"

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string filename = "data.json";

    if (rank == 0 && readIfSeqJson(filename)) {
        sequentialSolver(filename);
    }

    parallelSolver(filename);

    MPI_Finalize();

    return 0;
}

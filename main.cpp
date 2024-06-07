#include <mpi.h>
#include <omp.h>
#include <vector>

#include "LaplaceSolver.hpp"
#include "readJson.hpp"
#include "utilities.hpp"

int main(int argc, char** argv) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string filename = "input/data.json";

// This part is run only for the scalability test
#ifdef SCALABILITY_TEST

    std::vector<int> gridSizes;
    
    // Sace the pointer to stdout
    std::streambuf* coutbuf = std::cout.rdbuf();

    // Redirect stdout to /dev/null
    // This way the output of the parallelPerformance() function is not printed
    std::ostringstream ostr;
    std::cout.rdbuf(ostr.rdbuf());

    std::vector<std::vector<double>> times = parallelPerformance(filename, gridSizes);

    // Redirect stdout back to the original
    std::cout.rdbuf(coutbuf);

    if(rank==0){
        std::cout << "Grid size, Seq time, Par time" << std::endl;
        for(std::size_t i=0; i<gridSizes.size(); i++){
            std::cout << gridSizes[i] << ", " << times[0][i] << ", " << times[1][i] << std::endl;
        }
    }

#else

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

    if(readIfMpiJson(filename)){
        parallelSolver(filename);
    }

    if(readIfPerformanceJson(filename)){
        std::vector<int> gridSizes;
        std::vector<std::vector<double>> times = parallelPerformance(filename, gridSizes);
        if(rank==0){
            std::cout << "\n ------------ Scalability test ------------ "<< std::endl;
            std::cout << "Grid size, Seq time, Par time" << std::endl;
            for(std::size_t i=0; i<gridSizes.size(); i++){
                std::cout << gridSizes[i] << ", " << times[0][i] << ", " << times[1][i] << std::endl;
            }
        }
    }

#endif

    MPI_Finalize();

    return 0;
}

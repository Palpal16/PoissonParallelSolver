#include <fstream>

#include "readJson.hpp"
#include "json.hpp"
using json = nlohmann::json;


void readDataJson(const std::string & filename, std::size_t &N, double &tol, std::size_t &maxIter,
				std::function<double(double, double)> &force, bool &Dir, double &boundaryValue,
				std::function<double(double, double)> &boundaryFunction){
	std::ifstream f(filename);
	json data = json::parse(f);

	// Initializing parameters from JSON data
	N = data.value("GridDimension",1);
	tol = data.value("tolerance", 1e-3);
    maxIter = data.value("maxIterations",10);

	// Initializing forcing term
	std::string forceString = data.value("forceFunction","");
	MuparserFun Force(forceString, 2);
	force=Force;
	
	// Initializing boundary conditions
	Dir = data.value("homogeneousDirichlet",false);
	boundaryValue=data.value("boundaryValue",1.0);
	std::string boundaryFunctionString = data.value("boundaryFunction","");
	MuparserFun BoundaryFunction(boundaryFunctionString, 2);
	boundaryFunction=BoundaryFunction;
}

std::function<double(double, double)> readExactJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	// Initializing exact solution
	std::string exactSolString = data.value("exactSolution","");
	MuparserFun ExactSol(exactSolString, 2);
	return ExactSol;
}

int readDimJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	int dim = data.value("GridDimension",1);
	return dim;
}

bool readIfPrintJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	bool print = data.value("printSolution",false);
	return print;
}

bool readIfSeqJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	bool seq = data.value("solveSequential",false);
	return seq;
}

bool readIfMpiJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	bool mpi = data.value("solveParallel",false);
	return mpi;
}

bool readIfPerformanceJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	bool perf = data.value("performanceTest",false);
	return perf;
}

bool readIfVtkJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	bool vtk = data.value("WriteVTK",false);
	return vtk;
}

std::string readVtkNameJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	std::string vtkName = data.value("outputFile","output.vtk");
	return vtkName;
}
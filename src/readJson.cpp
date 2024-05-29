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
	N = data.value("GridDimension",10);
	tol = data.value("tolerance", 1e-7);
    maxIter = data.value("maxIterations",100);

	// Initializing forcing term
	std::string forceString = data.value("forceFunction","");
	MuparserFun Force(forceString, 2);
	force=Force;
	
	// Initializing boundary conditions
	Dir = data.value("homogeneousDirichlet",1);
	boundaryValue=data.value("boundaryValue",0.0);
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

bool readIfPrintJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	// Read the dimesion of the grid
	bool print = data.value("printSolution",true);
	return print;
}

bool readIfSeqJson(const std::string & filename){
	std::ifstream f(filename);
	json data = json::parse(f);

	// Read the dimesion of the grid
	bool seq = data.value("solveSequential",true);
	return seq;
}
#ifndef READJSON_HPP
#define READJSON_HPP

#include <string>
#include <functional>

#include "muparser_fun.hpp"

void readDataJson(const std::string & filename, std::size_t &N, double &tol, std::size_t &maxIter,
				std::function<double(double, double)> &force, bool &Dir, double &boundaryValue,
				std::function<double(double, double)> &boundaryFunction);

std::function<double(double, double)> readExactJson(const std::string & filename);
int readDimJson(const std::string & filename);

bool readIfPrintJson(const std::string & filename);
bool readIfSeqJson(const std::string & filename);
bool readIfMpiJson(const std::string & filename);
bool readIfPerformanceJson(const std::string & filename);
bool readIfVtkJson(const std::string & filename);
std::string readVtkNameJson(const std::string & filename);

#endif // READJSON_HPP
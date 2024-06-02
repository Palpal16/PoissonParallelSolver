# FirstChallengePacs

## Gradient-Based Optimization in C++ for minimum problem

This is the repository for the first assignement of the APSC course.\
The C++ program implements gradient-based optimization algorithms such as gradient descent, Adam, and Heavy Ball. The program reads input parameters from a JSON file, performs optimization based on the specified method, and outputs the final result.


## Features

- Supports gradient descent, Adam, and Heavy Ball optimization methods.
- Input parameters are read from a JSON file for flexibility and ease of use.
- Implements finite differences for gradient computation if specified in the input file.
- Provides options for different step size update strategies, including fixed step size, exponential decay, and Armijo rule.

## Requirements

- C++ compiler with C++11 support.
- [JSON for Modern C++](https://github.com/nlohmann/json) library for parsing JSON files.
- [muParser](https://beltoforion.de/article.php?a=muparser) library for parsing mathematical expressions and computing derivatives.

## Usage

1. Clone the repository:

```
git clone git@github.com:Palpal16/FirstChallengePacs.git
```

2. Modify the variable L_ROOT in the Makefile

3. Modify the file data.json with parameters of your choice

4. Compile and run the program


```
make
./main
```

## The program
The aim of this is program is to be able to run tests with different methods with only one compilation.\
This will increase the run time, but it avoids having to compile again for different tests.

To implement the function and the gradient I used the muParser library, since it's slighly more efficient than the muParserX and I didn't need the enhancements of the latter library.

Inside Data.hpp you can find the code that reads data from the json file and a few smaller functions. Instead in the main there are the three approximation methods.


## Parameters
To tune the program and run different tests you can modify the following variables:

1. Choose the minimization algorithm
```
"learningMethod": gradient descent = 0, Adam method = 1, heavy ball method = 2.
```
If using Adam or heavy ball remember to set the respective hyperparameters:\
"gamma1AD", "gamma2AD", "epsilonAD": 1e-8 or "etaHB": 0.01.



3. Choose the learning decay method
```
"alphaMethod":   Fixed learning = 0, Exponential decay = 1, Inverse decay = 2, Armijo rule = 3.
```
Remember not to use Armijo rule for heavy ball and Adam methods.\
Also set the starting value of the learning rate "alpha0" and the "decayParameter", this is usually referred to as mu for exponential deacy and sigma for Armijo rule.

5. Choose the gradient evaluation technique\
For FD set "finiteDifferences" to 1, instead for the exact gradient set it to any other integer.\
If using the finite difference scheme, you should select an appropriate step h and set it in "hFD"\

6. General parameters\
Starting point = "initial". Maximum number of iterations = "max_it". tolerance w.r.t. x and to the gradient = "tol_x" and "tol_df". The function to be minimized and its gradient = "fun" and "dfun".

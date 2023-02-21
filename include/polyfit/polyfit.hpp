#ifndef POLYFIT_POLYFIT_HPP
#define POLYFIT_POLYFIT_HPP


// std lib
#include <cstddef>
#include <string>


namespace PolyFit {


// Adapted from https://github.com/codeplea/incbeta
double incbeta(double a, double b, double x);

double invincbeta(double y,double alpha, double beta);


// Calculate the t value for a Student distribution
// Adapted from http://www.cplusplus.com/forum/beginner/216098/
// **************************************************************
double CalculateTValueStudent(const double nu, const double alpha);

// Cumulative distribution for Student-t
// **************************************************************
double cdfStudent(const double nu, const double t);

// Cumulative distribution for Fisher F
// **************************************************************
double cdfFisher(const double df1, const double df2, const double x);

// Initialize a 2D array
// **************************************************************
double **Make2DArray(const std::size_t rows, const std::size_t cols);

// Transpose a 2D array
// **************************************************************
double **MatTrans(double **array, const std::size_t rows, const std::size_t cols);

// Perform the multiplication of matrix A[m1,m2] by B[m2,m3]
// **************************************************************
double **MatMul(const std::size_t m1, const std::size_t m2, const std::size_t m3, double **A, double **B);

// Perform the multiplication of matrix A[m1,m2] by vector v[m2,1]
// **************************************************************
void MatVectMul(const std::size_t m1, const std::size_t m2, double **A, double *v, double *Av);


// Calculates the determinant of a matrix
// **************************************************************
double determinant(double **a, const std::size_t k);


// Perform the
// **************************************************************
void transpose(double **num, double **fac, double **inverse, const std::size_t r);

// Calculates the cofactors
// **************************************************************
void cofactor(double **num, double **inverse, const std::size_t f);


// Display a matrix
// **************************************************************
void displayMat(double **A, const std::size_t n, const std::size_t m);

// Calculate the residual sum of squares (RSS)
// **************************************************************
double CalculateRSS(const double *x, const double *y, const double *a, double **Weights,
const bool fixed, const std::size_t N, const std::size_t n);

// Calculate the total sum of squares (TSS)
// **************************************************************
double CalculateTSS(const double *x, const double *y, const double *a, double **Weights,
const bool fixed, const std::size_t N, const std::size_t n);

// Calculate coefficient R2 - COD
// **************************************************************
double CalculateR2COD(const double *x, const double *y, const double *a, double **Weights,
const bool fixed, const std::size_t N, const std::size_t n);

// Calculate the coefficient R2 - adjusted
// **************************************************************
double CalculateR2Adj(const double *x, const double *y, const double *a, double **Weights,
const bool fixed,const std::size_t N, const std::size_t n);

// Perform the fit of data n data points (x,y) with a polynomial of order k
// **************************************************************
void PolyFit(const double *x, double *y, const std::size_t n, const std::size_t k, const bool fixedinter,
const double fixedinterval, double *beta, double **Weights, double **XTWXInv);

// Calculate the polynomial at a given x value
// **************************************************************
double calculatePoly(const double x, const double *a, const std::size_t n);

// Calculate and write the confidence bands in a file
// **************************************************************
void WriteCIBands(std::string filename, const double *x, const double *coefbeta, double **XTXInv,
const double tstudentval, const double SE, const std::size_t n, const std::size_t k);

// Calculate the weights matrix
// **************************************************************
void CalculateWeights(const double *erry, double **Weights, const std::size_t n,
const int type);

// Calculate the standard error on the beta coefficients
// **************************************************************
void CalculateSERRBeta(const bool fixedinter, const double SE, std::size_t k, double *serbeta, double **XTWXInv);

// Display the polynomial
// **************************************************************
void DisplayPolynomial(const std::size_t k);

// Display the ANOVA test result
// **************************************************************
void DisplayANOVA(const std::size_t nstar, const std::size_t k, const double TSS, const double RSS);


// Display the coefficients of the polynomial
// **************************************************************
void DisplayCoefs(const std::size_t k, const std::size_t nstar, const double tstudentval, const double *coefbeta, const double *serbeta);

// Display some statistics values
// **************************************************************
void DisplayStatistics(const std::size_t n, const std::size_t nstar, const std::size_t k, const double RSS, const double R2,
const double R2Adj, const double SE);


// Display the covariance and correlation matrix
// **************************************************************
void DisplayCovCorrMatrix(const std::size_t k, const double sigma, const bool fixed, double **XTWXInv);


}  // PolyFit namespace


#endif  // POLYFIT_POLYFIT_HPP define

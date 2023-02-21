// project lib
#include "polyfit/polyfit.hpp"
using namespace PolyFit;


// std lib
#include <cmath>
#include <iostream>
using namespace std;


// The main program
// **************************************************************
int main(int argc, char *argv[]) {

    cout << "Polynomial fit!" << endl;

    // Input values
    // **************************************************************
    size_t k = 2;                                    // Polynomial order
    bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
    int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
    double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
    double alphaval = 0.05;                          // Critical apha value

    double x[] = {0., 0.5, 1.0, 2.0, 4.0, 6.0};
    double y[] = {0., 0.21723, 0.43445, 0.99924, 2.43292, 4.77895};
    double erry[] = {0.1, 0.3, 0.2, 0.4, 0.1, 0.3};       // Data points (err on y) (if applicable)

    // Definition of other variables
    // **************************************************************
    size_t n = 0;                                    // Number of data points (adjusted later)
    size_t nstar = 0;                                // equal to n (fixed intercept) or (n-1) not fixed
    double coefbeta[k+1];                            // Coefficients of the polynomial
    double serbeta[k+1];                             // Standard error on coefficients
    double tstudentval = 0.;                         // Student t value
    double SE = 0.;                                  // Standard error

    double **XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
    double **Weights;                                // Matrix Weights [n,n]


    // Initialize values
    // **************************************************************
    n = sizeof(x)/sizeof(double);
    nstar = n-1;
    if (fixedinter) nstar = n;

    cout << "Number of points: " << n << endl;
    cout << "Polynomial order: " << k << endl;
    if (fixedinter) {
        cout << "A0 is fixed!" << endl;
    } else {
        cout << "A0 is adjustable!" << endl;
    }

    if (k>nstar) {
        cout << "The polynomial order is too high. Max should be " << n << " for adjustable A0 ";
        cout << "and " << n-1 << " for fixed A0. ";
        cout << "Program stopped" << endl;
        return -1;
    }

    if (k==nstar) {
        cout << "The degree of freedom is equal to the number of points. ";
        cout << "The fit will be exact." << endl;
    }

    XTWXInv = Make2DArray(k+1,k+1);
    Weights = Make2DArray(n,n);

    // Build the weight matrix
    // **************************************************************
    CalculateWeights(erry, Weights, n, wtype);

    cout << "Weights" << endl;
    displayMat(Weights,n,n);

    if (determinant(Weights,n)==0.) {
        cout << "One or more points have 0 error. Review the errors on points or use no weighting. ";
        cout << "Program stopped" << endl;
        return -1;
    }

    // Calculate the coefficients of the fit
    // **************************************************************
    PolyFit::PolyFit(x,y,n,k,fixedinter,fixedinterval,coefbeta,Weights,XTWXInv);


    // Calculate related values
    // **************************************************************
    double RSS = CalculateRSS(x,y,coefbeta,Weights,fixed,n,k+1);
    double TSS = CalculateTSS(x,y,coefbeta,Weights,fixedinter,n,k+1);
    double R2 = CalculateR2COD(x,y,coefbeta,Weights,fixedinter,n,k+1);
    double R2Adj = CalculateR2Adj(x,y,coefbeta,Weights,fixedinter,n,k+1);

    if ((nstar-k)>0) {
        SE = sqrt(RSS/(nstar-k));
        tstudentval = fabs(CalculateTValueStudent(nstar-k, 1.-0.5*alphaval));
    }
    cout << "t-student value: " << tstudentval << endl << endl;

    // Calculate the standard errors on the coefficients
    // **************************************************************
    CalculateSERRBeta(fixedinter,SE,k,serbeta,XTWXInv);

    // Display polynomial
    // **************************************************************
    DisplayPolynomial(k);

    // Display polynomial coefficients
    // **************************************************************
    DisplayCoefs(k, nstar, tstudentval, coefbeta, serbeta);

    // Display statistics
    // **************************************************************
    DisplayStatistics(n,nstar,k,RSS,R2,R2Adj,SE);

    // Display ANOVA table
    // **************************************************************
    DisplayANOVA(nstar, k, TSS, RSS);

    // Write the prediction and confidence intervals
    // **************************************************************
    WriteCIBands("CIBands2.dat",x,coefbeta,XTWXInv,tstudentval,SE,n,k);

    // Display the covariance and correlation matrix
    // **************************************************************
    DisplayCovCorrMatrix(k, SE, fixedinter, XTWXInv);



}

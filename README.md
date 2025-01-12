# polyfit

A fork and refactor of [nasa/polyfit](https://github.com/nasa/polyfit) as a consumable C++ library.

The following changes were made:

- Separate header and source files
- CMake support

## Original Readme

```markdown
This is a stand alone C++ program to fit a curve with a polynomial.

To compile: g++ -o Polyfit Polyfit.cpp

Inputs:

k: Degree of the polynomial
fixedinter: Fixed the intercept (coefficient A0)
wtype: Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
alphaval: Critical apha value
x[]: Array of x values to be fitted
y[]: Array of y values to be fitted
erry[]: Array of error of y (if applicable)
```

## Original Licenses

```comment
// ********************************************************************
// * Code PolyFit                                                     *
// * Written by Ianik Plante                                          *
// *                                                                  *
// * KBR                                                              *
// * 2400 NASA Parkway, Houston, TX 77058                             *
// * Ianik.Plante-1@nasa.gov                                          *
// *                                                                  *
// * This code is used to fit a series of n points with a polynomial  *
// * of degree k, and calculation of error bars on the coefficients.  *
// * If error is provided on the y values, it is possible to use a    *
// * weighted fit as an option. Another option provided is to fix the *
// * intercept value, i.e. the first parameter.                       *
// *                                                                  *
// * This code has been written partially using data from publicly    *
// * available sources.                                               *
// *                                                                  *
// * The code works to the best of the author's knowledge, but some   *
// * bugs may be present. This code is provided as-is, with no        *
// * warranty of any kind. By using this code you agree that the      *
// * author, the company KBR or NASA are not responsible for possible *
// * problems related to the usage of this code.                      *
// *                                                                  *
// * The program has been reviewed and approved by export control for *
// * public release. However some export restriction exists. Please   *
// * respect applicable laws.                                         *
// *                                                                  *
// ********************************************************************


/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */
```

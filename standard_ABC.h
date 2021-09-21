/* ABC algorithm coded using C programming language */

/* Artificial Bee Colony (ABC) is one of the most recently defined algorithms by Dervis Karaboga in 2005, motivated by the intelligent behavior of honey bees. */

/* Reference Papers */

/* D. Karaboga, AN IDEA BASED ON HONEY BEE SWARM FOR NUMERICAL OPTIMIZATION,TECHNICAL REPORT-TR06, Erciyes University, Engineering Faculty, Computer Engineering Department 2005. */
/* D. Karaboga, B. Basturk, A powerful and Efficient Algorithm for Numerical Function Optimization: Artificial Bee Colony (ABC) Algorithm, Journal of Global Optimization, Volume:39, Issue:3,pp:459-171, November 2007, ISSN:0925-5001 , doi: 10.1007/s10898-007-9149-x */
/* D. Karaboga, B. Basturk, On The Performance Of Artificial Bee Colony (ABC) Algorithm, Applied Soft Computing,Volume 8, Issue 1, January 2008, Pages 687-697. */
/* D. Karaboga, B. Akay, A Comparative Study of Artificial Bee Colony Algorithm, Applied Mathematics and Computation, 214, 108-132, 2009. */

/* Copyright © 2009 Erciyes University, Intelligent Systems Research Group, The Dept. of Computer Engineering */

/* Contact:
   Dervis Karaboga (karaboga@erciyes.edu.tr )
   Bahriye Basturk Akay (bahriye@erciyes.edu.tr)
*/

// 2021-09-21
// modified by Choi, T. J.

#ifndef ABC_H
#define ABC_H

#include <math.h>
#include <stdlib.h>

/* a function pointer returning double and taking a D-dimensional array as argument */
typedef double (*FunctionCallback)(double *const solution, const int D);

void standard_ABC(const int gnty_s, const double *const intl_lwr_bnd, const double *const intl_uppr_bnd, const double *const srch_lwr_bnd, const double *const srch_uppr_bnd, double (*objective_function)(double *const solution, const int D), void (*callback_function)(const double GlobalMin, const double *const GlobalParams), const int evlt_mx, const int arg_NP, const int arg_FoodNumber, const int arg_limit);

#endif

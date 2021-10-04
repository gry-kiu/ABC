// ABC algorithm
// modified by Choi, T. J.

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

#ifndef ABC_H
#define ABC_H

#include <math.h>
#include <stdlib.h>

void standard_ABC(const int D, const double *const lb_init, const double *const ub_init, const double *const lb_iter, const double *const ub_iter, double (*objective)(const double *const solution, const int D), void (*callback)(const double GlobalMin, const double *const GlobalParams), const int max_function_evaluations, const int NP, const int FoodNumber, const double limit);

#endif

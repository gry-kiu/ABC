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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "standard_ABC.h"

// objective
double sphere(const double *const solution, const int D);

// callback
void callback(const double GlobalMin, const double *const GlobalParams);

int main(int argc, char *argv[]) {
  srand(time(NULL));
  const int D = 30;
  double *lb = (double *)malloc(sizeof(double) * D);
  double *ub = (double *)malloc(sizeof(double) * D);
  for (int i = 0; i < D; i++) {
    lb[i] = -100.0;
    ub[i] = 100.0;
  }

  standard_ABC(D, lb, ub, lb, ub, sphere, callback, 150000, 40, 20, 100);

  free(lb);
  free(ub);
  return 0;
}

double sphere(const double *const solution, const int D) {
  double top = 0;
  for (int j = 0; j < D; j++) {
    top = top + solution[j] * solution[j];
  }
  return top;
}

void callback(const double GlobalMin, const double *const GlobalParams) {
  printf("%e\n", GlobalMin);
}

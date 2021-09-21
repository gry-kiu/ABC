// 2021-09-21
// modified by Choi, T. J.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void standard_ABC(const int gnty_s, const double *const intl_lwr_bnd, const double *const intl_uppr_bnd, const double *const srch_lwr_bnd, const double *const srch_uppr_bnd, double (*objective_function)(double *const x, const int gnty_s), void (*callback_function)(const double GlobalMin, const double *const GlobalParams), const int evlt_mx, const int arg_NP, const int arg_FoodNumber, const int arg_limit);

// objective function
double sphere_function(double *const solution, const int D);

// callback function
void callback_function(const double GlobalMin, const double *const GlobalParams);

int main(int argc, char *argv[]) {
  srand(time(NULL));
  const int gnty_s = 30;
  double *lwr_bnd = (double *)malloc(sizeof(double) * gnty_s);
  double *uppr_bnd = (double *)malloc(sizeof(double) * gnty_s);
  for (int j = 0; j < gnty_s; j++) {
    lwr_bnd[j] = -100.0;
    uppr_bnd[j] = 100.0;
  }

  standard_ABC(gnty_s, lwr_bnd, uppr_bnd, lwr_bnd, uppr_bnd, sphere_function, callback_function, 150000, 40, 20, 100);

  free(lwr_bnd);
  free(uppr_bnd);
  return 0;
}

double sphere_function(double *const solution, const int D) {
  double top = 0;
  for (int j = 0; j < D; j++) {
    top = top + solution[j] * solution[j];
  }
  return top;
}

void callback_function(const double GlobalMin, const double *const GlobalParams) {
  printf("%e\n", GlobalMin);
}

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Control Parameters of ABC algorithm */
int NP; /* The number of colony size (employed bees+onlooker bees) */
int FoodNumber; /* The number of food sources equals the half of the colony size */
int limit; /* A food source which could not be improved through "limit" trials is abandoned by its employed bee */

/* Problem specific variables */
int D; /* The number of parameters of the problem to be optimized */
double *lb; /* lower bound of the parameters. */
double *ub; /* upper bound of the parameters. */

double **Foods; /* Foods is the population of food sources. */ // shape=(FoodNumber, D)
double *f; /* objective function values */ // shape=(FoodNumber,)
double *fitness; /* fitness (quality) values */ // shape=(FoodNumber,)
double *trial; /* trial numbers */ // shape=(FoodNumber,)
double *prob; /* probabilities of food sources (solutions) to be chosen */ // shape=(FoodNumber,)
double *solution; /* New solution (neighbor) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomly chosen solution different from i */ // shape=(D,)
double GlobalMin; /* Optimum solution obtained by ABC algorithm */
double *GlobalParams; /* Parameters of the optimum solution */ // shape=(D,)

/* a function pointer returning double and taking a D-dimensional array as argument */
typedef double (*FunctionCallback)(double *const solution, const int D);

/* Write your own objective function name instead of sphere */
FunctionCallback function;

/* Fitness function */
double CalculateFitness(double fun) {
  double result = 0;
  if (fun >= 0) {
    result = 1 / (fun + 1);
  } else {
    result = 1 + fabs(fun);
  }
  return result;
}

/* The best food source is memorized */
void MemorizeBestSource() {
  for (int i = 0; i < FoodNumber; i++) {
    if (f[i] < GlobalMin) {
      GlobalMin = f[i];
      for (int j = 0; j < D; j++)
        GlobalParams[j] = Foods[i][j];
    }
  }
}

/* Variables and counters are initialized */
void init(int index) {
  for (int j = 0; j < D; j++) {
    double r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    Foods[index][j] = r * (ub[j] - lb[j]) + lb[j];
    solution[j] = Foods[index][j];
  }
  f[index] = function(solution, D);
  fitness[index] = CalculateFitness(f[index]);
  trial[index] = 0;
}

/* All food sources are initialized */
void initial() {
  for (int i = 0; i < FoodNumber; i++) {
    init(i);
  }
  GlobalMin = f[0];
  for (int i = 0; i < D; i++)
    GlobalParams[i] = Foods[0][i];
}

void SendEmployedBees() {
  /* Employed Bee Phase */
  for (int i = 0; i < FoodNumber; i++) {
    /* The parameter to be changed is determined randomly */
    double r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    int param2change = (int)(r * D);

    /* A randomly chosen solution is used in producing a mutant solution of the solution i */
    r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    int neighbor = (int)(r * FoodNumber);

    /* Randomly selected solution must be different from the solution i */
    while (neighbor == i) {
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      neighbor = (int)(r * FoodNumber);
    }
    for (int j = 0; j < D; j++)
      solution[j] = Foods[i][j];

    /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
    r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    solution[param2change] = Foods[i][param2change] + (Foods[i][param2change] - Foods[neighbor][param2change]) * (r - 0.5) * 2;

    /* if generated parameter value is out of boundaries, it is shifted onto the boundaries */
    for (int j = 0; j < D; j++) {
      if (solution[param2change] < lb[j])
        solution[param2change] = lb[j];
      if (solution[param2change] > ub[j])
        solution[param2change] = ub[j];
    }
    double ObjValSol = function(solution, D);
    double FitnessSol = CalculateFitness(ObjValSol);

    /* a greedy selection is applied between the current solution i and its mutant */
    if (FitnessSol > fitness[i]) {
      /* If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i */
      trial[i] = 0;
      for (int j = 0; j < D; j++)
        Foods[i][j] = solution[j];
      f[i] = ObjValSol;
      fitness[i] = FitnessSol;
    } else { /* if the solution i can not be improved, increase its trial counter */
      trial[i] = trial[i] + 1;
    }
  }

  /* end of employed bee phase */
}

/* A food source is chosen with the probability which is proportional to its quality */
void CalculateProbabilities() {
  int i;
  double maxfit;
  maxfit = fitness[0];
  for (i = 1; i < FoodNumber; i++) {
    if (fitness[i] > maxfit)
      maxfit = fitness[i];
  }

  for (i = 0; i < FoodNumber; i++) {
    prob[i] = (0.9 * (fitness[i] / maxfit)) + 0.1;
  }
}

void SendOnlookerBees() {
  int i = 0;
  int t = 0;
  /* onlooker Bee Phase */
  while (t < FoodNumber) {
    double r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    /* choose a food source depending on its probability to be chosen */
    if (r < prob[i]) {
      t++;

      /* The parameter to be changed is determined randomly */
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      int param2change = (int)(r * D);

      /* A randomly chosen solution is used in producing a mutant solution of the solution i */
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      int neighbor = (int)(r * FoodNumber);

      /* Randomly selected solution must be different from the solution i */
      while (neighbor == i) {
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        neighbor = (int)(r * FoodNumber);
      }
      for (int j = 0; j < D; j++)
        solution[j] = Foods[i][j];

      /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      solution[param2change] = Foods[i][param2change] + (Foods[i][param2change] - Foods[neighbor][param2change]) * (r - 0.5) * 2;

      /* if generated parameter value is out of boundaries, it is shifted onto the boundaries */
      for (int j = 0; j < D; j++) {
        if (solution[param2change] < lb[j])
          solution[param2change] = lb[j];
        if (solution[param2change] > ub[j])
          solution[param2change] = ub[j];
      }
      double ObjValSol = function(solution, D);
      double FitnessSol = CalculateFitness(ObjValSol);

      /* a greedy selection is applied between the current solution i and its mutant */
      if (FitnessSol > fitness[i]) {
        /* If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i */
        trial[i] = 0;
        for (int j = 0; j < D; j++)
          Foods[i][j] = solution[j];
        f[i] = ObjValSol;
        fitness[i] = FitnessSol;
      } else { /* if the solution i can not be improved, increase its trial counter */
        trial[i] = trial[i] + 1;
      }
    } /* if */
    i++;
    if (i == FoodNumber)
      i = 0;
  } /* while */

  /* end of onlooker bee phase */
}

/* determine the food sources whose trial counter exceeds the "limit" value. */
void SendScoutBees() {
  int maxtrialindex = 0;
  for (int i = 1; i < FoodNumber; i++) {
    if (trial[i] > trial[maxtrialindex])
      maxtrialindex = i;
  }
  if (trial[maxtrialindex] >= limit) {
    init(maxtrialindex);
  }
}

void standard_ABC(const int gnty_s, const double *const intl_lwr_bnd, const double *const intl_uppr_bnd, const double *const srch_lwr_bnd, const double *const srch_uppr_bnd, double (*objective_function)(double *const solution, const int D), void (*callback_function)(const double GlobalMin, const double *const GlobalParams), const int evlt_mx, const int arg_NP, const int arg_FoodNumber, const int arg_limit) {
  // -------------------- initialization --------------------

  int evlt_c = 0;

  NP = arg_NP;
  FoodNumber = arg_FoodNumber;
  limit = arg_limit;

  function = objective_function;

  D = gnty_s;
  lb = (double *)malloc(sizeof(double) * D);
  ub = (double *)malloc(sizeof(double) * D);
  for (int i = 0; i < D; i++) {
    lb[i] = intl_lwr_bnd[i];
    ub[i] = intl_uppr_bnd[i];
  }

  Foods = (double **)malloc(sizeof(double *) * FoodNumber);
  for (int i = 0; i < FoodNumber; i++) {
    Foods[i] = (double *)malloc(sizeof(double) * D);
  }
  f = (double *)malloc(sizeof(double) * FoodNumber);
  fitness = (double *)malloc(sizeof(double) * FoodNumber);
  trial = (double *)malloc(sizeof(double) * FoodNumber);
  prob = (double *)malloc(sizeof(double) * FoodNumber);
  solution = (double *)malloc(sizeof(double) * D);
  GlobalParams = (double *)malloc(sizeof(double) * D);

  // -------------------- main procedure --------------------

  initial();
  MemorizeBestSource();

  for (int i = 0; i < D; i++) {
    lb[i] = srch_lwr_bnd[i];
    ub[i] = srch_uppr_bnd[i];
  }

  while (evlt_c < evlt_mx) {
    SendEmployedBees();
    evlt_c += FoodNumber;
    CalculateProbabilities();
    SendOnlookerBees();
    evlt_c += FoodNumber;
    MemorizeBestSource();
    SendScoutBees();
  }

  if (callback_function) {
    callback_function(GlobalMin, GlobalParams);
  }

  // -------------------- termination --------------------

  free(lb);
  free(ub);

  for (int i = 0; i < FoodNumber; i++) {
    free(Foods[i]);
  }
  free(Foods);
  free(f);
  free(fitness);
  free(trial);
  free(prob);
  free(solution);
  free(GlobalParams);
}

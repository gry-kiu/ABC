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

#include "standard_ABC.h"

/* Control Parameters of ABC algorithm */
static int g_NP; /* The number of colony size (employed bees+onlooker bees) */
static int g_FoodNumber; /* The number of food sources equals the half of the colony size */
static int g_limit; /* A food source which could not be improved through "g_limit" trials is abandoned by its employed bee */

/* Problem specific variables */
static int g_D; /* The number of parameters of the problem to be optimized */
static const double *g_lb; /* lower bound of the parameters. */
static const double *g_ub; /* upper bound of the parameters. */

static double **g_Foods; /* g_Foods is the population of food sources. */ // shape=(g_FoodNumber, g_D)
static double *g_f; /* objective function values */ // shape=(g_FoodNumber,)
static double *g_fitness; /* g_fitness (quality) values */ // shape=(g_FoodNumber,)
static double *g_trial; /* g_trial numbers */ // shape=(g_FoodNumber,)
static double *g_prob; /* probabilities of food sources (solutions) to be chosen */ // shape=(g_FoodNumber,)
static double *g_solution; /* New g_solution (neighbor) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomly chosen g_solution different from i */ // shape=(g_D,)
static double g_GlobalMin; /* Optimum g_solution obtained by ABC algorithm */
static double *g_GlobalParams; /* Parameters of the optmum g_solution */ // shape=(g_D,)

/* Write your own objective function name instead of sphere */
static double (*g_objective)(const double *const solution, const int D);

static double CalculateFitness(double fun);
static void MemorizeBestSource();
static void init(int index);
static void initial();
static void SendEmployedBees();
static void CalculateProbabilities();
static void SendOnlookerBees();
static void SendScoutBees();

void standard_ABC(const int D, const double *const lb_init, const double *const ub_init, const double *const lb_iter, const double *const ub_iter, double (*objective)(const double *const solution, const int D), void (*callback)(const double GlobalMin, const double *const GlobalParams), const int max_function_evaluations, const int NP, const int FoodNumber, const double limit) {
  // ------------------------------------------------------------
  // initialization ---------------------------------------------
  // ------------------------------------------------------------

  g_D = D;
  g_lb = lb_init;
  g_ub = ub_init;
  g_objective = objective;
  g_NP = NP;
  g_FoodNumber = FoodNumber;
  g_limit = limit;

  g_Foods = (double **)malloc(sizeof(double *) * FoodNumber);
  for (int i = 0; i < FoodNumber; i++) {
    g_Foods[i] = (double *)malloc(sizeof(double) * D);
  }
  g_f = (double *)malloc(sizeof(double) * FoodNumber);
  g_fitness = (double *)malloc(sizeof(double) * FoodNumber);
  g_trial = (double *)malloc(sizeof(double) * FoodNumber);
  g_prob = (double *)malloc(sizeof(double) * FoodNumber);
  g_solution = (double *)malloc(sizeof(double) * D);
  g_GlobalParams = (double *)malloc(sizeof(double) * D);

  initial();
  int function_evaluation = FoodNumber;
  MemorizeBestSource();

  // ------------------------------------------------------------
  // iteration --------------------------------------------------
  // ------------------------------------------------------------

  g_lb = lb_iter;
  g_ub = ub_iter;
  while (function_evaluation < max_function_evaluations) {
    SendEmployedBees();
    function_evaluation += FoodNumber;
    CalculateProbabilities();
    SendOnlookerBees();
    function_evaluation += FoodNumber;
    MemorizeBestSource();
    SendScoutBees();
  }

  if (callback) {
    callback(g_GlobalMin, g_GlobalParams);
  }

  for (int i = 0; i < FoodNumber; i++) {
    free(g_Foods[i]);
  }
  free(g_Foods);
  free(g_f);
  free(g_fitness);
  free(g_trial);
  free(g_prob);
  free(g_solution);
  free(g_GlobalParams);
}

/* Fitness function */
static double CalculateFitness(double fun) {
  double result = 0;
  if (fun >= 0) {
    result = 1 / (fun + 1);
  } else {
    result = 1 + fabs(fun);
  }
  return result;
}

/* The best food source is memorized */
static void MemorizeBestSource() {
  for (int i = 0; i < g_FoodNumber; i++) {
    if (g_f[i] < g_GlobalMin) {
      g_GlobalMin = g_f[i];
      for (int j = 0; j < g_D; j++)
        g_GlobalParams[j] = g_Foods[i][j];
    }
  }
}

/* Variables and counters are initialized */
static void init(int index) {
  for (int j = 0; j < g_D; j++) {
    const double r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    g_Foods[index][j] = r * (g_ub[j] - g_lb[j]) + g_lb[j];
    g_solution[j] = g_Foods[index][j];
  }
  g_f[index] = g_objective(g_solution, g_D);
  g_fitness[index] = CalculateFitness(g_f[index]);
  g_trial[index] = 0;
}

/* All food sources are initialized */
static void initial() {
  for (int i = 0; i < g_FoodNumber; i++) {
    init(i);
  }
  g_GlobalMin = g_f[0];
  for (int i = 0; i < g_D; i++)
    g_GlobalParams[i] = g_Foods[0][i];
}

static void SendEmployedBees() {
  /* Employed Bee Phase */
  for (int i = 0; i < g_FoodNumber; i++) {
    /* The parameter to be changed is determined randomly */
    double r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    const int param2change = (int)(r * g_D);

    /* A randomly chosen solution is used in producing a mutant solution of the solution i */
    r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    int neighbor = (int)(r * g_FoodNumber);

    /* Randomly selected solution must be different from the solution i */
    while (neighbor == i) {
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      neighbor = (int)(r * g_FoodNumber);
    }
    for (int j = 0; j < g_D; j++)
      g_solution[j] = g_Foods[i][j];

    /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
    r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    g_solution[param2change] = g_Foods[i][param2change] + (g_Foods[i][param2change] - g_Foods[neighbor][param2change]) * (r - 0.5) * 2;

    /* if generated parameter value is out of boundaries, it is shifted onto the boundaries */
    for (int j = 0; j < g_D; j++) {
      if (g_solution[param2change] < g_lb[j])
        g_solution[param2change] = g_lb[j];
      if (g_solution[param2change] > g_ub[j])
        g_solution[param2change] = g_ub[j];
    }
    const double ObjValSol = g_objective(g_solution, g_D);
    const double FitnessSol = CalculateFitness(ObjValSol);

    /* a greedy selection is applied between the current solution i and its mutant */
    if (FitnessSol > g_fitness[i]) {
      /* If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the g_trial counter of solution i */
      g_trial[i] = 0;
      for (int j = 0; j < g_D; j++)
        g_Foods[i][j] = g_solution[j];
      g_f[i] = ObjValSol;
      g_fitness[i] = FitnessSol;
    } else { /* if the solution i cannot be improved, increase its g_trial counter */
      g_trial[i] = g_trial[i] + 1;
    }
  }

  /* end of employed bee phase */
}

/* A food source is chosen with the probability which is proportional to its quality */
static void CalculateProbabilities() {
  double maxfit = g_fitness[0];
  for (int i = 1; i < g_FoodNumber; i++) {
    if (g_fitness[i] > maxfit)
      maxfit = g_fitness[i];
  }

  for (int i = 0; i < g_FoodNumber; i++) {
    g_prob[i] = (0.9 * (g_fitness[i] / maxfit)) + 0.1;
  }
}

static void SendOnlookerBees() {
  int i = 0;
  int t = 0;
  /* onlooker Bee Phase */
  while (t < g_FoodNumber) {
    double r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
    /* choose a food source depending on its probability to be chosen */
    if (r < g_prob[i]) {
      t++;

      /* The parameter to be changed is determined randomly */
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      const int param2change = (int)(r * g_D);

      /* A randomly chosen solution is used in producing a mutant solution of the solution i */
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      int neighbor = (int)(r * g_FoodNumber);

      /* Randomly selected solution must be different from the solution i */
      while (neighbor == i) {
        r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
        neighbor = (int)(r * g_FoodNumber);
      }
      for (int j = 0; j < g_D; j++)
        g_solution[j] = g_Foods[i][j];

      /* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
      r = ((double)rand() / ((double)(RAND_MAX) + (double)(1)));
      g_solution[param2change] = g_Foods[i][param2change] + (g_Foods[i][param2change] - g_Foods[neighbor][param2change]) * (r - 0.5) * 2;

      /* if generated parameter value is out of boundaries, it is shifted onto the boundaries */
      for (int j = 0; j < g_D; j++) {
        if (g_solution[param2change] < g_lb[j])
          g_solution[param2change] = g_lb[j];
        if (g_solution[param2change] > g_ub[j])
          g_solution[param2change] = g_ub[j];
      }
      const double ObjValSol = g_objective(g_solution, g_D);
      const double FitnessSol = CalculateFitness(ObjValSol);

      /* a greedy selection is applied between the current solution i and its mutant */
      if (FitnessSol > g_fitness[i]) {
        /* If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the g_trial counter of solution i */
        g_trial[i] = 0;
        for (int j = 0; j < g_D; j++)
          g_Foods[i][j] = g_solution[j];
        g_f[i] = ObjValSol;
        g_fitness[i] = FitnessSol;
      } else { /* if the solution i cannot be improved, increase its g_trial counter */
        g_trial[i] = g_trial[i] + 1;
      }
    } /* if */
    i++;
    if (i == g_FoodNumber)
      i = 0;
  } /* while */

  /* end of onlooker bee phase */
}

/* determine the food sources whose g_trial counter exceeds the "g_limit" value. */
static void SendScoutBees() {
  int maxtrialindex = 0;
  for (int i = 1; i < g_FoodNumber; i++) {
    if (g_trial[i] > g_trial[maxtrialindex])
      maxtrialindex = i;
  }
  if (g_trial[maxtrialindex] >= g_limit) {
    init(maxtrialindex);
  }
}

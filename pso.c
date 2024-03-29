#include "pso.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdlib.h>


static inline double uniform() {
    return ((double)rand())/ RAND_MAX;
}

static inline double rd(double lo, double hi) {
	return  uniform() * (hi - lo) + lo;
}

void settings(pso_options *opt, int nvars) {
	opt->dim = nvars;
	opt->c1 = 1.49;
	opt->c2 = 1.49;
	opt->MaxIter = nvars * 100;
	opt->MinFractionNeighbors = 0.25;
	opt->SwarmSize = MAX(300, 20 * nvars);
	opt->StallIterLimit = 20;
	opt->TolFun = 1e-6;
	opt->w_max = 1.1;
	opt->w_min = 0.1;
}

void pso(objFcn fun, void *params, pso_result *solution, pso_options *options) {
	int num = options->SwarmSize, nvars = options->dim;
	double *lb = options->lb, *ub = options->ub;
	double c1 = options->c1, c2 = options->c2;
	double w = options->w_max; //inertia
	int minNeighbors = MAX(2, floor(num * options->MinFractionNeighbors));

	double x[num][nvars]; // position matrix
	double v[num][nvars]; // velocity matrix
	double x_b[num][nvars]; // best position matrix
	double fit[num]; // particle fitness vector
	double fit_b[num]; // best fitness vector

	int neighbors = minNeighbors;
	double bestFval = DBL_MAX, newBest = DBL_MAX;
	double bestFvalsWindow[options->StallIterLimit];
	int exitFlag = 0, iter = 0, w_count = 0;
	int bp;

	int phs[num]; // the prediction horizon particles use
	
	for (int p = 0; p < num; p++) {
		for (int d = 0; d < nvars; d++) {
			x[p][d] = rd(lb[d], ub[d]);
			x_b[p][d] = x[p][d];
			v[p][d] = rd(lb[d] - ub[d], ub[d] - lb[d]);
		}
		fit[p] = fun(x[p], nvars, params, &(phs[p]));
		fit_b[p] = fit[p];
	}

	while (exitFlag == 0) {
		iter++;
		for (int p = 0; p < num; p++) {
			int neighIdx[neighbors - 1];
			int i = 0;
			for (int count = 0; count < neighbors - 1; i++) {
				if (i != p) {
					neighIdx[count] = i;
					count++;
				}
			}
			for (; i < num; i++) {
				if (i != p) {
					int j = rand() % i;
					if (j < neighbors - 1)
						neighIdx[j] = i;
				}
			}
			int neighBest = p;
			for (int k = 0; k < neighbors - 1; k++) {
				if (fit_b[neighIdx[k]] < fit_b[neighBest])
					neighBest = neighIdx[k];
			}

			for (int d = 0; d < nvars; d++) {
				v[p][d] = w * v[p][d]
						+ c1 * uniform() * (x_b[p][d] - x[p][d])
						+ c2 * uniform()
								* (x_b[neighBest][d] - x[p][d]);
				x[p][d] += v[p][d];
				if (x[p][d] < lb[d]) {
					x[p][d] = lb[d];
					v[p][d] = 0;
				}
				if (x[p][d] > ub[d]) {
					x[p][d] = ub[d];
					v[p][d] = 0;
				}
			}

			fit[p] = fun(x[p], nvars, params, &(phs[p]));
			// update personal best position
			if (fit[p] < fit_b[p]) {
				fit_b[p] = fit[p];
				memmove((void *) &x_b[p], (void *) &x[p],
						sizeof(double) * nvars);
			}
		}
		newBest = DBL_MAX;
		for (int p = 0; p < num; p++) {
			if (fit_b[p] < newBest)
				newBest = fit_b[p];
		}
		bestFvalsWindow[(iter - 1) % (options->StallIterLimit)] = newBest;
		if (newBest < bestFval) {
			bestFval = newBest;
			w_count = MAX(0, w_count - 1);
			neighbors = minNeighbors;
		} else {
			w_count++;
			neighbors = MIN(num, neighbors + minNeighbors);
		}
		if (w_count < 2)
			w = MAX(options->w_min, MIN(options->w_max, 2 * w));
		else if (w_count > 5)
			w = MAX(options->w_min, MIN(options->w_max, 0.5 * w));
		exitFlag = iter > options->MaxIter ? 1 : 0;
	}

	solution->error = DBL_MAX;
	for (int p = 0; p < num; p++) {
		if (fit_b[p] < solution->error) {
			solution->error = fit_b[p];
			bp = p;
			//solution->ph = phs[p];
		}
	}
	
	for (int d = 0; d < nvars; d++) {
		solution->gbest[d] = 0;
	}
	int all0_ph = 0;
	double all0_fit = fun(solution->gbest, nvars, params, &all0_ph);
	if (solution->error < all0_fit) {
		memmove((void *) solution->gbest, (void *) &x[bp], sizeof(double) * nvars);
	} else {
		solution->ph = 0;
		solution->error = all0_fit;
	}
}

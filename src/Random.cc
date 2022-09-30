/*
  This file is part of the FRED system.

  Copyright (c) 2010-2015, University of Pittsburgh, John Grefenstette,
  Shawn Brown, Roni Rosenfield, Alona Fyshe, David Galloway, Nathan
  Stone, Jay DePasse, Anuroop Sriram, and Donald Burke.

  Licensed under the BSD 3-Clause license.  See the file "LICENSE" for
  more information.
*/

//
//
// File: Random.cc
//
#include "Random.h"
#include <cmath>
#include <stdio.h>

#define repeat for(;;)

Thread_RNG Random::Random_Number_Generator;

Thread_RNG::Thread_RNG() {
  thread_rng = new RNG [fred::omp_get_max_threads()];
}

void Thread_RNG::set_seed(unsigned long metaseed) {
  std::mt19937_64 seed_generator;
  seed_generator.seed(metaseed);
  for(int t = 0; t < fred::omp_get_max_threads(); ++t) {
    unsigned long new_seed = seed_generator();
    thread_rng[t].set_seed(new_seed);
  }
}


void RNG::set_seed(unsigned long seed) {
  mt_engine.seed(seed);
}

int RNG::draw_from_distribution(int n, double* dist) {
  double r = random();
  int i = 0;
  while(i <= n && dist[i] < r) {
    i++;
  }
  if(i <= n) {
    return i;
  } else {
    printf("Help! draw from distribution failed.\n");
    printf("Is distribution properly formed? (should end with 1.0)\n");
    for(int i = 0; i <= n; i++) {
      printf("%f ", dist[i]);
    }
    printf("\n");
    return -1;
  }
}

double RNG::exponential(double lambda) {
  double u = random();
  return (-log(u) / lambda);
}

// The Gamma sampler is copied from the C++ code underlying rgamma https://github.com/wch/r-source/blob/trunk/src/nmath/rgamma.c
double RNG::gamma(double k, double theta) {
  double a = k;
  double scale = theta;
  /* Constants : */
  const static double sqrt32 = 5.656854;
  const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

  /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
   * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
   * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
   */
  const static double q1 = 0.04166669;
  const static double q2 = 0.02083148;
  const static double q3 = 0.00801191;
  const static double q4 = 0.00144121;
  const static double q5 = -7.388e-5;
  const static double q6 = 2.4511e-4;
  const static double q7 = 2.424e-4;

  const static double a1 = 0.3333333;
  const static double a2 = -0.250003;
  const static double a3 = 0.2000062;
  const static double a4 = -0.1662921;
  const static double a5 = 0.1423657;
  const static double a6 = -0.1367177;
  const static double a7 = 0.1233795;

  /* State variables [FIXME for threading!] :*/
  static double aa = 0.;
  static double aaa = 0.;
  static double s, s2, d;    /* no. 1 (step 1) */
  static double q0, b, si, c;/* no. 2 (step 4) */

  double e, p, q, r, t, u, v, w, x, ret_val;

  if (std::isnan(a) || std::isnan(scale))
    FRED_VERBOSE(1,"In Gamma function, but have NAN parameters");
  if (a <= 0.0 || scale <= 0.0) {
    if(scale == 0. || a == 0.) return 0.;
    FRED_VERBOSE(1,"In Gamma function, but have negative or zero parameters");
  }
  if(!std::isfinite(a) || !std::isfinite(scale)) return INFINITY;

  if (a < 1.) { /* GS algorithm for parameters a < 1 */
    e = 1.0 + exp_m1 * a;
    repeat {
      p = e * random();
      if (p >= 1.0) {
	x = -log((e - p) / a);
	if (exponential(1.0) >= (1.0 - a) * log(x))
	  break;
      } else {
	x = exp(log(p) / a);
	if (exponential(1.0) >= x)
	  break;
      }
    }
    return scale * x;
  }

  /* --- a >= 1 : GD algorithm --- */

  /* Step 1: Recalculations of s2, s, d if a has changed */
  if (a != aa) {
    aa = a;
    s2 = a - 0.5;
    s = sqrt(s2);
    d = sqrt32 - s * 12.0;
  }
  /* Step 2: t = standard normal deviate,
     x = (s,1/2) -normal deviate. */

  /* immediate acceptance (i) */
  t = normal(0.0,1.0);
  x = s + 0.5 * t;
  ret_val = x * x;
  if (t >= 0.0)
    return scale * ret_val;

  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
  u = random();
  if (d * u <= t * t * t)
    return scale * ret_val;

  /* Step 4: recalculations of q0, b, si, c if necessary */

  if (a != aaa) {
    aaa = a;
    r = 1.0 / a;
    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
	   + q2) * r + q1) * r;

    /* Approximation depending on size of parameter a */
    /* The constants in the expressions for b, si and c */
    /* were established by numerical experiments */

    if (a <= 3.686) {
      b = 0.463 + s + 0.178 * s2;
      si = 1.235;
      c = 0.195 / s - 0.079 + 0.16 * s;
    } else if (a <= 13.022) {
      b = 1.654 + 0.0076 * s2;
      si = 1.68 / s + 0.275;
      c = 0.062 / s + 0.024;
    } else {
      b = 1.77;
      si = 0.75;
      c = 0.1515 / s;
    }
  }
  /* Step 5: no quotient test if x not positive */

  if (x > 0.0) {
    /* Step 6: calculation of v and quotient q */
    v = t / (s + s);
    if (fabs(v) <= 0.25)
      q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
				+ a3) * v + a2) * v + a1) * v;
    else
      q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);


    /* Step 7: quotient acceptance (q) */
    if (log(1.0 - u) <= q)
      return scale * ret_val;
  }

  repeat {
    /* Step 8: e = standard exponential deviate
     *	u =  0,1 -uniform deviate
     *	t = (b,si)-double exponential (laplace) sample */
    e = exponential(1.0);
    u = random();
    u = u + u - 1.0;
    if (u < 0.0)
      t = b - si * e;
    else
      t = b + si * e;
    /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
    if (t >= -0.71874483771719) {
      /* Step 10:	 calculation of v and quotient q */
      v = t / (s + s);
      if (fabs(v) <= 0.25)
	q = q0 + 0.5 * t * t *
	  ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
	    + a2) * v + a1) * v;
      else
	q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
      /* Step 11:	 hat acceptance (h) */
      /* (if q not positive go to step 8) */
      if (q > 0.0) {
	w = expm1(q);
	/*  ^^^^^ original code had approximation with rel.err < 2e-7 */
	/* if t is rejected sample again at step 8 */
	if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
	  break;
      }
    }
  } /* repeat .. until  `t' is accepted */
  x = s + 0.5 * t;
  return scale * x * x;
}

int RNG::poisson(double lambda) {
  const int STEP = 500; // This parameter prevents underflow
  double p = 1, lambdaLeft = lambda;
  int k = 0;
  do
    {
      k += 1;
      p *= random();
      while (p < 1 && lambdaLeft > 0) {
	if (lambdaLeft > STEP) {
	  p *= std::exp(STEP);
	  lambdaLeft -= STEP;
	} else {
	  p *= std::exp(lambdaLeft);
	  lambdaLeft = 0;
	}
      }
    } while (p > 1);
  return k-1;
}

int RNG::negative_binomial(double mu, double r) {
  return (mu==0) ? 0 : poisson(gamma(r,mu/r));
}

int RNG::binomial(int n, double p) {
  std::vector<double> cdf;
  build_binomial_cdf(p,n,cdf);
  return draw_from_cdf_vector(cdf);
}

double RNG::normal(double mu, double sigma) {
  return mu + sigma * normal_dist(mt_engine);
}

double RNG::lognormal(double mu, double sigma) {
  double z = normal(0.0,1.0);
  return exp(mu + sigma * z);
}


int RNG::draw_from_cdf(double* v, int size) {
  double r = random();
  int top = size - 1;
  int bottom = 0;
  int s = top / 2;
  while(bottom <= top) {
    if(r <= v[s]) {
      if(s == 0 || r > v[s - 1]) {
        return s;
      } else {
        top = s - 1;
      }
    } else { // r > v[s]
      if(s == size - 1) {
        return s;
      }
      if(r < v[s + 1]) {
        return s + 1;
      } else {
        bottom = s + 1;
      }
    }
    s = bottom + (top - bottom)/2;
  }
  // assert(bottom <= top);
  return -1;
}

int RNG::draw_from_cdf_vector(const vector<double>& v) {
  int size = v.size();
  double r = random();
  int top = size - 1;
  int bottom = 0;
  int s = top / 2;
  while(bottom <= top) {
    if(r <= v[s]) {
      if(s == 0 || r > v[s - 1]) {
        return s;
      } else {
        top = s - 1;
      }
    } else { // r > v[s]
      if(s == size - 1) {
        return s;
      }
      if(r < v[s + 1]) {
        return s + 1;
      } else {
        bottom = s + 1;
      }
    }
    s = bottom + (top - bottom) / 2;
  }
  return -1;
}

double binomial_coefficient(int n, int k) {
  if(k < 0 ||  k > n) {
    return 0;
  }
  if(k > n - k)  {
    k = n - k;
  }
  double c = 1.0;
  for(int i = 0; i < k; ++i) {
    c = c * (n - (k - (i + 1)));
    c = c / (i + 1);
  }
  return c;
}

void RNG::build_binomial_cdf(double p, int n, std::vector<double> &cdf) {
  for(int i = 0; i <= n; ++i) {
    double prob = 0.0;
    for(int j = 0; j <= i; ++j) {
      prob += binomial_coefficient(n, i)
        * pow(10, ((i * log10(p)) + ((n - 1) * log10(1 - p))));
    }
    if(i > 0) {
      prob += cdf.back();
    }
    if(prob < 1) {
      cdf.push_back(prob);
    } else {
      cdf.push_back(1.0);
      break;
    }
  }
  cdf.back() = 1.0;
}

void RNG::build_lognormal_cdf(double mu, double sigma, std::vector<double> &cdf) {
  int maxval = -1;
  int count[1000];
  for(int i = 0; i < 1000; i++) {
    count[i] = 0;
  }
  for(int i = 0; i < 1000; i++) {
    double x = lognormal(mu, sigma);
    int j = (int) x + 0.5;
    if(j > 999) {
      j = 999;
    }
    count[j]++;
    if(j > maxval) {
      maxval = j;
    }
  }
  for(int i = 0; i <= maxval; ++i) {
    double prob = (double) count[i] / 1000.0;
    if(i > 0) {
      prob += cdf.back();
    }
    if(prob < 1.0) {
      cdf.push_back( prob );
    } else {
      cdf.push_back(1.0);
      break;
    }
  }
  cdf.back() = 1.0;
}

void RNG::sample_range_without_replacement(int N, int s, int* result) {
  std::vector<bool> selected(N, false);
  for(int n = 0; n < s; ++n) {
    int i = random_int(0, N - 1);
    if(selected[i]) {
      if(i < N - 1 && !(selected[i + 1])) {
        ++i;
      } else if(i > 0 && !(selected[i - 1])) {
        --i;
      } else {
        --n;
        continue;
      }
    }
    selected[i] = true;
    result[n] = i;
  }
}



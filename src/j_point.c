/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */

#include "boolean.h"
#include "helper.h"
#include "point.h"
#include "j_point.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <fcntl.h>
#include <time.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
#include <gmp.h>

/** Constructor & Destructor */
J_Point init_j_point(J_Point p) {
	mpz_init(p.X);
	mpz_init(p.Y);
	mpz_init(p.Z);
	p.isInf = false;
	return p;
}

J_Point clean_j_point(J_Point p) {
	mpz_clear(p.X);
	mpz_clear(p.Y);
	mpz_clear(p.Z);
	return p;
}

J_Point copy_j_point(J_Point p) {
	J_Point results;
	results = init_j_point(results);
	mpz_set(results.X, p.X);
	mpz_set(results.Y, p.Y);
	mpz_set(results.Z, p.Z);
	results.isInf = p.isInf;
	return p;
}

// Assumption Z is same
int compare_j_point(J_Point p, J_Point q) {
	if ((mpz_cmp(p.X, q.X) == 0) && (mpz_cmp(p.Y, q.Y) == 0)) return 0;
	else return 1;
}

int compare_j_point_negate(J_Point p, J_Point q) {
	mpz_neg(p.Y, q.Y);
	return compare_j_point(p, q);
}

J_Point jacobian_curve_addition(J_Point p, J_Point q, mpz_t a, mpz_t modulo) {
	J_Point results;

	return results;
}

J_Point jacobian_curve_doubling(J_Point p, mpz_t a, mpz_t modulo) {
	J_Point results;

	return results;
}

J_Point jacobian_curve_substraction(J_Point p, J_Point q, mpz_t a, mpz_t modulo) {
	J_Point results;

	return results;
}

/** Convert Affine <-> Jacobian */
J_Point affine_to_jacobian(Point p) {
	J_Point results;

	return results;
}

Point jacobian_to_affine(J_Point p) {
	Point results;

	return results;
}

/* Created by freedomofkeima - 2014 */

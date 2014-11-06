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

// Jacobian coordinate with c = 2 and d = 3

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

// Assumption: p.Z is same with q.Z
int compare_j_point(J_Point p, J_Point q) {
	if ((mpz_cmp(p.X, q.X) == 0) && (mpz_cmp(p.Y, q.Y) == 0)) return 0;
	else return 1;
}

int compare_j_point_negate(J_Point p, J_Point q) {
	mpz_neg(p.Y, q.Y);
	return compare_j_point(p, q);
}


J_Point jacobian_curve_addition(J_Point p, J_Point q, mpz_t a, mpz_t modulo) {
	// Case INF
	if (p.isInf) return q;
	if (q.isInf) return p;

	J_Point results;
	results = init_j_point(results);
	mpz_t one, two, U1, U2, S1, S2;
	mpz_init(one); mpz_init(two);
	mpz_init(U1); mpz_init(U2);
	mpz_init(S1); mpz_init(S2);
	mpz_set_ui(one, 1);
	mpz_set_ui(two, 2);

	if (mpz_cmp(q.Z, one) == 0) { // case Z = 1
		mpz_set(U1, p.X);
		mpz_set(U2, q.X);
		mpz_set(S1, p.Y);
		mpz_set(S2, q.Y);
	} else {
		mpz_t  three;
		mpz_t Z1, Z2; // precompute
		mpz_init(three);
		mpz_init(Z1); mpz_init(Z2);

		mpz_set_ui(three, 3);

		mpz_mul(Z2, q.Z, q.Z);
		mpz_mul(U1, Z2, p.X);

		mpz_mul(Z1, p.Z, p.Z);
		mpz_mul(U2, Z1, q.X);

		mpz_mul(S1, Z2, q.Z);
		mpz_mul(S1, S1, p.Y);

		mpz_mul(S2, Z1, p.Z);
		mpz_mul(S2, S2, q.Y);
	}

	if (mpz_cmp(U1, U2) == 0) {
		if (mpz_cmp(S1, S2) == 0) {
			return jacobian_curve_doubling(p, a, modulo);
		} else {
			results.isInf = true;
			return results;
		}
	}

	mpz_t H, H2, H3, R, R2, temp;
	mpz_init(H); mpz_init(H2); mpz_init(H3);
	mpz_init(R); mpz_init(R2); mpz_init(temp);

	mpz_sub(H, U2, U1);
	positive_modulo(H, H, modulo);
	mpz_sub(R, S2, S1);
	positive_modulo(R, R, modulo);

	// pre-compute
	mpz_mul(H2, H, H); // H^2
	positive_modulo(H2, H2, modulo);
	mpz_mul(H3, H2, H); // H^3
	positive_modulo(H3, H3, modulo);
	mpz_mul(R2, R, R); // R^2
	positive_modulo(R2, R2, modulo);

	mpz_mul(results.X, two, U1); // 2 * U1
	mpz_mul(results.X, results.X, H2); // 2 * U1 * H^2
	positive_modulo(results.X, results.X, modulo);
	mpz_add(results.X, H3, results.X); // H^3 + 2 * U1 * H2
	mpz_sub(results.X, R2, results.X); // R^2 - (H^3 + 2 * U1 * H2)
	positive_modulo(results.X, results.X, modulo);

	mpz_mul(temp, U1, H2); // U1 * H^2
	positive_modulo(temp, temp, modulo);
	mpz_sub(temp, temp, results.X); // U1 * H^2 - X3
	positive_modulo(temp, temp, modulo);
	mpz_mul(temp, temp, R); // R * (U1 * H^2 - X3)
	positive_modulo(temp, temp, modulo);
	mpz_mul(results.Y, S1, H3); // S1 * H^3
	positive_modulo(results.Y, results.Y, modulo);
	mpz_sub(results.Y, temp, results.Y); // R * (U1 * H^2 - X3) - S1 * H^3
	positive_modulo(results.Y, results.Y, modulo);

	if (mpz_cmp(q.Z, one) != 0) { // case Z != 1
		mpz_mul(results.Z, H, p.Z); // H * Z1
		mpz_mul(results.Z, results.Z, q.Z); // H * Z1 * Z2
		positive_modulo(results.Z, results.Z, modulo);
	} else mpz_set(results.Z, H); // H

	return results;
}

J_Point jacobian_curve_doubling(J_Point p, mpz_t a, mpz_t modulo) { // Point doubling
	J_Point results;
	results = init_j_point(results);

	// Case p is a point in infinity
	if (p.isInf) return p;

	// Case y = 0
	mpz_t zero_value;
	mpz_init(zero_value);
	if (mpz_cmp(p.Y, zero_value) == 0) {
		results.isInf = true;
		return results;
	}

	// Check whether a = -3
	mpz_t check, one, two, three, four, eight, X2, Y2, Y4, Z2, Z4;
	mpz_init(check); 
	mpz_init(one); mpz_init(two);
	mpz_init(three); mpz_init(four); mpz_init(eight);
	mpz_init(X2); 
	mpz_init(Y2); mpz_init(Y4); 
	mpz_init(Z2); mpz_init(Z4);
	mpz_set_si(check, -3);
	mpz_set_ui(one, 1);
	mpz_set_ui(two, 2);
	mpz_set_ui(three, 3);
	mpz_set_ui(four, 4);
	mpz_set_ui(eight, 8);
	mpz_set(Z2, p.Z); mpz_set(Z4, p.Z);

	// Pre-compute
	mpz_mul(X2, p.X, p.X);
	positive_modulo(X2, X2, modulo);
	mpz_mul(Y2, p.Y, p.Y);
	positive_modulo(Y2, Y2, modulo);
	mpz_mul(Y4, Y2, Y2);
	positive_modulo(Y4, Y4, modulo);
	if (mpz_cmp(Z4, one) != 0) { // Z^4
		mpz_mul(Z2, p.Z, p.Z);
		positive_modulo(Z2, Z2, modulo);
		if (mpz_cmp(a, check) != 0) {
			mpz_mul(Z4, Z2, Z2);
			positive_modulo(Z4, Z4, modulo);
		}
	}

	// TODO: Check both methods
	mpz_t S, M, temp, temp2;
	mpz_init(S); mpz_init(M); 
	mpz_init(temp); mpz_init(temp2);

	mpz_mul(S, four, p.X); // 4 * X
	mpz_mul(S, S, Y2); // 4 * X * Y^2
	positive_modulo(S, S, modulo);

	if (mpz_cmp(a, check) == 0) { 
		mpz_add(temp, p.X, Z2);
		mpz_sub(temp2, p.X, Z2);
		mpz_mul(M, temp, temp2);
		mpz_mul(M, three, M); // 3 * (X + Z^2) * (X - Z^2)
	} else {
		mpz_mul(X2, three, X2); // 3 * X^2
		mpz_mul(M, a, Z4); // a * Z^4
		mpz_add(M, X2, M); // 3 * X^2 + a * Z^4
	}
	positive_modulo(M, M, modulo);
	
	mpz_mul(results.X, two, S); // 2 * S
	mpz_mul(temp, M, M);
	positive_modulo(temp, temp, modulo);
	mpz_sub(results.X, temp, results.X); // M^2 - 2 * S
	positive_modulo(results.X, results.X, modulo);

	mpz_sub(temp, S, results.X); // S - X'
	mpz_mul(temp, M, temp); // M * (S - X')
	positive_modulo(temp, temp, modulo);
	mpz_mul(results.Y, eight, Y4); // 8 * Y^4
	mpz_sub(results.Y, temp, results.Y); // M * (S - X') - 8 * Y^4
	positive_modulo(results.Y, results.Y, modulo);

	mpz_mul(results.Z, two, p.Y); // 2 * Y
	mpz_mul(results.Z, results.Z, p.Z); // 2 * Y * Z
	positive_modulo(results.Z, results.Z, modulo);

	return results;
}

J_Point jacobian_curve_substraction(J_Point p, J_Point q, mpz_t a, mpz_t modulo) {
	mpz_neg(q.Y, q.Y);
	return jacobian_curve_addition(p, q, a, modulo);
}

/** Convert Affine <-> Jacobian */
J_Point affine_to_jacobian(Point p) {
	J_Point results;
	results = init_j_point(results);
	mpz_set(results.X, p.x);
	mpz_set(results.Y, p.y);
	mpz_set_ui(results.Z, 1);
	results.isInf = p.isInf;
	return results;
}

Point jacobian_to_affine(J_Point p, mpz_t modulo) {
	Point results;
	results = init_point(results);

	mpz_t Z, one;
	mpz_init(Z);
	mpz_set(Z, p.Z);

	mpz_mul(Z, Z, p.Z); // Z^2
	mpz_invert(results.x, Z, modulo); // Z^-2
	mpz_mul(results.x, results.x, p.X); // X.Z^-2
	positive_modulo(results.x, results.x, modulo);

	mpz_mul(Z, Z, p.Z); // Z^3
	mpz_invert(results.y, Z, modulo); // Z^-3
	mpz_mul(results.y, results.y, p.Y); // Y.Z^-3
	positive_modulo(results.y, results.y, modulo);

	results.isInf = p.isInf;

	return results;
}

/* Created by freedomofkeima - 2014 */
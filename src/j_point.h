/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */

#ifndef J_POINT_H
#define J_POINT_H

#include "boolean.h"
#include "helper.h"
#include "point.h"

typedef struct {
	mpz_t X;
	mpz_t Y;
	mpz_t Z;
	boolean isInf;
} J_Point; // Jacobian Coordinates

/** Constructor & Destructor */
J_Point init_j_point(J_Point p);
J_Point clean_j_point(J_Point p);
J_Point copy_j_point(J_Point p);

int compare_j_point(J_Point p, J_Point q);
int compare_j_point_negate(J_Point p, J_Point q);

J_Point jacobian_curve_addition(J_Point p, J_Point q, mpz_t a, mpz_t modulo);
J_Point jacobian_curve_doubling(J_Point p, mpz_t a, mpz_t modulo);
J_Point jacobian_curve_substraction(J_Point p, J_Point q, mpz_t a, mpz_t modulo);

/** Convert Affine <-> Jacobian */
J_Point affine_to_jacobian(Point p);
Point jacobian_to_affine(J_Point p, mpz_t modulo);

#endif
/* Created by freedomofkeima - 2014 */

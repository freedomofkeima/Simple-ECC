/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */

#ifndef ECC_H
#define ECC_H

#include "boolean.h"

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

typedef struct {
	mpz_t x;
	mpz_t y;
	boolean isInf;
} Point;

/** Constructor & Destructor */
Point init_point(Point p);
Point clean_point(Point p);
Point copy_point(Point p);

int compare_point(Point p, Point q);
int compare_point_negate(Point p, Point q);

Point affine_curve_addition(Point p, Point q, mpz_t a, mpz_t modulo);
Point affine_curve_doubling(Point p, mpz_t a, mpz_t modulo);
Point affine_curve_substraction(Point p, Point q, mpz_t a, mpz_t modulo);

/** Left-to-right binary algorithm */
Point left_to_right_binary(Point p, mpz_t a, mpz_t k, mpz_t modulo);
/** Right-to-left binary algorithm */
Point right_to_left_binary(Point p, mpz_t a, mpz_t k, mpz_t modulo);
/** Montgomery ladder algorithm */
Point montgomery_ladder(Point p, mpz_t a, mpz_t k, mpz_t modulo);

void positive_modulo(mpz_t results, mpz_t a, mpz_t modulo);

/** Encrypt & Decrypt */
Point encrypt_ECIES(mpz_t encrypted_message, char* message, Point public_key, Point p, mpz_t a, mpz_t modulo); // return chosen point
void decrypt_ECIES(mpz_t encrypted_message, Point chosen_point, mpz_t private_key, Point p, mpz_t a, mpz_t modulo); // return message

#endif
/* Created by freedomofkeima - 2014 */

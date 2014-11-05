/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */

#include "boolean.h"
#include "ecc.h"

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
Point init_point(Point p) {
	mpz_init(p.x);
	mpz_init(p.y);
	p.isInf = false;
	return p;
}

Point clean_point(Point p) {
	mpz_clear(p.x);
	mpz_clear(p.y);
	return p;
}

Point copy_point(Point p) {
	Point results;
	results = init_point(results);
	mpz_set(results.x, p.x);
	mpz_set(results.y, p.y);
	results.isInf = p.isInf;
	return p;
}

int compare_point(Point p, Point q) {
	if ((mpz_cmp(p.x, q.x) == 0) && (mpz_cmp(p.y, q.y) == 0)) return 0;
	else return 1;
}

int compare_point_negate(Point p, Point q) {
	mpz_neg(q.y, q.y);
	return compare_point(p, q);
}

Point affine_curve_addition(Point p, Point q, mpz_t a, mpz_t modulo) {
	// Case INF
	if (p.isInf) return q;
	if (q.isInf) return p;

	Point results;
	results = init_point(results);

	// Case P = -Q
	if (compare_point_negate(p, q) == 0) {
		results.isInf = true;
		return results;
	}
	// Case P = Q
	if (compare_point(p, q) == 0) return affine_curve_doubling(p, a, modulo);

	mpz_t lambda, lh, rh;
	mpz_init(lambda);
	mpz_init(lh); mpz_init(rh);

	/** Create (q.x - p.x)^-1 % m */
	mpz_sub(lh, q.x, p.x);
	mpz_invert(lh, lh, modulo);

	mpz_sub(rh, q.y, p.y);
	/** Calculate lambda (q.y - p.y) . (q.x - p.x)^-1 % m */
	mpz_mul(lambda, lh, rh);
	positive_modulo(lambda, lambda, modulo);

	/** new_x = lambda^2 - p.x - q.x */
	mpz_mul(results.x, lambda, lambda);
	mpz_sub(results.x, results.x, p.x);
	mpz_sub(results.x, results.x, q.x);
	positive_modulo(results.x, results.x, modulo);

	/** new_y = lambda * (p.x - new_x) - p.y */
	mpz_sub(results.y, p.x, results.x);
	mpz_mul(results.y, results.y, lambda);
	mpz_sub(results.y, results.y, p.y);
	positive_modulo(results.y, results.y, modulo);

	return results;
}

Point affine_curve_doubling(Point p, mpz_t a, mpz_t modulo) {
	Point results;
	results = init_point(results);

	// Case y = 0
	mpz_t zero_value;
	mpz_init(zero_value);
	if (mpz_cmp(p.y, zero_value) == 0) {
		results.isInf = true;
		return results;
	}

	mpz_t lambda, mul, lh, rh;
	mpz_init(lambda); mpz_init(mul);
	mpz_init(lh); mpz_init(rh);

	/** Create (2y)^-1 % m */
	mpz_set_str(mul, "2", 10);
	mpz_mul(lh, mul, p.y);
	// Assuming inverse exists (retval = 1)
	mpz_invert(lh, lh, modulo);

	mpz_set_str(mul, "3", 10);
	mpz_mul(rh, p.x, p.x); // x^2
	positive_modulo(rh, rh, modulo); // x^2 % m
	mpz_mul(rh, rh, mul); // (3x^2)
	mpz_add(rh, rh, a); // (3x^2 + a)
	positive_modulo(rh, rh, modulo); // (3x^2 + a) % m

	mpz_mul(lambda, lh, rh);
	positive_modulo(lambda, lambda, modulo); // (3x^2 + a) . (2y)^-1 % m

	/** new_x = lambda^2 - 2 * previous_x */
	mpz_mul(results.x, lambda, lambda);
	mpz_sub(results.x, results.x, p.x);
	mpz_sub(results.x, results.x, p.x);
	positive_modulo(results.x, results.x, modulo);

	/** new_y = lambda * (previous_x - new_x) - previous_y */
	mpz_sub(results.y, p.x, results.x);
	mpz_mul(results.y, results.y, lambda);
	mpz_sub(results.y, results.y, p.y);
	positive_modulo(results.y, results.y, modulo);

	return results;
}

Point affine_curve_substraction(Point p, Point q, mpz_t a, mpz_t modulo) {
	mpz_neg(q.y, q.y);
	return affine_curve_addition(p, q, a, modulo);
}

/** Left-to-right binary algorithm */
Point left_to_right_binary(Point p, mpz_t a, mpz_t k, mpz_t modulo) {
	Point results, temp;
	results = copy_point(p); // R0	
	temp = copy_point(p); // R1

	int i; // loop variable
	char * binary;
	int binary_size = 0;

	binary = mpz_get_str(NULL, 2, k);
	binary_size = (int) strlen(binary);


	for (i = binary_size - 2; i >= 0; i--) { // i = n-2 downto 0
		results = affine_curve_doubling(results, a, modulo);
		if (binary[binary_size - i - 1] == '1') results = affine_curve_addition(results, temp, a, modulo);
	}

	return results;
}

/** Right-to-left binary algorithm */
Point right_to_left_binary(Point p, mpz_t a, mpz_t k, mpz_t modulo) {
	Point results, temp;
	results = init_point(results);
	results.isInf = true;
	temp = copy_point(p);

	int i; // loop variable
	char * binary;
	int binary_size = 0;

	binary = mpz_get_str(NULL, 2, k);
	binary_size = (int) strlen(binary);

	for (i = 0; i < binary_size; i++) { // i = 0 to n-1 
		if(binary[binary_size - i - 1] == '1') results = affine_curve_addition(results, temp, a, modulo);
		temp = affine_curve_doubling(temp, a, modulo);
	}

	return results;
}

/** Montgomery ladder algorithm */
Point montgomery_ladder(Point p, mpz_t a, mpz_t k, mpz_t modulo) {
	Point results, temp;
	results = init_point(results);
	results.isInf = true;
	temp = copy_point(p);

	int i; // loop variable
	char * binary;
	int binary_size = 0;

	binary = mpz_get_str(NULL, 2, k);
	binary_size = (int) strlen(binary);

	for (i = binary_size - 1; i >= 0; i--) { // i = n-1 downto 0
		if (binary[binary_size - i - 1] == '1') {
			results = affine_curve_addition(results, temp, a, modulo);
			temp = affine_curve_doubling(temp, a, modulo);
		} else {
			temp = affine_curve_addition(results, temp, a, modulo);
			results = affine_curve_doubling(results, a, modulo);
		}
	}

	return results;
}

void positive_modulo(mpz_t results, mpz_t a, mpz_t modulo) {
	mpz_t zero_value;
	mpz_init(zero_value);

	mpz_tdiv_r(results, a, modulo);
	if (mpz_cmp(results, zero_value) < 0) mpz_add(results, results, modulo);
}

/** Encrypt & Decrypt */
Point encrypt_ECIES(mpz_t encrypted_message, char* message, Point public_key, Point p, mpz_t a, mpz_t modulo) { // return chosen point
	Point chosen_point, encoded_point;
	mpz_t k_random;
	chosen_point = init_point(chosen_point);
	encoded_point = init_point(encoded_point);
	mpz_init(k_random);
	get_random(k_random, 32); // choose random value k

	mpz_set_str(encrypted_message, message, 36);
	printf("[SIMPLIFIED ECIES] Message: ");
	mpz_out_str(stdout, 36, encrypted_message);
	printf("\n");

	// Encrypt for Simplified ECIES (Q = public_key_2)
	// Compute kP
	chosen_point = left_to_right_binary(p, a, k_random, modulo); // kP
	gmp_printf("[SIMPLIFIED ECIES] Chosen point [X Y]: %Zd %Zd\n", chosen_point.x, chosen_point.y);
	encoded_point = left_to_right_binary(public_key, a, k_random, modulo); // kQ
	gmp_printf("[SIMPLIFIED ECIES] Encoded point [X Y]: %Zd %Zd\n", encoded_point.x, encoded_point.y);
	mpz_mul(encrypted_message, encrypted_message, encoded_point.x);
	positive_modulo(encrypted_message, encrypted_message, modulo);

	return chosen_point;
}

void decrypt_ECIES(mpz_t encrypted_message, Point chosen_point, mpz_t private_key, Point p, mpz_t a, mpz_t modulo) { // return message
	Point decoded_point;
	decoded_point = init_point(decoded_point);

	// Decrypt for Simplified ECIES (using modulo inverse)
	decoded_point = left_to_right_binary(chosen_point, a, private_key, modulo);
	gmp_printf("[SIMPLIFIED ECIES] Decoded point [X Y]: %Zd %Zd\n", decoded_point.x, decoded_point.y);
	mpz_invert(decoded_point.x, decoded_point.x, modulo);
	mpz_mul(encrypted_message, encrypted_message, decoded_point.x);
	positive_modulo(encrypted_message, encrypted_message, modulo);
	printf("[SIMPLIFIED ECIES] Original message: ");
	mpz_out_str(stdout, 36, encrypted_message);
	printf("\n\n");
}

/* Created by freedomofkeima - 2014 */

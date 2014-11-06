/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */

#include "boolean.h"
#include "point.h"
#include "j_point.h"
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

/** Left-to-right binary algorithm */
Point affine_left_to_right_binary(Point p, mpz_t a, mpz_t k, mpz_t modulo) {
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

J_Point jacobian_left_to_right_binary(J_Point p, mpz_t a, mpz_t k, mpz_t modulo) {
	J_Point results, temp;
	results = copy_j_point(p); // R0	
	temp = copy_j_point(p); // R1

	int i; // loop variable
	char * binary;
	int binary_size = 0;

	binary = mpz_get_str(NULL, 2, k);
	binary_size = (int) strlen(binary);


	for (i = binary_size - 2; i >= 0; i--) { // i = n-2 downto 0
		results = jacobian_curve_doubling(results, a, modulo);
		if (binary[binary_size - i - 1] == '1') results = jacobian_curve_addition(results, temp, a, modulo);
	}

	return results;
}

/** Right-to-left binary algorithm */
Point affine_right_to_left_binary(Point p, mpz_t a, mpz_t k, mpz_t modulo) {
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
Point affine_montgomery_ladder(Point p, mpz_t a, mpz_t k, mpz_t modulo) {
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

J_Point jacobian_montgomery_ladder(J_Point p, mpz_t a, mpz_t k, mpz_t modulo) {
	J_Point results, temp;
	results = init_j_point(results);
	results.isInf = true;
	temp = copy_j_point(p);

	int i; // loop variable
	char * binary;
	int binary_size = 0;

	binary = mpz_get_str(NULL, 2, k);
	binary_size = (int) strlen(binary);

	for (i = binary_size - 1; i >= 0; i--) { // i = n-1 downto 0
		if (binary[binary_size - i - 1] == '1') {
			results = jacobian_curve_addition(results, temp, a, modulo);
			temp = jacobian_curve_doubling(temp, a, modulo);
		} else {
			temp = jacobian_curve_addition(results, temp, a, modulo);
			results = jacobian_curve_doubling(results, a, modulo);
		}
	}

	return results;
}

/** Encrypt & Decrypt */
Point encrypt_ECIES(mpz_t encrypted_message, char* message, Point public_key, Point p, mpz_t a, mpz_t modulo) { // return chosen point
	Point chosen_point, encoded_point;
	J_Point j_chosen_point;
	J_Point j_encoded_point;
	J_Point j_p = affine_to_jacobian(p); // G in Jacobian
	J_Point j_public_key = affine_to_jacobian(public_key); // Public key receiver in Jacobian
	mpz_t k_random;
	mpz_init(k_random);
	get_random(k_random, 32); // choose random value k

	mpz_set_str(encrypted_message, message, 36);
	printf("[SIMPLIFIED ECIES] Message: ");
	mpz_out_str(stdout, 36, encrypted_message);
	printf("\n");

	// Encrypt for Simplified ECIES (Q = public_key_2)
	// Compute kP
	j_chosen_point = jacobian_left_to_right_binary(j_p, a, k_random, modulo); // kP
	gmp_printf("[SIMPLIFIED ECIES] Chosen point - Jacobian [X Y Z]: %Zd %Zd %Zd\n", j_chosen_point.X, j_chosen_point.Y, j_chosen_point.Z);
	chosen_point = jacobian_to_affine(j_chosen_point, modulo);
	gmp_printf("[SIMPLIFIED ECIES] Chosen point - Affine [X Y]: %Zd %Zd\n", chosen_point.x, chosen_point.y);
	j_encoded_point = jacobian_left_to_right_binary(j_public_key, a, k_random, modulo); // kQ
	gmp_printf("[SIMPLIFIED ECIES] Encoded point - Jacobian [X Y Z]: %Zd %Zd %Zd\n", j_encoded_point.X, j_encoded_point.Y, j_encoded_point.Z);
	encoded_point = jacobian_to_affine(j_encoded_point, modulo);
	gmp_printf("[SIMPLIFIED ECIES] Encoded point - Affine [X Y]: %Zd %Zd\n", encoded_point.x, encoded_point.y);
	mpz_mul(encrypted_message, encrypted_message, encoded_point.x);
	positive_modulo(encrypted_message, encrypted_message, modulo);

	return chosen_point;
}

void decrypt_ECIES(mpz_t encrypted_message, Point chosen_point, mpz_t private_key, Point p, mpz_t a, mpz_t modulo) { // return message
	Point decoded_point;
	J_Point j_decoded_point;
	J_Point j_chosen_point = affine_to_jacobian(chosen_point);

	// Decrypt for Simplified ECIES (using modulo inverse)
	j_decoded_point = jacobian_left_to_right_binary(j_chosen_point, a, private_key, modulo);
	gmp_printf("[SIMPLIFIED ECIES] Decoded point - Jacobian [X Y Z]: %Zd %Zd %Zd\n", j_decoded_point.X, j_decoded_point.Y, j_decoded_point.Z);
	decoded_point = jacobian_to_affine(j_decoded_point, modulo);
	gmp_printf("[SIMPLIFIED ECIES] Decoded point - Affine [X Y]: %Zd %Zd\n", decoded_point.x, decoded_point.y);
	mpz_invert(decoded_point.x, decoded_point.x, modulo);
	mpz_mul(encrypted_message, encrypted_message, decoded_point.x);
	positive_modulo(encrypted_message, encrypted_message, modulo);
	printf("[SIMPLIFIED ECIES] Original message: ");
	mpz_out_str(stdout, 36, encrypted_message);
	printf("\n\n");
}

/* Created by freedomofkeima - 2014 */

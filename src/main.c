/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */

#include "boolean.h"
#include "helper.h"
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

// ECC Parameters (P-256 NIST)
char*a_v="-3";
char*b_v="5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b";
char*p_v="115792089210356248762697446949407573530086143415290314195533631308867097853951";
char*r_v="115792089210356248762697446949407573529996955224135760342422259061068512044369";
char*gx_v="5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b";
char*gy_v="4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
// What is s and c?
char*s_v="c49d360886e704936a6678e1139d26b7819f7e90";
char*c_v="7efba1662985be9403cb055c75d4f7e0ce8d84a9c5114abcaf3177680104fa0d";

long long max_iteration = 100;

static struct timeval tm1;

static inline void start() {
	gettimeofday(&tm1, NULL);
}

static inline void stop() {
	struct timeval tm2;
	gettimeofday(&tm2, NULL);

	double t = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
	printf("Average elapsed time: %.2f ms\n", (double) t / max_iteration);
}

int main() {
	/** Fm, where modulo = m */
	mpz_t a, b, k, r, modulo;
	int i = 0; // loop variable
	Point p, next_p;
	p = init_point(p);
	mpz_init(a);
	mpz_init(b);
	mpz_init(k);
	mpz_init(r); // order
	mpz_init(modulo);

	/** Initialize parameters of ECC (F2p) */
	mpz_set_str(a, a_v, 10);
	mpz_set_str(b, b_v, 16);
	mpz_set_str(modulo, p_v, 10);
	mpz_set_str(r, r_v, 10);
	mpz_set_str(p.x, gx_v, 16);
	mpz_set_str(p.y, gy_v, 16);

	get_random(k, 32); // generate random test

	// Convert Affine coordinate to Jacobian coordinate
	J_Point j_p, j_next_p;
	j_next_p = init_j_point(j_next_p);
	j_p = affine_to_jacobian(p); // Generator point

	/** Test Left-to-right binary algorithm */
	start(); // start operation
	while (i < max_iteration) {
		next_p = affine_left_to_right_binary(p, a, k, modulo); // Q = [k]P
		// gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
		i++;
	}
	printf("--[AFFINE] Left to right binary algorithm--\n");
	stop(); // stop operation

	start(); // start operation
	i = 0;
	while (i < max_iteration) {
		j_next_p = jacobian_left_to_right_binary(j_p, a, k, modulo); // Q = [k]P
		// gmp_printf("%Zd %Zd\n", j_next_p.X, j_next_p.Y);
		next_p = jacobian_to_affine(j_next_p, modulo);
		// gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
		i++;
	}
	printf("--[JACOBIAN] Left to right binary algorithm--\n");
	stop(); // stop operation

	/** Test Right-to-left binary algorithm */
	start(); // start operation
	i = 0;
	while (i < max_iteration) {
		next_p = affine_right_to_left_binary(p, a, k, modulo); // Q = [k]P
		// gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
		i++;
	}
	printf("--[AFFINE] Right to left binary algorithm--\n");
	stop(); // stop operation

	/** Test Montgomery ladder algorithm (Against time-based attack) */
	start(); // start operation
	i = 0;
	while (i < max_iteration) {
		next_p = affine_montgomery_ladder(p, a, k, modulo); // Q = [k]P
		//gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
		i++;
	}
	printf("--[AFFINE] Montgomery ladder algorithm--\n");
	stop(); // stop operation

	/** -------------------------------------------------------------------------*/
	Point public_key_1, public_key_2, shared_key;
	mpz_t private_key_1, private_key_2;
	mpz_init(private_key_1); mpz_init(private_key_2);
	// TODO : Key should be padded to fixed size (serializable)
	// Note: (2^-256 chance of failure, can be ignored)
	get_random(private_key_1, 32); // 256 bit
	get_random(private_key_2, 32); // 256 bit

	gmp_printf("Private key [A B]: %Zd %Zd\n\n", private_key_1, private_key_2);
	public_key_1 = affine_left_to_right_binary(p, a, private_key_1, modulo);
	public_key_2 = affine_left_to_right_binary(p, a, private_key_2, modulo);

	gmp_printf("Public key 1 [X Y]: %Zd %Zd\n", public_key_1.x, public_key_1.y);
	gmp_printf("Public key 2 [X Y]: %Zd %Zd\n\n", public_key_2.x, public_key_2.y);

	/** -------------------------------------------------------------------------*/
	// ElGamal Encrypt - Decrypt (Map message to chunk of points in EC)
	Point message, chosen_point, encoded_point, decoded_point;
	mpz_t k_message;
	mpz_init(k_message);
	mpz_set_ui(k_message, 123456789);
	message = affine_left_to_right_binary(p, a, k_message, modulo);
	gmp_printf("[Encrypt] Message [X Y]: %Zd %Zd\n", message.x, message.y);
	get_random(k_message, 32);
	// Encrypt example
	chosen_point = affine_left_to_right_binary(p, a, k_message, modulo); // chosen point (r)
	gmp_printf("[Encrypt] Chosen point [X Y]: %Zd %Zd\n", chosen_point.x, chosen_point.y);
	encoded_point = affine_left_to_right_binary(public_key_2, a, k_message, modulo); // r * Pu2
	encoded_point = affine_curve_addition(message, encoded_point, a, modulo);
	// TODO : chosen_point & encoded_point should be padded to P-bit
	gmp_printf("[Decrypt] Encoded point [X Y]: %Zd %Zd\n", encoded_point.x, encoded_point.y);
	
	// Decrypt example (encoded_point - private_key * chosen_point)
	decoded_point = affine_left_to_right_binary(chosen_point, a, private_key_2, modulo);
	decoded_point = affine_curve_substraction(encoded_point, decoded_point, a, modulo);
	gmp_printf("[Decrypt] Original message [X Y]: %Zd %Zd\n\n", decoded_point.x, decoded_point.y);

	/** -------------------------------------------------------------------------*/
	// Simplified ECIES (Ref: Page 256 Cryptography Theory & Practice 2nd Ed. - Douglas)
	char* message_string = "hello"; // 0..9, a..z (base 36)
	mpz_t encrypted_message;
	mpz_init(encrypted_message);
	chosen_point = encrypt_ECIES(encrypted_message, message_string, public_key_2, p, a, modulo);
	gmp_printf("[SIMPLIFIED ECIES] Encrypted message: %Zd\n", encrypted_message);
	decrypt_ECIES(encrypted_message, chosen_point, private_key_2, p, a, modulo);

	/** -------------------------------------------------------------------------*/
	// TODO : Public key validation!
	// Shared key (ECDH) - key secure exchange
	shared_key = affine_left_to_right_binary(public_key_2, a, private_key_1, modulo);
	gmp_printf("Shared key [X Y]: %Zd %Zd\n", shared_key.x, shared_key.y);

	// TODO : ECDSA - digital signature algorithm (Consider + ECDH)

	/** Cleaning up */
	p = clean_point(p);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(k);
	mpz_clear(r);
	mpz_clear(modulo);
	public_key_1 = clean_point(public_key_1);
	public_key_2 = clean_point(public_key_2);
	mpz_clear(private_key_1);
	mpz_clear(private_key_2);

	return EXIT_SUCCESS;
}
/* Created by freedomofkeima - 2014 */

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

// Define your TEST here
#define TEST_MODULAR_OPERATION true
#define TEST_SCALAR_OPERATION true
#define TEST_SCALAR_ALGORITHM true
#define TEST_ENCRYPT_DECRYPT true
#define TEST_SIMPLIFIED_ECIES true


long long max_iteration;

// ECC Parameters (P-256 NIST)
char*a_v="-3";
char*b_v="5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b";
char*p_v="115792089210356248762697446949407573530086143415290314195533631308867097853951";
char*r_v="115792089210356248762697446949407573529996955224135760342422259061068512044369";
char*gx_v="6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
char*gy_v="4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
char*s_v="c49d360886e704936a6678e1139d26b7819f7e90";
char*c_v="7efba1662985be9403cb055c75d4f7e0ce8d84a9c5114abcaf3177680104fa0d";

static struct timeval tm1;

static inline void start() {
	gettimeofday(&tm1, NULL);
}

static inline void stop() {
	struct timeval tm2;
	gettimeofday(&tm2, NULL);

	double t = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
	printf("Average elapsed time: %.8f ms\n", (double) t / max_iteration);
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

	mpz_t zero_value, k2;
	mpz_init(zero_value);
	mpz_init(k2);

	/** Compare ADDITION, MULTIPLICATION, and INVERSION */
	if (TEST_MODULAR_OPERATION) {
		max_iteration = 10000;

		while (mpz_cmp(k, zero_value) == 0) {
			get_random(k, 32); // generate random test (256 bits)
			positive_modulo(k, k, modulo);
		}
		printf("Random k (in Binary): ");
		mpz_out_str(stdout, 2, k);
		printf("\n");

		while (mpz_cmp(k2, zero_value) == 0) {
			get_random(k2, 32); // generate random test (256 bits)
			positive_modulo(k2, k2, modulo);
		}
		printf("Random k2 (in Binary): ");
		mpz_out_str(stdout, 2, k2);
		printf("\n");

		/** Addition */
		start(); // start operation
		i = 0;
		while (i < max_iteration) {
			mpz_add(k, k, k2);
			positive_modulo(k, k, modulo);
			i++;
		}
		printf("--[ADDITION]--\n");
		stop(); // stop operation

		/** Multiplication */
		start(); // start operation
		i = 0;
		mpz_t two;
		mpz_init(two);
		mpz_set_si(two, 2);
		while (i < max_iteration) {
			mpz_mul(k, k, two);
			positive_modulo(k, k, modulo);
			i++;
		}
		printf("--[MULTIPLICATION k * 2]--\n");
		stop(); // stop operation

		start(); // start operation
		i = 0;
		while (i < max_iteration) {
			mpz_mul(k, k, k2);
			positive_modulo(k, k, modulo);
			i++;
		}
		printf("--[MULTIPLICATION k * k]--\n");
		stop(); // stop operation

		/** Inversion */
		start(); // start operation
		i = 0;
		while (i < max_iteration) {
			mpz_invert(k, k, modulo);
			positive_modulo(k, k, modulo);
			i++;
		}
		printf("--[INVERSION]--\n");
		stop(); // stop operation
	}

	/** -------------------------------------------------------------------------*/
	// Convert Affine coordinate to Jacobian coordinate
	J_Point j_p, j_next_p;
	j_next_p = init_j_point(j_next_p);
	j_p = affine_to_jacobian(p); // Generator point

	if (TEST_SCALAR_OPERATION) {
		max_iteration = 100;
		/** Affine addition */

		/** Affine doubling */

		/** Jacobian addition */

		/** Jacobian doubling */

		/** Affine-Jacobian addition */

		/** Affine-Jacobian doubling */
	}

	/** -------------------------------------------------------------------------*/
	if (TEST_SCALAR_ALGORITHM) {
		max_iteration = 100;
		while (mpz_cmp(k, zero_value) == 0) {
			get_random(k, 32); // generate random test (256 bits)
			positive_modulo(k, k, modulo);
		}
		printf("\nRandom k (in Binary): ");
		mpz_out_str(stdout, 2, k);
		printf("\n");

		/** Test Left-to-right binary algorithm */
		start(); // start operation
		i = 0;
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

		start(); // start operation
		i = 0;
		while (i < max_iteration) {
			j_next_p = jacobian_affine_left_to_right_binary(j_p, p, a, k, modulo); // Q = [k]P
			// gmp_printf("%Zd %Zd\n", j_next_p.X, j_next_p.Y);
			next_p = jacobian_to_affine(j_next_p, modulo);
			// gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
			i++;
		}
		printf("--[JACOBIAN-AFFINE] Left to right binary algorithm--\n");
		stop(); // stop operation

		start(); // start operation
		int w = 4; // windows size
		i = 0;
		while (i < max_iteration) {
			j_next_p = jacobian_affine_sliding_NAF(j_p, p, a, k, modulo, w); // Q = [k]P
			// gmp_printf("%Zd %Zd\n", j_next_p.X, j_next_p.Y);
			next_p = jacobian_to_affine(j_next_p, modulo);
			// gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
			i++;
		}
		printf("--[JACOBIAN-AFFINE] Sliding NAF Left to right binary algorithm (w = 4)--\n");
		stop(); // stop operation

		start(); // start operation
		w = 5; // windows size
		i = 0;
		while (i < max_iteration) {
			j_next_p = jacobian_affine_sliding_NAF(j_p, p, a, k, modulo, w); // Q = [k]P
			// gmp_printf("%Zd %Zd\n", j_next_p.X, j_next_p.Y);
			next_p = jacobian_to_affine(j_next_p, modulo);
			// gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
			i++;
		}
		printf("--[JACOBIAN-AFFINE] Sliding NAF Left to right binary algorithm (w = 5)--\n");
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
			j_next_p = jacobian_montgomery_ladder(j_p, a, k, modulo); // Q = [k]P
			// gmp_printf("%Zd %Zd\n", j_next_p.X, j_next_p.Y);
			next_p = jacobian_to_affine(j_next_p, modulo);
			// gmp_printf("%Zd %Zd\n", next_p.x, next_p.y);
			i++;
		}
		printf("--[JACOBIAN] Montgomery ladder algorithm--\n");
		stop(); // stop operation
	}

	/** -------------------------------------------------------------------------*/
	J_Point public_key_1, public_key_2, shared_key;
	mpz_t private_key_1, private_key_2;
	mpz_init(private_key_1); mpz_init(private_key_2);
	// TODO : Key should be padded to fixed size (serializable)
	// Note: (2^-256 chance of failure, can be ignored)
	while (mpz_cmp(private_key_1, zero_value) == 0) {
		get_random(private_key_1, 32); // 256 bit
		positive_modulo(private_key_1, private_key_1, modulo);
	}
	while (mpz_cmp(private_key_2, zero_value) == 0) {
		get_random(private_key_2, 32); // 256 bit
		positive_modulo(private_key_2, private_key_2, modulo);
	}

	gmp_printf("Private key [A B]: %Zd %Zd\n\n", private_key_1, private_key_2);
	public_key_1 = jacobian_left_to_right_binary(j_p, a, private_key_1, modulo);
	public_key_2 = jacobian_left_to_right_binary(j_p, a, private_key_2, modulo);

	gmp_printf("Public key 1 - Jacobian [X Y Z]: %Zd %Zd %Zd\n", public_key_1.X, public_key_1.Y, public_key_1.Z);
	gmp_printf("Public key 2 - Jacobian [X Y Z]: %Zd %Zd %Zd\n", public_key_2.X, public_key_2.Y, public_key_2.Z);

	Point public_key_1_decoded = jacobian_to_affine(public_key_1, modulo);
	Point public_key_2_decoded = jacobian_to_affine(public_key_2, modulo);

	gmp_printf("Public key 1 - Affine [X Y]: %Zd %Zd\n", public_key_1_decoded.x, public_key_1_decoded.y);
	gmp_printf("Public key 2 - Affine [X Y]: %Zd %Zd\n\n", public_key_2_decoded.x, public_key_2_decoded.y);

	/** -------------------------------------------------------------------------*/
	if (TEST_ENCRYPT_DECRYPT) {
		// ElGamal Encrypt - Decrypt (Map message to chunk of points in EC)
		J_Point message, chosen_point, encoded_point, decoded_point;
		mpz_t k_message;
		mpz_init(k_message);
		mpz_set_ui(k_message, 123456789);
		message = jacobian_left_to_right_binary(j_p, a, k_message, modulo);

		Point message_decoded = jacobian_to_affine(message, modulo);
		gmp_printf("[Encrypt] Message - Affine [X Y] %Zd %Zd\n", message_decoded.x, message_decoded.y);
		gmp_printf("[Encrypt] Message - Jacobian [X Y Z]: %Zd %Zd %Zd\n", message.X, message.Y, message.Z);
		while (mpz_cmp(k_message, zero_value) == 0) {
			get_random(k_message, 32);
			positive_modulo(k_message, k_message, modulo);
		}
		// Encrypt example
		chosen_point = jacobian_left_to_right_binary(j_p, a, k_message, modulo); // chosen point (r)
		gmp_printf("[Encrypt] Chosen point - Jacobian [X Y Z]: %Zd %Zd %Zd\n", chosen_point.X, chosen_point.Y, chosen_point.Z);
		encoded_point = jacobian_left_to_right_binary(public_key_2, a, k_message, modulo); // r * Pu2
		encoded_point = jacobian_curve_addition(message, encoded_point, a, modulo);
		// TODO : chosen_point & encoded_point should be padded to P-bit
		gmp_printf("[Decrypt] Encoded point - Jacobian [X Y Z]: %Zd %Zd %Zd\n", encoded_point.X, encoded_point.Y, encoded_point.Z);
	
		// Decrypt example (encoded_point - private_key * chosen_point)
		decoded_point = jacobian_left_to_right_binary(chosen_point, a, private_key_2, modulo);
		decoded_point = jacobian_curve_substraction(encoded_point, decoded_point, a, modulo);
		gmp_printf("[Decrypt] Original message - Jacobian [X Y Z]: %Zd %Zd %Zd\n", decoded_point.X, decoded_point.Y, decoded_point.Z);
		message_decoded = jacobian_to_affine(decoded_point, modulo);
		gmp_printf("[Decrypt] Original message - Affine [X Y] %Zd %Zd\n\n", message_decoded.x, message_decoded.y);
	}
	/** -------------------------------------------------------------------------*/
	if (TEST_SIMPLIFIED_ECIES) {
		// Simplified ECIES (Ref: Page 256 Cryptography Theory & Practice 2nd Ed. - Douglas)
		char* message_string = "hello"; // 0..9, a..z (base 36)
		mpz_t encrypted_message;
		mpz_init(encrypted_message);

		int partition = strlen(message_string) / 24;
		int partition_modulo = strlen(message_string) % 24;
		if (partition_modulo != 0) partition++;

		for (i = 0; i < partition; i++) {
			// 24 characters from message_string + 1 null-terminate
			char* chunked_message_string = (char*) malloc(25 * sizeof(char));
			int size = 24;
			if ((i == partition - 1) && (partition_modulo != 0)) size = partition_modulo;
			strncpy(chunked_message_string, message_string + i*24, size);
			chunked_message_string[size] = '\0'; // null-terminate

			Point c_point = encrypt_ECIES(encrypted_message, chunked_message_string, public_key_2_decoded, p, a, modulo);
			gmp_printf("[SIMPLIFIED ECIES] Encrypted message: %Zd\n", encrypted_message);
			decrypt_ECIES(encrypted_message, c_point, private_key_2, p, a, modulo);
		}
	}
	/**-------------------------------------------------------------------------*/
	// TODO : Public key validation!
	// Shared key (ECDH) - key secure exchange
	shared_key = jacobian_left_to_right_binary(public_key_2, a, private_key_1, modulo);
	gmp_printf("Shared key - Jacobian [X Y Z]: %Zd %Zd %Zd\n", shared_key.X, shared_key.Y, shared_key.Z);
	Point shared_key_decoded = jacobian_to_affine(shared_key, modulo);
	gmp_printf("Shared key - Affine [X Y]: %Zd %Zd\n", shared_key_decoded.x, shared_key_decoded.y);

	// TODO : ECDSA - digital signature algorithm

	/** Cleaning up */
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(k);
	mpz_clear(r);
	mpz_clear(modulo);
	mpz_clear(private_key_1);
	mpz_clear(private_key_2);

	return EXIT_SUCCESS;
}
/* Created by freedomofkeima - 2014 */

/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */

#include "boolean.h"
#include "random.h"

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

void get_random(mpz_t results, int num_bytes) { // multiple of 8
	unsigned long long *data = (unsigned long long*) malloc((num_bytes / 8) * sizeof(long long));
	int ret_val;
	FILE *fp;
	fp = fopen("/dev/urandom", "r");
	if (fp == NULL) {
		fprintf(stderr, "cannot open random number device!");
		return;
	}
	ret_val = fread(data, 8, num_bytes / 8, fp);
	fclose(fp);

	mpz_init(results);
	mpz_import(results, (num_bytes / 8), 1, sizeof(data[0]), 0, 0, data);

	mpz_t zero_value;
	mpz_init(zero_value);
}

/* Created by freedomofkeima - 2014 */

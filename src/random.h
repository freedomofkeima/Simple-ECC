/**
  * JAIST - Visiting Student 2014
  * Iskandar Setiadi s1416051@jaist.ac.jp
  *
  */


#ifndef RANDOM_H
#define RANDOM_H

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

void get_random(mpz_t results, int num_bytes); // multiple of 8

#endif
/* Created by freedomofkeima - 2014 */

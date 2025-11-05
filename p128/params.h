#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define LIMBS 25
#define MAX_GS 15
#define MAX_K 1109

typedef struct uint_custom { uint64_t c[LIMBS]; } uint_custom;
typedef struct uint_custom2 { uint_custom a, b; } uint_custom2;
typedef struct fp { uint64_t c[LIMBS]; } fp;
typedef struct fp2 {fp a, b;} fp2;
typedef struct proj {fp2 x, z;} proj;

#define NUM_PRIMES 93
#define KEY_STORAGE_SIZE 24 /* 24byte=192bit > {0,1,2}^93 */

static const unsigned A_primes[NUM_PRIMES] = {
      4,   5,  11,  17,  23,  31,  41,  47,  59,  67,  73,  83,  97, 103, 109, 127, 
    137, 149, 157, 167, 179, 191, 197, 211, 227, 233, 241, 257, 269, 277, 283, 307, 
    313, 331, 347, 353, 367, 379, 389, 401, 419, 431, 439, 449, 461, 467, 487, 499, 
    509, 523, 547, 563, 571, 587, 599, 607, 617, 631, 643, 653, 661, 677, 691, 709, 
    727, 739, 751, 761, 773, 797, 811, 823, 829, 853, 859, 877, 883, 907, 919, 937, 
    947, 967, 977, 991, 1009, 1019, 1031, 1039, 1051, 1063, 1087, 1093, 1103,
};

static const unsigned B_primes[NUM_PRIMES] = {
      3,   7,  13,  19,  29,  37,  43,  53,  61,  71,  79,  89, 101, 107, 113, 131, 
    139, 151, 163, 173, 181, 193, 199, 223, 229, 239, 251, 263, 271, 281, 293, 311, 
    317, 337, 349, 359, 373, 383, 397, 409, 421, 433, 443, 457, 463, 479, 491, 503, 
    521, 541, 557, 569, 577, 593, 601, 613, 619, 641, 647, 659, 673, 683, 701, 719, 
    733, 743, 757, 769, 787, 809, 821, 827, 839, 857, 863, 881, 887, 911, 929, 941, 
    953, 971, 983, 997, 1013, 1021, 1033, 1049, 1061, 1069, 1091, 1097, 1109,
};

#include "constants.h"

#endif

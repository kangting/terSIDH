#ifndef terSIDH_H
#define terSIDH_H

#include <stdbool.h>

#include "params.h"

typedef struct private_key {
    int8_t e[NUM_PRIMES]; // 24byte=192bit
} private_key;

typedef struct public_key {
    uint_custom2 A; 
    proj xR, xS; /* A : Montgomery coefficient: represents y^2 = x^3 + Ax^2 + x */
} public_key;

extern const public_key base;

void random_generators_in_subgroup(proj *P, proj *Q, const proj *A, bool want_A_side);

void setup(proj **points, bool Alice);
void tersidh_private(private_key *priv);
bool validate_basic(public_key const *in);

void isogeny(uint8_t number, proj *A, proj **points, int8_t es[NUM_PRIMES], const unsigned primes[NUM_PRIMES]);
void isogeny_reduced(uint8_t number, proj *A, proj **points, int8_t es[NUM_PRIMES], const unsigned primes[NUM_PRIMES]);
void isogeny_constant(uint8_t number, proj *A, proj **points, int8_t es[NUM_PRIMES], const unsigned primes[NUM_PRIMES]);
bool keygen(public_key *out, proj **points, public_key const *in, private_key *priv, bool Alice);
bool shared(fp2 *out, public_key *in, private_key const *priv, bool Alice);

#endif

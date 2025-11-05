#ifndef UINT_custom_H
#define UINT_custom_H

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#include "params.h"

extern const uint_custom uint_custom_0;
extern const uint_custom uint_custom_1;
extern const uint_custom uint_custom_2;

bool uint_custom_eq(uint_custom const *x, uint_custom const *y);

void uint_custom_set(uint_custom *x, uint64_t y);
void uint_custom_set_jh(uint_custom *x, uint_custom const *y);

size_t uint_custom_len(uint_custom const *x);
bool uint_custom_bit(uint_custom const *x, uint64_t k);

bool uint_custom_add3(uint_custom *x, uint_custom const *y, uint_custom const *z); /* returns carry */
bool uint_custom_sub3(uint_custom *x, uint_custom const *y, uint_custom const *z); /* returns borrow */

void uint_custom_mul3_64(uint_custom *x, uint_custom const *y, uint64_t z);

/* q = floor(x / y), r = x mod y, where y fits in 64 bits */
void uint_custom_divmod_64(uint_custom *q, uint64_t *r, uint_custom const *x, uint64_t y);

bool uint_custom_eq(uint_custom const *x, uint_custom const *y);

void uint_custom_random(uint_custom *x, uint_custom const *m);   /* uniform in the interval [0;m) */

void uint_custom_copy(uint_custom *x, uint_custom const *y);

#endif

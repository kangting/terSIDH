#ifndef FP_H
#define FP_H

#include <stdbool.h>

#include "params.h"

extern const fp fp_0;
extern const fp fp_1;

/* fp */
bool fp_eq(fp const *x, fp const *y);
void fp_set(fp *x, uint64_t y);

void fp_enc(fp *x, uint_custom const *y); /* encode to Montgomery representation */
void fp_dec(uint_custom *x, fp const *y); /* decode from Montgomery representation */

void fp_add3(fp *x, fp const *y, fp const *z);
void fp_add2(fp *x, fp const *y);

void fp_sub3(fp *x, fp const *y, fp const *z);
void fp_sub2(fp *x, fp const *y);

void fp_enc(fp *x, uint_custom const *y);
void fp_dec(uint_custom *x, fp const *y);

void fp_mul3(fp *x, fp const *y, fp const *z);
void fp_mul2(fp *x, fp const *y);

void fp_sqrt(fp *x);
bool fp_issquare(fp *x);

void fp_sqr2(fp *x, fp const *y);
void fp_sqr1(fp *x);

void fp_pow(fp *x, uint_custom const *e);

void fp_inv(fp *x);
void fp_sqrt(fp *x);

bool fp_issquare(fp *x);

void fp_random(fp *x);

#endif
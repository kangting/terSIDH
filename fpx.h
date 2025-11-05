#ifndef FPX_H
#define FPX_H

#include <stdbool.h>

#include "params.h"
#include "fp.h"

extern const fp fp_0;
extern const fp fp_1;
extern const fp2 fp2_0;
extern const fp2 fp2_1;

/* fp2 */
bool fp2_eq(fp2 const *x, fp2 const *y);
void fp2_enc(fp2 *x, uint_custom2 const *y);
void fp2_dec(uint_custom2 *x, fp2 const *y);
void fp2_set(fp2 *x, uint64_t y, uint64_t z);
void fp2_copy(fp2* x, const fp2* y);
void fp2_cswap(fp2 *x, fp2 *y, uint64_t c);

void fp2_add3(fp2 *x, fp2 const *y, fp2 const *z);
void fp2_add2(fp2 *x, fp2 const *y);
void fp2_sub3(fp2 *x, fp2 const *y, fp2 const *z);
void fp2_sub2(fp2 *x, fp2 const *y);
void fp2_mul3(fp2 *x, fp2 const *y, fp2 const *z);
void fp2_mul2(fp2 *x, fp2 const *y);
void fp2_sqr2(fp2 *x, fp2 const *y);
void fp2_sqr1(fp2 *x);
void fp2_pow(fp2 *x, uint_custom const *e);
void fp2_inv(fp2 *x);
void fp2_neg1(fp2 *a);
void fp2_neg2(fp2 *c, const fp2 *a);

bool fp2_issquare(fp2 *x); /* destroys input! */
int fp2_sqrt(fp2 *y, fp2 const *x);

void fp2_random(fp2 *x);

/* for test*/
void print_fp(fp x); 
void print_fp2(fp2 x);
void copy_point(proj* P, proj const* Q);

#endif
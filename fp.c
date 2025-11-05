#include <stdio.h>

#include <stddef.h>
#include <string.h>
#include <assert.h>

#include "params.h"
#include "uint_custom.h"
#include "fp.h"

uint64_t *fp_mul_counter = NULL;
uint64_t *fp_sqr_counter = NULL;
uint64_t *fp_inv_counter = NULL;
uint64_t *fp_sqrt_counter = NULL;


/*fp*/
__attribute__((visibility("default")))
bool fp_eq(fp const *x, fp const *y)
{
    uint64_t r = 0;
    for (size_t k = 0; k < LIMBS; ++k)
        r |= x->c[k] ^ y->c[k];
    return !r;
}

void fp_set(fp *x, uint64_t y)
{
    uint_custom_set((uint_custom *) x, y);
    fp_enc(x, (uint_custom *) x);
}

static void reduce_once(uint_custom *x)
{
    uint_custom t;
    if (!uint_custom_sub3(&t, x, &p))
        *x = t;
}

void fp_add3(fp *x, fp const *y, fp const *z)
{
    bool c = uint_custom_add3((uint_custom *) x, (uint_custom *) y, (uint_custom *) z);
    (void) c; assert(!c);
    reduce_once((uint_custom *) x);
}

void fp_add2(fp *x, fp const *y)
{
    fp_add3(x, x, y);
}

void fp_sub3(fp *x, fp const *y, fp const *z)
{
    if (uint_custom_sub3((uint_custom *) x, (uint_custom *) y, (uint_custom *) z))
        uint_custom_add3((uint_custom *) x, (uint_custom *) x, &p);
}

void fp_sub2(fp *x, fp const *y)
{
    fp_sub3(x, x, y);
}


/* Montgomery arithmetic */


void fp_enc(fp *x, uint_custom const *y)
{
    fp_mul3(x, (fp *) y, &r_squared_mod_p);
}

void fp_dec(uint_custom *x, fp const *y)
{
    fp_mul3((fp *) x, y, (fp *) &uint_custom_1);
}

void fp_mul3(fp *x, fp const *y, fp const *z)
{
    uint64_t t[2 * LIMBS] = {0};
    for (size_t k = 0; k < LIMBS; ++k) {

        uint64_t *r = t + k;

        uint64_t temp = y->c[k] * z->c[0] + r[0];
        uint64_t m = inv_min_p_mod_r * temp;

        bool c = 0, o = 0;
        for (size_t i = 0; i < LIMBS; ++i) {
            __uint128_t u = (__uint128_t) m * p.c[i];
            o = __builtin_add_overflow(r[i], o, &r[i]);
            o |= __builtin_add_overflow(r[i], (uint64_t) u, &r[i]);
            c = __builtin_add_overflow(r[i + 1], c, &r[i + 1]);
            c |= __builtin_add_overflow(r[i + 1], (uint64_t)(u >> 64), &r[i + 1]);
        }
        o = __builtin_add_overflow(r[LIMBS], o, &r[LIMBS]);
        assert(!c && !o);

        c = o = 0;
        for (size_t i = 0; i < LIMBS; ++i) {
            __uint128_t u = (__uint128_t) y->c[k] * z->c[i];
            o = __builtin_add_overflow(r[i], o, &r[i]);
            o |= __builtin_add_overflow(r[i], (uint64_t) u, &r[i]);
            c = __builtin_add_overflow(r[i + 1], c, &r[i + 1]);
            c |= __builtin_add_overflow(r[i + 1], (uint64_t)(u >> 64), &r[i + 1]);
        }
        o = __builtin_add_overflow(r[LIMBS], o, &r[LIMBS]);
        assert(!c && !o);
    }

    for (size_t i = 0; i < LIMBS; ++i) {
        x->c[i] = t[LIMBS + i];
    }

    reduce_once((uint_custom *) x);
}

void fp_mul2(fp *x, fp const *y)
{
    fp_mul3(x, x, y);
}

void fp_sqrt(fp *x)
{
    fp_pow(x, &p_plus_1_f);
}

bool fp_issquare(fp *x)
{
    fp tmp;
    fp_pow(&tmp, x, &p_minus_1_halves);
    return fp_eq(&tmp, &fp_1);
}

void fp_sqr2(fp *x, fp const *y)
{
    if (fp_sqr_counter) ++*fp_sqr_counter;
    uint64_t *mulcnt = fp_mul_counter;
    fp_mul_counter = NULL;
    fp_mul3(x, y, y);
    fp_mul_counter = mulcnt;
}

void fp_sqr1(fp *x)
{
    fp_sqr2(x, x);
}

/* (obviously) not constant time in the exponent */

void fp_pow(fp *x, uint_custom const *e)
{
    fp y = *x;
    *x = fp_1;
    for (size_t k = 0; k < LIMBS; ++k) {
        uint64_t t = e->c[k];
        for (size_t i = 0; i < 64; ++i, t >>= 1) {
            if (t & 1)
                fp_mul2(x, &y);
            fp_sqr1(&y);
        }
    }
}

void fp_inv(fp *x)
{
    if (fp_inv_counter) ++*fp_inv_counter;
    uint64_t *mulcnt = fp_mul_counter;
    fp_mul_counter = NULL;
    fp_pow(x, &p_minus_2);
    fp_mul_counter = mulcnt;
}

void fp_random(fp *x)
{
    uint_custom_random((uint_custom *) x, &p);
}
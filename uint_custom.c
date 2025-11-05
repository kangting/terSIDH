#include <stdio.h>

#include <string.h>
#include <stddef.h>
#include <assert.h>

#include "params.h"
#include "uint_custom.h"
#include "rng.h"

/* assumes little-endian throughout. */

bool uint_custom_eq(uint_custom const *x, uint_custom const *y)
{
    uint64_t r = 0;
    for (size_t k = 0; k < LIMBS; ++k)
        r |= x->c[k] ^ y->c[k];
    return !r;
}

void uint_custom_set(uint_custom *x, uint64_t y)
{
    x->c[0] = y;
    for (size_t i = 1; i < LIMBS; ++i)
        x->c[i] = 0;
}

void uint_custom_set_jh(uint_custom *x, uint_custom const *y)
{
    for (size_t i = 0; i < LIMBS; ++i)
        x->c[i] = y->c[i];
}

size_t uint_custom_len(uint_custom const *x)
{
    for (size_t i = LIMBS - 1, j; i < LIMBS; --i) {
        uint64_t v = x->c[i];
        if (!v) continue;
        for (j = 0; v; ++j, v >>= 1);
        return 64*i+j;
    }
    return 0;
}

bool uint_custom_bit(uint_custom const *x,  uint64_t k)
{
    return 1 & (x->c[k / 64] >> k % 64);
}


bool uint_custom_add3(uint_custom *x, uint_custom const *y, uint_custom const *z)
{
    bool c = 0;
    for (size_t i = 0; i < LIMBS; ++i) {
        uint64_t t;
        c = __builtin_add_overflow(y->c[i], c, &t);
        c |= __builtin_add_overflow(t, z->c[i], &x->c[i]);
    }
    return c;
}

bool uint_custom_sub3(uint_custom *x, uint_custom const *y, uint_custom const *z)
{
    bool b = 0;
    for (size_t i = 0; i < LIMBS; ++i) {
        uint64_t t;
        b = __builtin_sub_overflow(y->c[i], b, &t);
        b |= __builtin_sub_overflow(t, z->c[i], &x->c[i]);
    }
    return b;
}

void uint_custom_mul3_64(uint_custom *x, uint_custom const *y, uint64_t z)
{
    uint64_t c = 0;
    for (size_t i = 0; i < LIMBS; ++i) {
        __uint128_t t = y->c[i] * (__uint128_t) z + c;
        c = t >> 64;
        x->c[i] = t;
    }
}


void uint_custom_divmod_64(uint_custom *q, uint64_t *r, uint_custom const *x, uint64_t y)
{
    // Long division by 64-bit divisor y (non-zero). Assumes little-endian limbs.
    __uint128_t rem = 0;
    for (int i = (int)LIMBS - 1; i >= 0; --i) {
        __uint128_t cur = (rem << 64) | ( __uint128_t ) x->c[i];
        uint64_t qdigit = (uint64_t)(cur / y);
        rem = cur - ( (__uint128_t)qdigit * y );
        q->c[i] = qdigit;
    }
    if (r) *r = (uint64_t)rem;
}


void uint_custom_random(uint_custom *x, uint_custom const *m)
{
    if (!m) {
        randombytes(x, sizeof(*x));
        return;
    }

    assert(memcmp(m, &uint_custom_0, sizeof(*m)));

    const size_t bits = uint_custom_len(m);
    const uint64_t mask = ((uint64_t) 1 << bits % 64) - 1;

    while (true) {

        *x = uint_custom_0;
        randombytes(x, (bits+7)/8);

        if (bits < 64*LIMBS)
            x->c[bits/64] &= mask;

        for (size_t i = LIMBS-1; i < LIMBS; --i) {
            if (x->c[i] < m->c[i])
                return;
            else if (x->c[i] > m->c[i])
                break;
        }
    }
}

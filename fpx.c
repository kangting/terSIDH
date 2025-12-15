#include <stdio.h>

#include <stddef.h>
#include <string.h>
#include <assert.h>

#include "params.h"
#include "uint_custom.h"
#include "fp.h"
#include "fpx.h"


/* fp2 */
__attribute__((visibility("default")))
bool fp2_eq(fp2 const *x, fp2 const *y) {
    return fp_eq(&x->a, &y->a) && fp_eq(&x->b, &y->b);
}

__attribute__((visibility("default")))
void fp2_enc(fp2 *x, uint_custom2 const *y) {
    fp_enc(&x->a, &(y->a));
    fp_enc(&x->b, &(y->b));
}

__attribute__((visibility("default")))
void fp2_dec(uint_custom2 *x, fp2 const *y) {
    fp_dec(&(x->a), &y->a);
    fp_dec(&(x->b), &y->b);
}

void fp2_set(fp2 *x, uint64_t y, uint64_t z)
{
    fp_set(&x->a, y);
    fp_set(&x->b, z);
}

void fp2_copy(fp2* x, const fp2* y)
{
    memcpy(&x->a, &y->a, sizeof(fp));
    memcpy(&x->b, &y->b, sizeof(fp));
}

void fp2_cswap(fp2 *x, fp2 *y, uint64_t c)
{
    // c=1 → mask=0xFFFFFFFFFFFFFFFF, c=0 → mask=0x0000000000000000
    uint64_t mask = (uint64_t)(-(int64_t)c);

    for (size_t i = 0; i < LIMBS; i++) {
        uint64_t t = (x->a.c[i] ^ y->a.c[i]) & mask;
        x->a.c[i] ^= t;
        y->a.c[i] ^= t;
    }

    for (size_t i = 0; i < LIMBS; i++) {
        uint64_t t = (x->b.c[i] ^ y->b.c[i]) & mask;
        x->b.c[i] ^= t;
        y->b.c[i] ^= t;
    }
}

void fp2_add3(fp2 *x, fp2 const *y, fp2 const *z) 
{
    fp_add3(&x->a, &y->a, &z->a);
    fp_add3(&x->b, &y->b, &z->b);
}

void fp2_add2(fp2 *x, fp2 const *y) 
{
    fp_add2(&x->a, &y->a);  // x.a += y.a
    fp_add2(&x->b, &y->b);  // x.b += y.b
}

void fp2_sub3(fp2 *x, fp2 const *y, fp2 const *z) 
{
    fp_sub3(&x->a, &y->a, &z->a);
    fp_sub3(&x->b, &y->b, &z->b);
}

void fp2_sub2(fp2 *x, fp2 const *y) 
{
    fp_sub2(&x->a, &y->a);
    fp_sub2(&x->b, &y->b);
}

void fp2_mul3(fp2 *x, fp2 const *y, fp2 const *z) 
{
    fp t1, t2, t3, t4;

    fp_mul3(&t1, &y->a, &z->a);  // t1 = y.a * z.a
    fp_mul3(&t2, &y->b, &z->b);  // t2 = y.b * z.b

    fp_add3(&t3, &y->a, &y->b);  // t3 = y.a + y.b
    fp_add3(&t4, &z->a, &z->b);  // t4 = z.a + z.b

    fp_mul3(&x->b, &t3, &t4);    // x.b = (y.a + y.b) * (z.a + z.b)
    fp_sub2(&x->b, &t1);         // x.b -= t1
    fp_sub2(&x->b, &t2);         // x.b -= t2

    fp_sub3(&x->a, &t1, &t2);    // x.a = t1 - t2
}

void fp2_mul2(fp2 *x, fp2 const *y) 
{
    fp2 temp;
    fp2_mul3(&temp, x, y);
    *x = temp;
}

void fp2_sqr2(fp2 *x, fp2 const *y) 
{
    fp t1, t2, t3;

    fp_add3(&t1, &y->a, &y->b);  // t1 = y.a + y.b
    fp_sub3(&t2, &y->a, &y->b);  // t2 = y.a - y.b
    fp_mul3(&t3, &y->a, &y->b);  // t3 = y.a * y.b

    fp_mul3(&x->a, &t1, &t2);    // x.a = t1 * t2
    fp_add3(&x->b, &t3, &t3);    // x.b = 2 * t3
}

void fp2_sqr1(fp2 *x)
{
    fp2_sqr2(x, x);
}

void fp2_pow(fp2 *x, uint_custom const *e)
{
    fp2 y = *x;
    *x = fp2_1;
    for (size_t k = 0; k < LIMBS; ++k) {
        uint64_t t = e->c[k];
        for (size_t i = 0; i < 64; ++i, t >>= 1) {
            if (t & 1)
                fp2_mul2(x, &y);
            fp2_sqr1(&y);
        }
    }
}

void fp2_inv(fp2 *x) 
{
    fp t1, t2, denominator;

    fp_mul3(&t1, &x->a, &x->a);  // t1 = a^2
    fp_mul3(&t2, &x->b, &x->b);  // t2 = b^2
    fp_add3(&denominator, &t1, &t2);  // denominator = a^2 + b^2

    fp_inv(&denominator);  // denominator = 1/(a^2 + b^2)

    fp a_temp = x->a;
    fp b_temp = x->b;

    fp_mul3(&x->a, &a_temp, &denominator);  // x->a = a / (a^2 + b^2)
    fp_mul3(&x->b, &b_temp, &denominator);  // x->b = b / (a^2 + b^2)
    fp_sub3(&x->b, &fp_0, &x->b); // neg
}

void fp2_neg1(fp2 *a)
{   
    fp_sub3(&a->a, &fp_0, &a->a);
    fp_sub3(&a->b, &fp_0, &a->b);


}

void fp2_neg2(fp2 *c, const fp2 *a)
{   
    fp_sub3(&c->a, &fp_0, &a->a);
    fp_sub3(&c->b, &fp_0, &a->b);
}

bool fp2_issquare(fp2 *x)
{
    fp2_pow(x, &p_minus_1_halves);
    return fp2_eq(x, &fp2_1);
}

int fp2_sqrt(fp2 *y, const fp2 *x){

    fp a = x->a;
    fp b = x->b;

    fp r2, sqrtval, check;
    fp_set(&r2, 2);
    fp_inv(&r2);

    fp c, d;
    fp_sqr2(&c, &a);
    fp_sqr2(&d, &b);
    fp_add2(&c, &d);

    sqrtval = c;

    fp_sqrt(&sqrtval); // sqrtval = sqrt(a^2 + b^2)

    fp_add3(&c, &a, &sqrtval);
    fp_mul2(&c, &r2); // c = (a + sqrt(a^2+b^2)) / 2

    check = c;
    if(!fp_issquare(&check)){
        fp_sub3(&c, &a, &sqrtval);
        fp_mul2(&c, &r2); // c = (a - sqrt(a^2+b^2)) / 2
    }
    
    fp_sqrt(&c); // c = sqrt({a + sqrt(a^2 + b^2)} / 2)

    y->a = c;

    fp_mul3(&d, &c, &fp_2);
    fp_inv(&d);
    fp_mul2(&d, &b); // d = b / 2c

    y->b = d;

    return 1;
}

void fp2_random(fp2 *x)
{
    fp_random(&x->a);
    fp_random(&x->b);
}

// for test, print fp and fp2
__attribute__((visibility("default")))
void print_fp(fp x) {
    for (int i = 0; i < LIMBS; i++) {
        printf("%016lx ", x.c[i]);
    }
    printf("\n");
}

__attribute__((visibility("default")))
void print_fp2(fp2 x) {
    printf("Real part: ");
    print_fp(x.a);
    printf("Imaginary part: ");
    print_fp(x.b);
}

__attribute__((visibility("default")))
void copy_point(proj* P, proj const* Q)
{
    fp2_copy(&(P->x), &(Q->x));
    fp2_copy(&(P->z), &(Q->z));
}

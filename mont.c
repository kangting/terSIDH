
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

#include "params.h"
#include "uint_custom.h"
#include "fpx.h"
#include "mont.h"
#include "rng.h"
#include "poly.h"
#include "steps.h"


uint64_t *xmul_counters = NULL;  /* array indexed by bit length */
uint64_t *isog_counters = NULL;  /* array indexed by degree */

// simultaneous square-and-multiply, computes x^exp and y^exp 
void exp_by_squaring(fp2* x, fp2* y, uint64_t exp)
{
        fp2 result1, result2;
        fp2_set(&result1, 1, 0);
        fp2_set(&result2, 1, 0);

    while (exp)
    {
        if (exp & 1){
          fp2_mul2(&result1, x);
          fp2_mul2(&result2, y);
        }
        
        fp2_sqr1(x);
        fp2_sqr1(y);
        exp >>= 1;
    }

    fp2_cswap(&result1, x, 1);
    fp2_cswap(&result2, y, 1);

}

/* biquad_1 and biquad_minus1 */
void biquad_pm1_opt_fp2(fp2 *coeff, fp2 *outplus,fp2 *outminus,const proj *P,const fp2 *A, const fp2 *C)
{
  fp2 Pplus; fp2_add3(&coeff[3],&P->x,&P->z);
  fp2 Pminus; fp2_sub3(&coeff[4],&P->x,&P->z);
  fp2_sqr2(&Pplus, &coeff[3]);
  fp2_sqr2(&Pminus, &coeff[4]);

  fp2_mul3(&outplus[0],&Pminus,C);
  outplus[2] = outplus[0];
  fp2_mul3(&outminus[0],&Pplus,C);
  outminus[2] = outminus[0];

  fp2_add3(&coeff[0],&outplus[0],&outminus[0] );
  fp2_sub3(&coeff[1], &outminus[0],&outplus[0] );

  fp2 u, fp2zero;
  fp2zero.a = fp_0;
  fp2zero.b = fp_0;
  fp2_sub3(&coeff[2], &Pplus, &Pminus);
  fp2_mul2(&coeff[2],A);
  fp2_sub3(&u, &fp2zero,  &coeff[2]);

  fp2 t;
  fp2_add3(&t,&outminus[0],&outminus[0]);
  fp2_sub3(&outplus[1],&u,&t);

  fp2_add3(&t,&outplus[0],&outplus[0]);
  fp2_sub3(&outminus[1],&t,&u);
}

/* biquad and biquad_inv */
void biquad_both_fp2(fp2 *out,fp2 *outinv, fp2 *coeff, fp2 *coeffQ, const fp2 *C)
{
  fp2 t0, t1, t2;

  fp2_mul3(&t0, &coeff[3], &coeffQ[1]);
  fp2_mul3(&t1, &coeff[4], &coeffQ[0]);

  fp2_add3(&out[0], &t0, &t1);
  fp2_sqr1(&out[0]);
  fp2_mul2(&out[0],C);


  fp2_sub3(&out[2], &t1, &t0);
  fp2_sqr1(&out[2]);
  fp2_mul2(&out[2],C);

  fp2_mul3(&t0, &coeff[1], &coeffQ[3]);
  fp2_mul3(&t1, &coeff[0], &coeffQ[2]);
  fp2_mul3(&t2, &coeff[2], &coeffQ[2]);

  fp2 fp2zero;
  fp2zero.a = fp_0;
  fp2zero.b = fp_0;


  fp2_add3(&t0,&t0,&t1);
  fp2_add3(&t0,&t0,&t2);
  fp2_sub3(&out[1], &fp2zero, &t0);

 
  outinv[1] = out[1];
  outinv[2] = out[0];
  outinv[0] = out[2];
}

bool is_infinity(proj const *P)
{
    return fp2_eq(&P->z, &fp2_0);
}

bool is_affine(proj const *P)
{
    return fp2_eq(&P->z, &fp2_1);
}

void affinize(proj *P, proj *Q)
{
    if (!Q) {
        if (is_infinity(P)) {
            P->x = fp2_1;
            P->z = fp2_0;
        }
        else if (!is_affine(P)) {
            fp2_inv(&P->z);
            fp2_mul2(&P->x, &P->z);
            P->z = fp2_1;
        }
    }
    else if (is_affine(P))
        affinize(Q, NULL);
    else if (is_affine(Q))
        affinize(P, NULL);
    else if (is_infinity(P) || is_infinity(Q)) {
        affinize(P, NULL);
        affinize(Q, NULL);
    }
    else {
        /* batch inversions */
        fp2 t;
        fp2_mul3(&t, &P->z, &Q->z); // Pz * Qz          // C^2
        fp2_inv(&t);                // Pz^-1 * Qz^-1    // 1/C^2
        fp2_mul2(&Q->z, &t);        // Pz^-1            // Qz <- 1/C
        fp2_mul2(&P->x, &Q->z);     // Px * Pz^-1       // Px <- A/C
        fp2_mul2(&P->z, &t);        // Qz^-1            // Pz <- 1/C
        fp2_mul2(&Q->x, &P->z);     // Qx * Qz^-1       // Qx <- A/C 
        P->z = Q->z = fp2_1;
    }
}


void xDBL(proj *Q, proj const *P, proj const *A)
{
    fp2 a, b, c;
    fp2_add3(&a, &P->x, &P->z); // a=x+z
    fp2_sqr1(&a); // a=(x+z)^2
    fp2_sub3(&b, &P->x, &P->z); // b=x-z
    fp2_sqr1(&b); // b=(x-z)^2
    fp2_sub3(&c, &a, &b); // c=4xz
    fp2_add2(&b, &b); fp2_add2(&b, &b); // 16xz
    if (!is_affine(A))
        fp2_mul2(&b, &A->z); // b=16x
    fp2_mul3(&Q->x, &a, &b); // 4Ct0^2 * (t1^2)
    fp2_add3(&a, &A->z, &A->z); // a = 2C
    fp2_add2(&a, &A->x); // A+2C
    fp2_mul2(&a, &c); // (t1^2-t0^2)*(A+2C)
    fp2_add2(&a, &b); // (t1^2-t0^2)*(A+2C) + 4Ct0^2
    fp2_mul3(&Q->z, &a, &c); // 
}

void xDBLv2(proj* Q, proj const* P, proj const* A24)
{
    // This version receives the coefficient value A24 = (A+2C:4C) 
    fp2 t0, t1, t2;

    fp2_add3(&t0, &P->x, &P->z);
    fp2_sqr2(&t0, &t0);
    fp2_sub3(&t1, &P->x, &P->z);
    fp2_sqr2(&t1, &t1);
    fp2_sub3(&t2, &t0, &t1);
    fp2_mul3(&t1, &t1, &A24->z);
    fp2_mul3(&Q->x, &t0, &t1);
    fp2_mul3(&t0, &t2, &A24->x);
    fp2_add3(&t0, &t0, &t1);
    fp2_mul3(&Q->z, &t0, &t2);

}

void swap_points(proj* P, proj* Q, const uint64_t option)
{ // Swap points
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    uint64_t temp;

    for (int i = 0; i < LIMBS; i++) {
        temp = option & (P->x.a.c[i] ^ Q->x.a.c[i]);
        P->x.a.c[i] = temp ^ P->x.a.c[i];
        Q->x.a.c[i] = temp ^ Q->x.a.c[i];
        temp = option & (P->x.b.c[i] ^ Q->x.b.c[i]);
        P->x.b.c[i] = temp ^ P->x.b.c[i];
        Q->x.b.c[i] = temp ^ Q->x.b.c[i];
        temp = option & (P->z.a.c[i] ^ Q->z.a.c[i]);
        P->z.a.c[i] = temp ^ P->z.a.c[i];
        Q->z.a.c[i] = temp ^ Q->z.a.c[i];
        temp = option & (P->z.b.c[i] ^ Q->z.b.c[i]);
        P->z.b.c[i] = temp ^ P->z.b.c[i];
        Q->z.b.c[i] = temp ^ Q->z.b.c[i];
    }
}

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp2 a, b, c, d;
    fp2_add3(&a, &P->x, &P->z);
    fp2_sub3(&b, &P->x, &P->z);
    fp2_add3(&c, &Q->x, &Q->z);
    fp2_sub3(&d, &Q->x, &Q->z);
    fp2_mul2(&a, &d);
    fp2_mul2(&b, &c);
    fp2_add3(&S->x, &a, &b);
    fp2_sub3(&S->z, &a, &b);
    fp2_sqr1(&S->x);
    fp2_sqr1(&S->z);
    if (!is_affine(PQ))
        fp2_mul2(&S->x, &PQ->z);
    fp2_mul2(&S->z, &PQ->x);
}

void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    fp2 a, b, c, d;

    fp2_add3(&a, &Q->x, &Q->z);
    fp2_sub3(&b, &Q->x, &Q->z);
    fp2_add3(&c, &P->x, &P->z);
    fp2_sub3(&d, &P->x, &P->z);
    fp2_sqr2(&R->x, &c);
    fp2_sqr2(&S->x, &d);
    fp2_mul2(&c, &b);
    fp2_mul2(&d, &a);
    fp2_sub3(&b, &R->x, &S->x);
    fp2_add3(&a, &A->z, &A->z);
    if (!is_affine(A))
        fp2_mul3(&R->z, &S->x, &a);
    else
        fp2_add3(&R->z, &S->x, &S->x);
    fp2_add3(&S->x, &A->x, &a);
    fp2_add2(&R->z, &R->z);
    fp2_mul2(&R->x, &R->z);
    fp2_mul2(&S->x, &b);
    fp2_sub3(&S->z, &c, &d);
    fp2_add2(&R->z, &S->x);
    fp2_add3(&S->x, &c, &d);
    fp2_mul2(&R->z, &b);
    fp2_sqr1(&S->x);
    fp2_sqr1(&S->z);
    if (!is_affine(PQ))
        fp2_mul2(&S->x, &PQ->z);
    fp2_mul2(&S->z, &PQ->x);
}


/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
void xMUL(proj *Q, proj const *A, proj const *P, uint_custom const *k)
{
    size_t i = uint_custom_len(k);
    if (xmul_counters) ++xmul_counters[i];

    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp2_1;
    Q->z = fp2_0;

    while (i --> 0) {

        const bool bit = uint_custom_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */

        xDBLADD(Q, &R, Q, &R, &Pcopy, A);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */

    }
}

// Fixed-iteration Montgomery ladder; runs exactly max_bits rounds
void xMUL_fixedbits(proj *Q, proj const *A, proj const *P, uint_custom const *k, size_t max_bits)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp2_1;
    Q->z = fp2_0;

    for (size_t i = max_bits; i-- > 0; ) {
        uint64_t bit = (uint64_t)uint_custom_bit(k, i);
        uint64_t mask = (uint64_t)0 - bit;  // 0 or all-ones
        swap_points(Q, &R, mask);
        xDBLADD(Q, &R, Q, &R, &Pcopy, A);
        swap_points(Q, &R, mask);
    }
}

/* Montgomery <-> Edwards conversions */

static void mont2edw(proj *d, proj const *A)
{
    /* d = (A-2)/(A+2) */
    fp2_sub3(&d->x, &A->x, &A->z);
    fp2_sub2(&d->x, &A->z);
    fp2_add3(&d->z, &A->x, &A->z);
    fp2_add2(&d->z, &A->z);
}

static void edw2mont(proj *A, proj const *d)
{
    /* A = 2*(1+d)/(1-d) */
    fp2_add3(&A->x, &d->z, &d->x);
    fp2_add2(&A->x, &A->x);
    fp2_sub3(&A->z, &d->z, &d->x);
}

/* Montgomery isogeny evaluation [Costello-Hisil, Renes] */

static void montisog_eval_init(proj *Q, const proj *P, const proj *K)
{
    fp2 t;

    fp2_mul3(&Q->x, &P->x, &K->x);
    fp2_mul3(&t, &P->z, &K->z);
    fp2_sub2(&Q->x, &t);

    fp2_mul3(&Q->z, &P->x, &K->z);
    fp2_mul3(&t, &P->z, &K->x);
    fp2_sub2(&Q->z, &t);
}

static void montisog_eval_consume(proj *Q, const proj *P, const proj *M)
{
    fp2 t0, t1, t2;

    fp2_sub3(&t0, &P->x, &P->z);
    fp2_add3(&t1, &M->x, &M->z);
    fp2_mul2(&t0, &t1);

    fp2_sub3(&t1, &M->x, &M->z);
    fp2_add3(&t2, &P->x, &P->z);
    fp2_mul2(&t1, &t2);

    fp2_add3(&t2, &t0, &t1);
    fp2_mul2(&Q->x, &t2);

    fp2_sub3(&t2, &t0, &t1);
    fp2_mul2(&Q->z, &t2);
}

static void montisog_eval_finish(proj *P, proj *Q)  /* destroys Q */
{
    fp2_sqr1(&Q->x);
    fp2_sqr1(&Q->z);
    fp2_mul2(&P->x, &Q->x);
    fp2_mul2(&P->z, &Q->z);
}

/* Edwards isogeny codomain [Moody-Shumow] */

static void edwisog_curve_init(proj *prod, proj const *K)
{
    fp2_sub3(&prod->x, &K->x, &K->z);
    fp2_add3(&prod->z, &K->x, &K->z);
}

static void edwisog_curve_consume(proj *prod, proj const *M)
{
    fp2 Y, Z;

    fp2_sub3(&Y, &M->x, &M->z);
    fp2_add3(&Z, &M->x, &M->z);

    fp2_mul2(&prod->x, &Y);
    fp2_mul2(&prod->z, &Z);
}

static void edwisog_curve_finish(proj *d, proj const *prod, uint64_t l)
{

    proj newd = *prod;
    fp2_sqr1(&newd.x); fp2_sqr1(&newd.z);   /* ^2 */
    fp2_sqr1(&newd.x); fp2_sqr1(&newd.z);   /* ^4 */
    fp2_sqr1(&newd.x); fp2_sqr1(&newd.z);   /* ^8 */

    /* square-and-multiply to compute d^l */
    for (uint64_t k = l; k; k >>= 1) {
        if (k & 1) {
            fp2_mul2(&newd.x, &d->x);
            fp2_mul2(&newd.z, &d->z);
        }
        fp2_sqr1(&d->x);
        fp2_sqr1(&d->z);
    }

    *d = newd;
}


void eval_4_isog(proj *P, fp2 *K1, fp2 *K2, fp2 *K3)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined 
  // by the 3 coefficients in coeff (computed in the function get_4_isog()).
  // Inputs: the coefficients defining the isogeny, and the projective point P = (X:Z).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    fp2 t0, t1;
    
    fp2_add3(&t0, &P->x, &P->z);                    // t0 = X+Z
    fp2_sub3(&t1, &P->x, &P->z);                    // t1 = X-Z
    fp2_mul3(&P->x, &t0, K2);                       // X = (X+Z)*K2
    fp2_mul3(&P->z, &t1, K3);                       // Z = (X-Z)*K3
    fp2_mul3(&t0, &t0, &t1);                        // t0 = (X+Z)*(X-Z)
    fp2_mul3(&t0, &t0, K1);                         // t0 = K1*(X+Z)*(X-Z)
    fp2_add3(&t1, &P->x, &P->z);                    // t1 = (X-Z)*K3 + (X+Z)*K2
    fp2_sub3(&P->z, &P->x, &P->z);                  // Z = (X-Z)*K2 - (X+Z)*K1
    fp2_sqr1(&t1);                                  // t1 = [(X-Z)*K3 + (X+Z)*K2^2
    fp2_sqr1(&P->z);                                // Z = [(X-Z)*K3 - (X+Z)*K2]^2
    fp2_add3(&P->x, &t0, &t1);                      // X = K1*(X+Z)*(X-Z) + [(X-Z)*K3 + (X+Z)*K2]^2
    fp2_sub3(&t0, &P->z, &t0);                      // t0 = [(X-Z)*K3 - (X+Z)*K2]^2 -K1*(X+Z)*(X-Z)
    fp2_mul3(&P->x, &P->x, &t1);                    // Xfinal
    fp2_mul3(&P->z, &P->z, &t0);                    // Zfinal
}

// Degree-4 isogeny with kernel generated by P such that [2]P != (0 ,0)
// Outputs the curve coefficient in the form A24=(A+2C:4C)
void xISOG_4(proj *A, proj **points, proj *K, int num)
{   
    /* computing curve */
    fp2 K1, K2, K3;
    fp2_sub3(&K2, &K->x, &K->z);                // K2 = X-Z
    fp2_add3(&K3, &K->x, &K->z);                // K3 = X+Z
    fp2_sqr2(&K1, &K->z);                       // K1 = Z^2
    fp2_add2(&K1, &K1);                         // K1 = K1+K1
    fp2_sqr2(&A->z, &K1);                       // C = K1^2
    fp2_add2(&K1, &K1);                         // K1 = K1+K1
    fp2_sqr2(&A->x, &K->x);                     // A = X^2
    fp2_add2(&A->x, &A->x);                     // A = A+A
    fp2_sqr1(&A->x);                            // A = A^2

    /* evaluation */
    int start = (num % 2 == 1) ? 1 : 0;
    for (int i = start; i < start + num; i++) {
        eval_4_isog(points[i], &K1, &K2, &K3);
    }

    fp2_add2(&A->x, &A->x);
    fp2_sub2(&A->x, &A->z);
    fp2_add2(&A->x, &A->x);
}

void xISOG_4pts(proj *A, proj **points, proj const *K, long long k)
{
    assert(k >= 3);
    assert(k % 2 == 1);

    int sq4tvelu = 0;
    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);

    if (bs) {
        sq4tvelu = 1;
        assert(bs > 0);
        assert(gs > 0);
        assert(!(bs & 1));
    }

    proj Aed;
    fp2_add3(&Aed.z, &A->z, &A->z);   // tmp = 2C
    fp2_add3(&Aed.x, &A->x, &Aed.z);  // Aed.x = A + 2C
    fp2_sub3(&Aed.z, &A->x, &Aed.z);  // Aed.z = A - 2C

    fp2 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13;
    fp2 Psum1, Pdif1, Psum2, Pdif2, Psum3, Pdif3, Psum4, Pdif4;
    fp2 coeffQ1[4], coeffQ2[4], coeffQ3[4], coeffQ4[4];
    fp2_add3(&Psum1, &points[0]->x, &points[0]->z);   //precomputations
    fp2_sub3(&Pdif1, &points[0]->x, &points[0]->z);

    fp2_add3(&Psum2, &points[1]->x, &points[1]->z);   //precomputations
    fp2_sub3(&Pdif2, &points[1]->x, &points[1]->z);

    fp2_add3(&Psum3, &points[2]->x, &points[2]->z);   //precomputations
    fp2_sub3(&Pdif3, &points[2]->x, &points[2]->z);

    fp2_add3(&Psum4, &points[3]->x, &points[3]->z);   //precomputations
    fp2_sub3(&Pdif4, &points[3]->x, &points[3]->z);

    fp2_add3(&coeffQ1[0], &points[0]->x, &points[0]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ1[1], &points[0]->x, &points[0]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ1[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ1[1]); //(x-z)^2
    fp2_add3(&coeffQ1[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ1[2], &tmp0, &tmp1); // 


    fp2_add3(&coeffQ2[0], &points[1]->x, &points[1]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ2[1], &points[1]->x, &points[1]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ2[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ2[1]); //(x-z)^2
    fp2_add3(&coeffQ2[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ2[2], &tmp0, &tmp1); // 


    fp2_add3(&coeffQ3[0], &points[2]->x, &points[2]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ3[1], &points[2]->x, &points[2]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ3[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ3[1]); //(x-z)^2
    fp2_add3(&coeffQ3[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ3[2], &tmp0, &tmp1); // 

    fp2_add3(&coeffQ4[0], &points[3]->x, &points[3]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ4[1], &points[3]->x, &points[3]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ4[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ4[1]); //(x-z)^2
    fp2_add3(&coeffQ4[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ4[2], &tmp0, &tmp1); // 

    int Minit[k];
    proj M[k];

    for (long long s = 0; s < k; s++) {
        Minit[s] = 0;
    }
    M[1] = *K;
    Minit[1] = 1;

    xDBL(&M[2], K, A);
    Minit[2] = 1;

    if (sq4tvelu) {
        for (long long s = 3; s < k; ++s) {
            if (s & 1) {
                long long i = s / 2;  
                assert(s == 2*i + 1);

                if (i < bs) {
                    if (s == 3) {
                        assert(Minit[1]);
                        assert(Minit[2]);
                        xADD(&M[s],
                             &M[2],
                             &M[1],
                             &M[1]);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            } 
            else {
                long long i = s/2 - 1; 
                assert(s == 2*i + 2);

                if (i < (k-1)/2 - 2*bs*gs) {
                    if (s == 4) {
                        assert(Minit[2]);
                        xDBL(&M[s], &M[2], A);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            }
            if (bs > 0) {
                if (s == 2*bs) {
                    assert(Minit[bs-1]);
                    assert(Minit[bs+1]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[bs+1],
                         &M[bs-1],
                         &M[2]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 4*bs) {
                    assert(Minit[2*bs]);
                    xDBL(&M[s], &M[2*bs], A);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 6*bs) {
                    assert(Minit[2*bs]);
                    assert(Minit[4*bs]);
                    xADD(&M[s],
                         &M[4*bs],
                         &M[2*bs],
                         &M[2*bs]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s % (4*bs) == 2*bs) {
                    long long j = s / (4*bs);
                    assert(s == 2*bs*(2*j + 1));
                    if (j < gs) {
                        assert(Minit[s-4*bs]);
                        assert(Minit[s-8*bs]);
                        assert(Minit[4*bs]);
                        xADD(&M[s],
                             &M[s-4*bs],
                             &M[4*bs],
                             &M[s-8*bs]);
                        Minit[s] = 1;
                        continue;
                    }
                }
            }
        }
    }
    else {
        for (long long i = 3; i <= (k-1)/2; i++) {
            Minit[i] = 1;
            xADD(&M[i],
                 &M[i-1],
                 K,
                 &M[i-2]);
        }
    }

    proj Abatch;
    Abatch.x = fp2_1; 
    Abatch.z = fp2_1;

    proj Q1, Q2, Q3, Q4;
    Q1.x = fp2_1;
    Q1.z = fp2_1;

    Q2.x = fp2_1;
    Q2.z = fp2_1;

    Q3.x = fp2_1;
    Q3.z = fp2_1;

    Q4.x = fp2_1;
    Q4.z = fp2_1;

    if (sq4tvelu) {
        long long TIlen = 2*bs + poly_tree1size(bs);
        fp2 TI[TIlen];
        for (long long i = 0; i < bs; ++i) {
            assert(Minit[2*i+1]);
            fp2_neg2(&TI[2*i], &M[2*i+1].x);
            fp2_copy(&TI[2*i+1], &M[2*i+1].z);
        }
        poly_tree1(TI + 2*bs, TI, bs);

        fp2 TP1[3*gs];
        fp2 TPinv1[3*gs];
        fp2 TP2[3*gs];
        fp2 TPinv2[3*gs];
        fp2 TP3[3*gs];
        fp2 TPinv3[3*gs];
        fp2 TP4[3*gs];
        fp2 TPinv4[3*gs];
        fp2 T1[3*gs];
        fp2 Tminus1[3*gs];
        fp2 coeff[5];        

        for (long long j = 0;j < gs;++j) {
            assert(Minit[2*bs*(2*j+1)]);
            biquad_pm1_opt_fp2(coeff, T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],&A->x,&A->z);
            biquad_both_fp2(TP1+3*j,TPinv1+3*j, coeff, coeffQ1, &A->z);
            biquad_both_fp2(TP2+3*j,TPinv2+3*j, coeff, coeffQ2, &A->z);
            biquad_both_fp2(TP3+3*j,TPinv3+3*j, coeff, coeffQ3, &A->z);
            biquad_both_fp2(TP4+3*j,TPinv4+3*j, coeff, coeffQ4, &A->z);

        }

        poly_multiprod2(TP1,gs);
        poly_multiprod2(TPinv1,gs);


        poly_multiprod2(TP2,gs);
        poly_multiprod2(TPinv2,gs);


        poly_multiprod2(TP3,gs);
        poly_multiprod2(TPinv3,gs);


        poly_multiprod2(TP4,gs);
        poly_multiprod2(TPinv4,gs);

        poly_multiprod2_selfreciprocal(T1,gs);
        poly_multiprod2_selfreciprocal(Tminus1,gs);

        

        long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
        fp2 precomp[precompsize];
        poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

        fp2 v[bs];

        poly_multieval_postcompute(v,bs,TP1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.x,&v[i]);



        poly_multieval_postcompute(v,bs,TP2,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q2.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q2.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv2,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q2.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q2.x,&v[i]);



        poly_multieval_postcompute(v,bs,TP3,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q3.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q3.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv3,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q3.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q3.x,&v[i]);



        poly_multieval_postcompute(v,bs,TP4,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q4.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q4.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv4,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q4.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q4.x,&v[i]);



        poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.x,&v[i]);

        poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.z,&v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.z,&v[i]);

        for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
            assert(Minit[2*i+2]);
            fp2_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
            fp2_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
            fp2_mul2(&Abatch.x,&tmp1);
            fp2_mul2(&Abatch.z,&tmp0);


            fp2_mul3(&tmp2, &tmp1, &Psum1);
            fp2_mul3(&tmp3, &tmp0, &Pdif1);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q1.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q1.z, &tmp4);

            fp2_mul3(&tmp2, &tmp1, &Psum2);
            fp2_mul3(&tmp3, &tmp0, &Pdif2);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q2.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q2.z, &tmp4);

            fp2_mul3(&tmp2, &tmp1, &Psum3);
            fp2_mul3(&tmp3, &tmp0, &Pdif3);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q3.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q3.z, &tmp4);

            fp2_mul3(&tmp2, &tmp1, &Psum4);
            fp2_mul3(&tmp3, &tmp0, &Pdif4);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q4.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q4.z, &tmp4);
        }
    } 
    else {
        for (long long i = 1;i <= (k-1)/2;++i) {
            assert(Minit[i]);
            fp2_sub3(&tmp1, &M[i].x, &M[i].z);
            fp2_add3(&tmp0, &M[i].x, &M[i].z);
            if (i > 1) {
                fp2_mul2(&Abatch.x,&tmp1);
                fp2_mul2(&Abatch.z,&tmp0);
            } else {
                Abatch.x = tmp1;
                Abatch.z = tmp0;
            }
            fp2_mul3(&tmp2, &tmp1, &Psum1); 
            fp2_mul3(&tmp3, &tmp0, &Pdif1);
            fp2_add3(&tmp4, &tmp3, &tmp2);

            fp2_mul3(&tmp5, &tmp1, &Psum2); 
            fp2_mul3(&tmp6, &tmp0, &Pdif2);
            fp2_add3(&tmp7, &tmp6, &tmp5);

            fp2_mul3(&tmp8, &tmp1, &Psum3); 
            fp2_mul3(&tmp9, &tmp0, &Pdif3);
            fp2_add3(&tmp10, &tmp9, &tmp8);

            fp2_mul3(&tmp11, &tmp1, &Psum4); 
            fp2_mul3(&tmp12, &tmp0, &Pdif4);
            fp2_add3(&tmp13, &tmp11, &tmp12);

            if (i > 1) {
                fp2_mul2(&Q1.x, &tmp4);
                fp2_mul2(&Q2.x, &tmp7);
                fp2_mul2(&Q3.x, &tmp10);
                fp2_mul2(&Q4.x, &tmp13);
            } else {
                Q1.x = tmp4;
                Q2.x = tmp7;
                Q3.x = tmp10;
                Q4.x = tmp13;
            }
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_sub3(&tmp7, &tmp6, &tmp5);
            fp2_sub3(&tmp10, &tmp9, &tmp8);
            fp2_sub3(&tmp13, &tmp11, &tmp12);
            if (i > 1) {
                fp2_mul2(&Q1.z, &tmp4);
                fp2_mul2(&Q2.z, &tmp7);
                fp2_mul2(&Q3.z, &tmp10);
                fp2_mul2(&Q4.z, &tmp13);
            } 
            else {
                Q1.z = tmp4;
                Q2.z = tmp7;
                Q3.z = tmp10;
                Q4.z = tmp13;
            }
        }
    }
    // point evaluation
    fp2_sqr1(&Q1.x);
    fp2_sqr1(&Q1.z);
    fp2_mul2(&points[0]->x, &Q1.x);
    fp2_mul2(&points[0]->z, &Q1.z);

    fp2_sqr1(&Q2.x);
    fp2_sqr1(&Q2.z);
    fp2_mul2(&points[1]->x, &Q2.x);
    fp2_mul2(&points[1]->z, &Q2.z);

    fp2_sqr1(&Q3.x);
    fp2_sqr1(&Q3.z);
    fp2_mul2(&points[2]->x, &Q3.x);
    fp2_mul2(&points[2]->z, &Q3.z);

    fp2_sqr1(&Q4.x);
    fp2_sqr1(&Q4.z);
    fp2_mul2(&points[3]->x, &Q4.x);
    fp2_mul2(&points[3]->z, &Q4.z);    

    exp_by_squaring(&Aed.x, &Aed.z, k);

    fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); // ^8
    fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); // ^8

    fp2_mul2(&Aed.z, &Abatch.x);
    fp2_mul2(&Aed.x, &Abatch.z);

    fp2_add3(&A->x, &Aed.x, &Aed.z);  // A = Aed.x + Aed.z = 2A'
    fp2_sub3(&A->z, &Aed.x, &Aed.z);  // A = Aed.x - Aed.z = 4C'
    fp2_add2(&A->x, &A->x);          // A = 4A'
    
}

void xISOG_3pts(proj *A, proj **points, proj const *K, long long k)
{
    assert(k >= 3);
    assert(k % 2 == 1);

    int sq4tvelu = 0;
    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);

    if (bs) {
        sq4tvelu = 1;
        assert(bs > 0);
        assert(gs > 0);
        assert(!(bs & 1));
    }

    proj Aed;
    fp2_add3(&Aed.z, &A->z, &A->z);   // tmp = 2C
    fp2_add3(&Aed.x, &A->x, &Aed.z);  // Aed.x = A + 2C
    fp2_sub3(&Aed.z, &A->x, &Aed.z);  // Aed.z = A - 2C

    fp2 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10;
    fp2 Psum1, Pdif1, Psum2, Pdif2, Psum3, Pdif3;
    fp2 coeffQ1[4], coeffQ2[4], coeffQ3[4];
    fp2_add3(&Psum1, &points[1]->x, &points[1]->z);   //precomputations
    fp2_sub3(&Pdif1, &points[1]->x, &points[1]->z);

    fp2_add3(&Psum2, &points[2]->x, &points[2]->z);   //precomputations
    fp2_sub3(&Pdif2, &points[2]->x, &points[2]->z);

    fp2_add3(&Psum3, &points[3]->x, &points[3]->z);   //precomputations
    fp2_sub3(&Pdif3, &points[3]->x, &points[3]->z);


    fp2_add3(&coeffQ1[0], &points[1]->x, &points[1]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ1[1], &points[1]->x, &points[1]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ1[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ1[1]); //(x-z)^2
    fp2_add3(&coeffQ1[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ1[2], &tmp0, &tmp1); // 


    fp2_add3(&coeffQ2[0], &points[2]->x, &points[2]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ2[1], &points[2]->x, &points[2]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ2[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ2[1]); //(x-z)^2
    fp2_add3(&coeffQ2[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ2[2], &tmp0, &tmp1); // 


    fp2_add3(&coeffQ3[0], &points[3]->x, &points[3]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ3[1], &points[3]->x, &points[3]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ3[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ3[1]); //(x-z)^2
    fp2_add3(&coeffQ3[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ3[2], &tmp0, &tmp1); // 
 

    int Minit[k];
    proj M[k];

    for (long long s = 0; s < k; s++) {
        Minit[s] = 0;
    }
    M[1] = *K;
    Minit[1] = 1;

    xDBL(&M[2], K, A);
    Minit[2] = 1;

    if (sq4tvelu) {
        for (long long s = 3; s < k; ++s) {
            if (s & 1) {
                long long i = s / 2;  
                assert(s == 2*i + 1);

                if (i < bs) {
                    if (s == 3) {
                        assert(Minit[1]);
                        assert(Minit[2]);
                        xADD(&M[s],
                             &M[2],
                             &M[1],
                             &M[1]);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            } 
            else {
                long long i = s/2 - 1; 
                assert(s == 2*i + 2);

                if (i < (k-1)/2 - 2*bs*gs) {
                    if (s == 4) {
                        assert(Minit[2]);
                        xDBL(&M[s], &M[2], A);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            }
            if (bs > 0) {
                if (s == 2*bs) {
                    assert(Minit[bs-1]);
                    assert(Minit[bs+1]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[bs+1],
                         &M[bs-1],
                         &M[2]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 4*bs) {
                    assert(Minit[2*bs]);
                    xDBL(&M[s], &M[2*bs], A);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 6*bs) {
                    assert(Minit[2*bs]);
                    assert(Minit[4*bs]);
                    xADD(&M[s],
                         &M[4*bs],
                         &M[2*bs],
                         &M[2*bs]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s % (4*bs) == 2*bs) {
                    long long j = s / (4*bs);
                    assert(s == 2*bs*(2*j + 1));
                    if (j < gs) {
                        assert(Minit[s-4*bs]);
                        assert(Minit[s-8*bs]);
                        assert(Minit[4*bs]);
                        xADD(&M[s],
                             &M[s-4*bs],
                             &M[4*bs],
                             &M[s-8*bs]);
                        Minit[s] = 1;
                        continue;
                    }
                }
            }
        }
    }
    else {
        for (long long i = 3; i <= (k-1)/2; i++) {
            Minit[i] = 1;
            xADD(&M[i],
                 &M[i-1],
                 K,
                 &M[i-2]);
        }
    }

    proj Abatch;
    Abatch.x = fp2_1; 
    Abatch.z = fp2_1;

    proj Q1, Q2, Q3;
    Q1.x = fp2_1;
    Q1.z = fp2_1;

    Q2.x = fp2_1;
    Q2.z = fp2_1;

    Q3.x = fp2_1;
    Q3.z = fp2_1;

    if (sq4tvelu) {
        long long TIlen = 2*bs + poly_tree1size(bs);
        fp2 TI[TIlen];
        for (long long i = 0; i < bs; ++i) {
            assert(Minit[2*i+1]);
            fp2_neg2(&TI[2*i], &M[2*i+1].x);
            fp2_copy(&TI[2*i+1], &M[2*i+1].z);
        }
        poly_tree1(TI + 2*bs, TI, bs);

        fp2 TP1[3*gs];
        fp2 TPinv1[3*gs];
        fp2 TP2[3*gs];
        fp2 TPinv2[3*gs];
        fp2 TP3[3*gs];
        fp2 TPinv3[3*gs];
        fp2 T1[3*gs];
        fp2 Tminus1[3*gs];
        fp2 coeff[5];        

        for (long long j = 0;j < gs;++j) {
            assert(Minit[2*bs*(2*j+1)]);
            biquad_pm1_opt_fp2(coeff, T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],&A->x,&A->z);
            biquad_both_fp2(TP1+3*j,TPinv1+3*j, coeff, coeffQ1, &A->z);
            biquad_both_fp2(TP2+3*j,TPinv2+3*j, coeff, coeffQ2, &A->z);
            biquad_both_fp2(TP3+3*j,TPinv3+3*j, coeff, coeffQ3, &A->z);

        }

        poly_multiprod2(TP1,gs);
        poly_multiprod2(TPinv1,gs);


        poly_multiprod2(TP2,gs);
        poly_multiprod2(TPinv2,gs);


        poly_multiprod2(TP3,gs);
        poly_multiprod2(TPinv3,gs);


        poly_multiprod2_selfreciprocal(T1,gs);
        poly_multiprod2_selfreciprocal(Tminus1,gs);

        

        long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
        fp2 precomp[precompsize];
        poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

        fp2 v[bs];

        poly_multieval_postcompute(v,bs,TP1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.x,&v[i]);



        poly_multieval_postcompute(v,bs,TP2,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q2.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q2.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv2,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q2.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q2.x,&v[i]);



        poly_multieval_postcompute(v,bs,TP3,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q3.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q3.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv3,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q3.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q3.x,&v[i]);




        poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.x,&v[i]);

        poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.z,&v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.z,&v[i]);

        for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
            assert(Minit[2*i+2]);
            fp2_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
            fp2_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
            fp2_mul2(&Abatch.x,&tmp1);
            fp2_mul2(&Abatch.z,&tmp0);


            fp2_mul3(&tmp2, &tmp1, &Psum1);
            fp2_mul3(&tmp3, &tmp0, &Pdif1);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q1.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q1.z, &tmp4);

            fp2_mul3(&tmp2, &tmp1, &Psum2);
            fp2_mul3(&tmp3, &tmp0, &Pdif2);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q2.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q2.z, &tmp4);

            fp2_mul3(&tmp2, &tmp1, &Psum3);
            fp2_mul3(&tmp3, &tmp0, &Pdif3);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q3.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q3.z, &tmp4);

        }
    } 
    else {
        for (long long i = 1;i <= (k-1)/2;++i) {
            assert(Minit[i]);
            fp2_sub3(&tmp1, &M[i].x, &M[i].z);
            fp2_add3(&tmp0, &M[i].x, &M[i].z);
            if (i > 1) {
                fp2_mul2(&Abatch.x,&tmp1);
                fp2_mul2(&Abatch.z,&tmp0);
            } else {
                Abatch.x = tmp1;
                Abatch.z = tmp0;
            }
            fp2_mul3(&tmp2, &tmp1, &Psum1); 
            fp2_mul3(&tmp3, &tmp0, &Pdif1);
            fp2_add3(&tmp4, &tmp3, &tmp2);

            fp2_mul3(&tmp5, &tmp1, &Psum2); 
            fp2_mul3(&tmp6, &tmp0, &Pdif2);
            fp2_add3(&tmp7, &tmp6, &tmp5);

            fp2_mul3(&tmp8, &tmp1, &Psum3); 
            fp2_mul3(&tmp9, &tmp0, &Pdif3);
            fp2_add3(&tmp10, &tmp9, &tmp8);


            if (i > 1) {
                fp2_mul2(&Q1.x, &tmp4);
                fp2_mul2(&Q2.x, &tmp7);
                fp2_mul2(&Q3.x, &tmp10);
            } else {
                Q1.x = tmp4;
                Q2.x = tmp7;
                Q3.x = tmp10;
            }
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_sub3(&tmp7, &tmp6, &tmp5);
            fp2_sub3(&tmp10, &tmp9, &tmp8);
            
            if (i > 1) {
                fp2_mul2(&Q1.z, &tmp4);
                fp2_mul2(&Q2.z, &tmp7);
                fp2_mul2(&Q3.z, &tmp10);
            } 
            else {
                Q1.z = tmp4;
                Q2.z = tmp7;
                Q3.z = tmp10;
            }
        }
    }
    // point evaluation
    fp2_sqr1(&Q1.x);
    fp2_sqr1(&Q1.z);
    fp2_mul2(&points[1]->x, &Q1.x);
    fp2_mul2(&points[1]->z, &Q1.z);

    fp2_sqr1(&Q2.x);
    fp2_sqr1(&Q2.z);
    fp2_mul2(&points[2]->x, &Q2.x);
    fp2_mul2(&points[2]->z, &Q2.z);

    fp2_sqr1(&Q3.x);
    fp2_sqr1(&Q3.z);
    fp2_mul2(&points[3]->x, &Q3.x);
    fp2_mul2(&points[3]->z, &Q3.z);   

    exp_by_squaring(&Aed.x, &Aed.z, k);

    fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); // ^8
    fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); // ^8

    fp2_mul2(&Aed.z, &Abatch.x);
    fp2_mul2(&Aed.x, &Abatch.z);

    fp2_add3(&A->x, &Aed.x, &Aed.z);  // A = Aed.x + Aed.z = 2A'
    fp2_sub3(&A->z, &Aed.x, &Aed.z);  // A = Aed.x - Aed.z = 4C'
    fp2_add2(&A->x, &A->x);          // A = 4A'
    
}

void xISOG_2pts(proj *A, proj **points, proj const *K, long long k)
{
    assert(k >= 3);
    assert(k % 2 == 1);

    int sq4tvelu = 0;
    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);

    if (bs) {
        sq4tvelu = 1;
        assert(bs > 0);
        assert(gs > 0);
        assert(!(bs & 1));
    }

    proj Aed;
    fp2_add3(&Aed.z, &A->z, &A->z);   // tmp = 2C
    fp2_add3(&Aed.x, &A->x, &Aed.z);  // Aed.x = A + 2C
    fp2_sub3(&Aed.z, &A->x, &Aed.z);  // Aed.z = A - 2C

    fp2 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    fp2 Psum1, Pdif1, Psum2, Pdif2;
    fp2 coeffQ1[4], coeffQ2[4];
    fp2_add3(&Psum1, &points[0]->x, &points[0]->z);   //precomputations
    fp2_sub3(&Pdif1, &points[0]->x, &points[0]->z);

    fp2_add3(&Psum2, &points[1]->x, &points[1]->z);   //precomputations
    fp2_sub3(&Pdif2, &points[1]->x, &points[1]->z);

    fp2_add3(&coeffQ1[0], &points[0]->x, &points[0]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ1[1], &points[0]->x, &points[0]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ1[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ1[1]); //(x-z)^2
    fp2_add3(&coeffQ1[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ1[2], &tmp0, &tmp1); // 


    fp2_add3(&coeffQ2[0], &points[1]->x, &points[1]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ2[1], &points[1]->x, &points[1]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ2[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ2[1]); //(x-z)^2
    fp2_add3(&coeffQ2[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ2[2], &tmp0, &tmp1); // 

    int Minit[k];
    proj M[k];

    for (long long s = 0; s < k; s++) {
        Minit[s] = 0;
    }
    M[1] = *K;
    Minit[1] = 1;

    xDBL(&M[2], K, A);
    Minit[2] = 1;

    if (sq4tvelu) {
        for (long long s = 3; s < k; ++s) {
            if (s & 1) {
                long long i = s / 2;  
                assert(s == 2*i + 1);

                if (i < bs) {
                    if (s == 3) {
                        assert(Minit[1]);
                        assert(Minit[2]);
                        xADD(&M[s],
                             &M[2],
                             &M[1],
                             &M[1]);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            } 
            else {
                long long i = s/2 - 1; 
                assert(s == 2*i + 2);

                if (i < (k-1)/2 - 2*bs*gs) {
                    if (s == 4) {
                        assert(Minit[2]);
                        xDBL(&M[s], &M[2], A);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            }
            if (bs > 0) {
                if (s == 2*bs) {
                    assert(Minit[bs-1]);
                    assert(Minit[bs+1]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[bs+1],
                         &M[bs-1],
                         &M[2]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 4*bs) {
                    assert(Minit[2*bs]);
                    xDBL(&M[s], &M[2*bs], A);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 6*bs) {
                    assert(Minit[2*bs]);
                    assert(Minit[4*bs]);
                    xADD(&M[s],
                         &M[4*bs],
                         &M[2*bs],
                         &M[2*bs]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s % (4*bs) == 2*bs) {
                    long long j = s / (4*bs);
                    assert(s == 2*bs*(2*j + 1));
                    if (j < gs) {
                        assert(Minit[s-4*bs]);
                        assert(Minit[s-8*bs]);
                        assert(Minit[4*bs]);
                        xADD(&M[s],
                             &M[s-4*bs],
                             &M[4*bs],
                             &M[s-8*bs]);
                        Minit[s] = 1;
                        continue;
                    }
                }
            }
        }
    }
    else {
        for (long long i = 3; i <= (k-1)/2; i++) {
            Minit[i] = 1;
            xADD(&M[i],
                 &M[i-1],
                 K,
                 &M[i-2]);
        }
    }

    proj Abatch;
    Abatch.x = fp2_1; 
    Abatch.z = fp2_1;

    proj Q1, Q2;
    Q1.x = fp2_1;
    Q1.z = fp2_1;

    Q2.x = fp2_1;
    Q2.z = fp2_1;

    if (sq4tvelu) {
        long long TIlen = 2*bs + poly_tree1size(bs);
        fp2 TI[TIlen];
        for (long long i = 0; i < bs; ++i) {
            assert(Minit[2*i+1]);
            fp2_neg2(&TI[2*i], &M[2*i+1].x);
            TI[2*i+1] = M[2*i+1].z;
        }
        poly_tree1(TI + 2*bs, TI, bs);

        fp2 TP1[3*gs];
        fp2 TPinv1[3*gs];
        fp2 TP2[3*gs];
        fp2 TPinv2[3*gs];
        fp2 T1[3*gs];
        fp2 Tminus1[3*gs];
        fp2 coeff[5];        

        for (long long j = 0;j < gs;++j) {
            assert(Minit[2*bs*(2*j+1)]);
            biquad_pm1_opt_fp2(coeff, T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],&A->x,&A->z);
            biquad_both_fp2(TP1+3*j,TPinv1+3*j, coeff, coeffQ1, &A->z);
            biquad_both_fp2(TP2+3*j,TPinv2+3*j, coeff, coeffQ2, &A->z);
        }

        poly_multiprod2(TP1,gs);
        poly_multiprod2(TPinv1,gs);


        poly_multiprod2(TP2,gs);
        poly_multiprod2(TPinv2,gs);


        poly_multiprod2_selfreciprocal(T1,gs);
        poly_multiprod2_selfreciprocal(Tminus1,gs);

        

        long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
        fp2 precomp[precompsize];
        poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

        fp2 v[bs];

        poly_multieval_postcompute(v,bs,TP1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.x,&v[i]);



        poly_multieval_postcompute(v,bs,TP2,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q2.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q2.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv2,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q2.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q2.x,&v[i]);

        poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.x,&v[i]);

        poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.z,&v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.z,&v[i]);

        for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);
        fp2_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
        fp2_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
        fp2_mul2(&Abatch.x,&tmp1);
        fp2_mul2(&Abatch.z,&tmp0);


        fp2_mul3(&tmp2, &tmp1, &Psum1);
        fp2_mul3(&tmp3, &tmp0, &Pdif1);


        fp2_add3(&tmp4, &tmp3, &tmp2);
        fp2_mul2(&Q1.x, &tmp4);
        fp2_sub3(&tmp4, &tmp3, &tmp2);
        fp2_mul2(&Q1.z, &tmp4);

        fp2_mul3(&tmp2, &tmp1, &Psum2);
        fp2_mul3(&tmp3, &tmp0, &Pdif2);


        fp2_add3(&tmp4, &tmp3, &tmp2);
        fp2_mul2(&Q2.x, &tmp4);
        fp2_sub3(&tmp4, &tmp3, &tmp2);
        fp2_mul2(&Q2.z, &tmp4);

        }
    } 
    else {
        for (long long i = 1;i <= (k-1)/2;++i) {
            assert(Minit[i]);
            fp2_sub3(&tmp1, &M[i].x, &M[i].z);
            fp2_add3(&tmp0, &M[i].x, &M[i].z);
            if (i > 1) {
                fp2_mul2(&Abatch.x,&tmp1);
                fp2_mul2(&Abatch.z,&tmp0);
            } else {
                Abatch.x = tmp1;
                Abatch.z = tmp0;
            }
            fp2_mul3(&tmp2, &tmp1, &Psum1); 
            fp2_mul3(&tmp3, &tmp0, &Pdif1);
            fp2_add3(&tmp4, &tmp3, &tmp2);

            fp2_mul3(&tmp5, &tmp1, &Psum2); 
            fp2_mul3(&tmp6, &tmp0, &Pdif2);
            fp2_add3(&tmp7, &tmp6, &tmp5);

            if (i > 1) {
                fp2_mul2(&Q1.x, &tmp4);
                fp2_mul2(&Q2.x, &tmp7);
            } else {
                Q1.x = tmp4;
                Q2.x = tmp7;
            }
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_sub3(&tmp7, &tmp6, &tmp5);
            if (i > 1) {
                fp2_mul2(&Q1.z, &tmp4);
                fp2_mul2(&Q2.z, &tmp7);
            } 
            else {
                Q1.z = tmp4;
                Q2.z = tmp7;
            }
        }
    }
    // point evaluation
    fp2_sqr1(&Q1.x);
    fp2_sqr1(&Q1.z);
    fp2_mul2(&points[1]->x, &Q1.x);
    fp2_mul2(&points[1]->z, &Q1.z);

    fp2_sqr1(&Q2.x);
    fp2_sqr1(&Q2.z);
    fp2_mul2(&points[2]->x, &Q2.x);
    fp2_mul2(&points[2]->z, &Q2.z);
  

    exp_by_squaring(&Aed.x, &Aed.z, k);

    fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); // ^8
    fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); // ^8

    fp2_mul2(&Aed.z, &Abatch.x);
    fp2_mul2(&Aed.x, &Abatch.z);

    fp2_add3(&A->x, &Aed.x, &Aed.z);  // A = Aed.x + Aed.z = 2A'
    fp2_sub3(&A->z, &Aed.x, &Aed.z);  // A = Aed.x - Aed.z = 4C'
    fp2_add2(&A->x, &A->x);          // A = 4A'
    
}

void xISOG_1pt(proj *A, proj **points, proj const *K, long long k)
{
    assert(k >= 3);
    assert(k % 2 == 1);

    int sq4tvelu = 0;
    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);

    if (bs) {
        sq4tvelu = 1;
        assert(bs > 0);
        assert(gs > 0);
        assert(!(bs & 1));
    }

    proj Aed;
    fp2_add3(&Aed.z, &A->z, &A->z);   // tmp = 2C
    fp2_add3(&Aed.x, &A->x, &Aed.z);  // Aed.x = A + 2C
    fp2_sub3(&Aed.z, &A->x, &Aed.z);  // Aed.z = A - 2C

    fp2 tmp0, tmp1, tmp2, tmp3, tmp4;
    fp2 Psum1, Pdif1;
    fp2 coeffQ1[4];
    fp2_add3(&Psum1, &points[1]->x, &points[1]->z);   //precomputations
    fp2_sub3(&Pdif1, &points[1]->x, &points[1]->z);


    fp2_add3(&coeffQ1[0], &points[1]->x, &points[1]->z);   //precomputations //Psum=x+z
    fp2_sub3(&coeffQ1[1], &points[1]->x, &points[1]->z);   // Pdif= x-z
    fp2_sqr2(&tmp0, &coeffQ1[0]); //(x+z)^2
    fp2_sqr2(&tmp1, &coeffQ1[1]); //(x-z)^2
    fp2_add3(&coeffQ1[3], &tmp0, &tmp1); // 2(x^2+Z^2)
    fp2_sub3(&coeffQ1[2], &tmp0, &tmp1); // 
 

    int Minit[k];
    proj M[k];

    for (long long s = 0; s < k; s++) {
        Minit[s] = 0;
    }
    M[1] = *K;
    Minit[1] = 1;

    xDBL(&M[2], K, A);
    Minit[2] = 1;

    if (sq4tvelu) {
        for (long long s = 3; s < k; ++s) {
            if (s & 1) {
                long long i = s / 2;  
                assert(s == 2*i + 1);

                if (i < bs) {
                    if (s == 3) {
                        assert(Minit[1]);
                        assert(Minit[2]);
                        xADD(&M[s],
                             &M[2],
                             &M[1],
                             &M[1]);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            } 
            else {
                long long i = s/2 - 1; 
                assert(s == 2*i + 2);

                if (i < (k-1)/2 - 2*bs*gs) {
                    if (s == 4) {
                        assert(Minit[2]);
                        xDBL(&M[s], &M[2], A);
                        Minit[s] = 1;
                        continue;
                    }
                    assert(Minit[s-2]);
                    assert(Minit[s-4]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[s-2],
                         &M[2],
                         &M[s-4]);
                    Minit[s] = 1;
                    continue;
                }
            }
            if (bs > 0) {
                if (s == 2*bs) {
                    assert(Minit[bs-1]);
                    assert(Minit[bs+1]);
                    assert(Minit[2]);
                    xADD(&M[s],
                         &M[bs+1],
                         &M[bs-1],
                         &M[2]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 4*bs) {
                    assert(Minit[2*bs]);
                    xDBL(&M[s], &M[2*bs], A);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s == 6*bs) {
                    assert(Minit[2*bs]);
                    assert(Minit[4*bs]);
                    xADD(&M[s],
                         &M[4*bs],
                         &M[2*bs],
                         &M[2*bs]);
                    Minit[s] = 1;
                    continue;
                } 
                else if (s % (4*bs) == 2*bs) {
                    long long j = s / (4*bs);
                    assert(s == 2*bs*(2*j + 1));
                    if (j < gs) {
                        assert(Minit[s-4*bs]);
                        assert(Minit[s-8*bs]);
                        assert(Minit[4*bs]);
                        xADD(&M[s],
                             &M[s-4*bs],
                             &M[4*bs],
                             &M[s-8*bs]);
                        Minit[s] = 1;
                        continue;
                    }
                }
            }
        }
    }
    else {
        for (long long i = 3; i <= (k-1)/2; i++) {
            Minit[i] = 1;
            xADD(&M[i],
                 &M[i-1],
                 K,
                 &M[i-2]);
        }
    }

    proj Abatch;
    Abatch.x = fp2_1; 
    Abatch.z = fp2_1;

    proj Q1;
    Q1.x = fp2_1;
    Q1.z = fp2_1;

    if (sq4tvelu) {
        long long TIlen = 2*bs + poly_tree1size(bs);
        fp2 TI[TIlen];
        for (long long i = 0; i < bs; ++i) {
            assert(Minit[2*i+1]);
            fp2_neg2(&TI[2*i], &M[2*i+1].x);
            fp2_copy(&TI[2*i+1], &M[2*i+1].z);
        }
        poly_tree1(TI + 2*bs, TI, bs);

        fp2 TP1[3*gs];
        fp2 TPinv1[3*gs];
        fp2 T1[3*gs];
        fp2 Tminus1[3*gs];
        fp2 coeff[5];        

        for (long long j = 0;j < gs;++j) {
            assert(Minit[2*bs*(2*j+1)]);
            biquad_pm1_opt_fp2(coeff, T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],&A->x,&A->z);
            biquad_both_fp2(TP1+3*j,TPinv1+3*j, coeff, coeffQ1, &A->z);

        }

        poly_multiprod2(TP1,gs);
        poly_multiprod2(TPinv1,gs);


        poly_multiprod2_selfreciprocal(T1,gs);
        poly_multiprod2_selfreciprocal(Tminus1,gs);

        

        long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
        fp2 precomp[precompsize];
        poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

        fp2 v[bs];

        poly_multieval_postcompute(v,bs,TP1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.z, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.z,&v[i]);

        poly_multieval_postcompute(v,bs,TPinv1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Q1.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Q1.x,&v[i]);


        poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.x, &v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.x,&v[i]);

        poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
        fp2_copy(&Abatch.z,&v[0]);
        for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.z,&v[i]);

        for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
            assert(Minit[2*i+2]);
            fp2_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
            fp2_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
            fp2_mul2(&Abatch.x,&tmp1);
            fp2_mul2(&Abatch.z,&tmp0);


            fp2_mul3(&tmp2, &tmp1, &Psum1);
            fp2_mul3(&tmp3, &tmp0, &Pdif1);


            fp2_add3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q1.x, &tmp4);
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            fp2_mul2(&Q1.z, &tmp4);

        }
    } 
    else {
        for (long long i = 1;i <= (k-1)/2;++i) {
            assert(Minit[i]);
            fp2_sub3(&tmp1, &M[i].x, &M[i].z);
            fp2_add3(&tmp0, &M[i].x, &M[i].z);
            if (i > 1) {
                fp2_mul2(&Abatch.x,&tmp1);
                fp2_mul2(&Abatch.z,&tmp0);
            } else {
                Abatch.x = tmp1;
                Abatch.z = tmp0;
            }
            fp2_mul3(&tmp2, &tmp1, &Psum1); 
            fp2_mul3(&tmp3, &tmp0, &Pdif1);
            fp2_add3(&tmp4, &tmp3, &tmp2);


            if (i > 1) {
                fp2_mul2(&Q1.x, &tmp4);
            } else {
                Q1.x = tmp4;
            }
            fp2_sub3(&tmp4, &tmp3, &tmp2);
            if (i > 1) {
                fp2_mul2(&Q1.z, &tmp4);
            } 
            else {
                Q1.z = tmp4;
            }
        }
    }
    // point evaluation
    fp2_sqr1(&Q1.x);
    fp2_sqr1(&Q1.z);
    fp2_mul2(&points[1]->x, &Q1.x);
    fp2_mul2(&points[1]->z, &Q1.z);

    exp_by_squaring(&Aed.x, &Aed.z, k);

    fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); fp2_sqr1(&Abatch.x); // ^8
    fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); fp2_sqr1(&Abatch.z); // ^8

    fp2_mul2(&Aed.z, &Abatch.x);
    fp2_mul2(&Aed.x, &Abatch.z);

    fp2_add3(&A->x, &Aed.x, &Aed.z);  // A = Aed.x + Aed.z = 2A'
    fp2_sub3(&A->z, &Aed.x, &Aed.z);  // A = Aed.x - Aed.z = 4C'
    fp2_add2(&A->x, &A->x);          // A = 4A'
    
}

void xISOG_origin(proj *A, proj *P, proj *K, uint64_t l, bool want_multiple)
{   

    assert(l >= 3);
    assert(l % 2 == 1);

    if (isog_counters) ++isog_counters[l];

    proj prod, Q;

    montisog_eval_init(&Q, P, K);
    edwisog_curve_init(&prod, K);


    proj M[3] = {*K};
    xDBL(&M[1], K, A);

    for (uint64_t i = 1; i < l / 2; ++i) {

        if (i >= 2)
            xADD(&M[i%3], &M[(i-1)%3], K, &M[(i-2)%3]);

        montisog_eval_consume(&Q, P, &M[i%3]);
        edwisog_curve_consume(&prod, &M[i%3]);
    }

    /* "dummy isogeny" trick [Meyer-Campos-Reith] */
    if (want_multiple) {
        if (l > 3)
            xADD(&M[l/2%3], &M[(l/2-1)%3], K, &M[(l/2-2)%3]);
        xADD(K, &M[l/2%3], &M[(l/2-1)%3], K);
    }

    montisog_eval_finish(P, &Q);

    proj d;
    mont2edw(&d, A);
    edwisog_curve_finish(&d, &prod, l);
    edw2mont(A, &d);
}


bool is_twist(fp2 const *x, fp2 const *A)
{
    fp2 t;
    fp2_add3(&t, x, A);  /* x + A */
    fp2_mul2(&t, x);     /* x^2 + Ax */
    fp2_add2(&t, &fp2_1); /* x^2 + Ax + 1 */
    fp2_mul2(&t, x);     /* x^3 + Ax^2 + x */
    return !fp2_issquare(&t);
}

__attribute__((visibility("default")))
void j_inv(const fp2 *A, const fp2 *C, fp2 *jinv)
{   
    fp2 t0, t1;
    
    fp2_sqr2(jinv, A);
    fp2_sqr2(&t1, C);
    fp2_add3(&t0, &t1, &t1);
    fp2_sub3(&t0, jinv, &t0);
    fp2_sub3(&t0, &t0, &t1);
    fp2_sub3(jinv, &t0, &t1);
    fp2_sqr2(&t1, &t1);
    fp2_mul3(jinv, jinv, &t1);
    fp2_add3(&t0, &t0, &t0);
    fp2_add3(&t0, &t0, &t0);
    fp2_sqr2(&t1, &t0);
    fp2_mul3(&t0, &t1, &t0);
    fp2_add3(&t0, &t0, &t0);
    fp2_add3(&t0, &t0, &t0);
    fp2_inv(jinv);
    fp2_mul3(jinv, &t0, jinv);
}

void ec_recover_y(fp2 *y, fp2 const *x, proj const *A){

    fp2 rhs, tmp;

    fp2_sqr2(&rhs, x);
    fp2_mul3(&tmp, &A->x, &rhs); // tmp = Ax^2

    fp2_mul2(&rhs, x); // rhs = x^3

    fp2_add2(&rhs, &tmp); // rhs = x^3 + Ax^2
    fp2_add2(&rhs, x); // rhs = x^3 + Ax^2 + x

    fp2_sqrt(y, &rhs); // y = sqrt(x^3 + Ax^2 + x)

}

// output only x3 (P + Q = (x3, y3))
void affine_add(fp2 *x3, fp2 const *x1, fp2 const *y1, fp2 const *x2, fp2 const *y2, fp2 const *A){

    fp2 dx, dy, lam, tmp;

    fp2_sub3(&dx, x2, x1);
    fp2_sub3(&dy, y2, y1);
    fp2_inv(&dx);
    fp2_mul3(&lam, &dy, &dx); 
    fp2_sqr1(&lam); // lam^2 = {(y2-y1)/(x2-x1)}^2

    fp2_add3(&tmp, x1, x2);
    fp2_add2(&tmp, A);
    fp2_sub3(x3, &lam, &tmp); // x3 = lam^2 - A - x1 - x2

}

void ADD(proj *S, proj const *P, proj const *Q, proj const *A){

    fp2 x1, x2, x3, y1, y2;
    proj P_aff = *P, Q_aff = *Q;
    affinize(&P_aff, &Q_aff);

    x1 = P_aff.x; x2 = Q_aff.x;

    ec_recover_y(&y1, &x1, A); ec_recover_y(&y2, &x2, A);

    affine_add(&x3, &x1, &y1, &x2, &y2, &A->x);

    fp2_copy(&S->x, &x3); fp2_copy(&S->z, &fp2_1);

}

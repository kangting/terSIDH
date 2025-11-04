#ifndef MONT_H
#define MONT_H

#include "params.h"
#include "poly.h"
#include "constants.h"
#include "fpx.h"

bool is_infinity(proj const *P);
bool is_affine(proj const *P);
void affinize(proj *P, proj *Q);
void swap_points(proj* P, proj* Q, const uint64_t option);

void biquad_pm1_opt_fp2(fp2 *coeff, fp2 *outplus,fp2 *outminus,const proj *P,const fp2 *A, const fp2 *C);
void biquad_both_fp2(fp2 *out,fp2 *outinv, fp2 *coeff, fp2 *coeffQ, const fp2 *C);

void xDBL(proj *Q, proj const *P, proj const *A);
void xDBLv2(proj* Q, proj const* P, proj const* A24);

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ);
void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xMUL(proj *Q, proj const *A, proj const *P, uint_custom const *k);
// Fixed-iterations Montgomery ladder (constant-time wrt bitlength)
void xMUL_fixedbits(proj *Q, proj const *A, proj const *P, uint_custom const *k, size_t max_bits);

void eval_4_isog(proj *P, fp2 *K1, fp2 *K2, fp2 *K3);
void xISOG_4(proj *A, proj **points, proj *K, int num);

void xISOG_4pts(proj *A, proj **points, proj const *K, long long k);
void xISOG_3pts(proj *A, proj **points, proj const *K, long long k);
void xISOG_2pts(proj *A, proj **points, proj const *K, long long k);
void xISOG_1pt(proj *A, proj **points, proj const *K, long long k);

void xISOG_origin(proj *A, proj *P, proj *K, uint64_t l, bool want_multiple);

bool is_twist(fp2 const *x, fp2 const *A);

void j_inv(const fp2 *A, const fp2 *C, fp2 *jinv);

void ec_recover_y(fp2 *y, fp2 const *x, proj const *A);
void affine_add(fp2 *x3, fp2 const *x1, fp2 const *y1, fp2 const *x2, fp2 const *y2, fp2 const *A);
void ADD(proj *S, proj const *P, proj const *Q, proj const *A);


#endif


#include <string.h>
#include <assert.h>

#include "uint_custom.h"
#include "fpx.h"
#include "mont.h"
#include "tersidh.h"
#include "rng.h"
#include "params.h"
#include "constants.h"
#include "steps.h"

#include <time.h>
#include <stdio.h>

#define BATCH_SIZE 93

__attribute__((visibility("default")))
/* A = 6 */
const public_key base = {
    {{{6}}, {{0}}},
    {{{{0}}, {{0}}},{{{0}}, {{0}}}}, {{{{0}}, {{0}}},{{{0}}, {{0}}}}
};

__attribute__((visibility("default")))
void setup(proj **points, bool Alice)
{
    proj A;
    fp2_enc(&A.x, &base.A);
    A.z = fp2_1;

    if (Alice) {
        fp2_enc(&points[0]->x, &PA);
        fp2_enc(&points[1]->x, &QA);
        fp2_enc(&points[2]->x, &PB);
        fp2_enc(&points[3]->x, &QB);
    } else {
        fp2_enc(&points[0]->x, &PB);
        fp2_enc(&points[1]->x, &QB);
        fp2_enc(&points[2]->x, &PA);
        fp2_enc(&points[3]->x, &QA);
    }
    points[0]->z = points[1]->z = points[2]->z = points[3]->z = fp2_1;
    return;
}

__attribute__((visibility("default")))
void tersidh_private(private_key *priv)
{
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; ; ) {
        uint8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            for (int k = 0; k < 4; ++k) {
                int8_t v = (buf[j] >> (k * 2)) & 0x3;
                if (v == 3) continue;
                v = v - 1;  // private key : -1, 0, 1

                priv->e[i++] = v;

                if (i >= NUM_PRIMES)
                    return;
            }
        }
    }
}

bool validate_basic(public_key const *in)
{
    /* make sure A < p */
    uint_custom dummy;
    if (!uint_custom_sub3(&dummy, &in->A.a, &p))
        return false;
    if (!uint_custom_sub3(&dummy, &in->A.b, &p))
        return false;

    /* make sure the curve is nonsingular: A != 2 */
    uint_custom pm2;
    uint_custom_set(&pm2, 2);
    if (uint_custom_eq(&in->A.a, &pm2) && uint_custom_eq(&in->A.b, &uint_custom_0))
        return false;

    /* make sure the curve is nonsingular: A != -2 */
    uint_custom_sub3(&pm2, &uint_custom_0, &pm2);
    if (uint_custom_eq(&in->A.a, &pm2) && uint_custom_eq(&in->A.b, &uint_custom_0))
        return false;

    return true;
}

void isogeny(uint8_t number, proj *A, proj **points, int8_t es[NUM_PRIMES], const unsigned primes[NUM_PRIMES]) {

    uint8_t kernel_index = 0;

    for (bool twist = false; ; ) {
        #define BATCH_SIZE 93

        uint8_t batch[(NUM_PRIMES+7)/8] = {0};
        size_t sz = 0;

        for (size_t i = 0; i < NUM_PRIMES && sz < BATCH_SIZE; ++i) {
            if (twist ? (es[i] > 0) : (es[i] < 0)) {
                batch[i/8] |= 1<<i%8;
                ++sz;
            }
        }

        if (!sz) {
            if (twist) break;
            twist = true;
            kernel_index++;
            number--;
            continue;
        }

        uint_custom cofactor_zero = p_cofactor;
        for (size_t i = 0; i < NUM_PRIMES; ++i)
            if (es[i] == 0)
                uint_custom_mul3_64(&cofactor_zero, &cofactor_zero, primes[i]);

        uint_custom k = cofactor_zero;
        for (size_t i = 0; i < NUM_PRIMES; ++i)
            if ((~batch[i/8] & 1<<i%8) && es[i] != 0)
                uint_custom_mul3_64(&k, &k, primes[i]);

        proj K;
        xMUL(points[kernel_index], A, points[kernel_index], &k);

        for (size_t i = 0; i < NUM_PRIMES; ++i) {
            if (~batch[i/8] & (1 << (i % 8))) continue;

            uint_custom cof = uint_custom_1;

            for (size_t j = i+1; j < NUM_PRIMES; ++j)
                if (batch[j/8] & (1 << (j % 8)))
                    uint_custom_mul3_64(&cof, &cof, primes[j]);

            if (uint_custom_len(&cof) > (cost_ratio_inv_mul >> !is_affine(A)))
                affinize(A, &K);

            xMUL(&K, A, points[kernel_index], &cof);

            if(primes[i] == 4) {
                xISOG_4(A, points, &K, number);
                es[i] -= twist ? 1 : -1;
                continue;
            }

            if(number == 4) xISOG_4pts(A, points, &K, primes[i]);
            else if(number == 3) xISOG_3pts(A, points, &K, primes[i]);
            else if(number == 2) xISOG_2pts(A, points, &K, primes[i]);
            else xISOG_1pt(A, points, &K, primes[i]);

            es[i] -= twist ? 1 : -1;
        }

        assert(!is_infinity(A));
    }
}

void isogeny_reduced(uint8_t number, proj *A, proj **points, int8_t es[NUM_PRIMES], const unsigned primes[NUM_PRIMES]) {

    uint8_t batch[(NUM_PRIMES+7)/8] = {0};
    size_t sz = 0;

    for (size_t i = 0; i < NUM_PRIMES && sz < BATCH_SIZE; ++i) {
        if (es[i] != 0) {
            batch[i/8] |= 1<<i%8;
            ++sz;
        }
    }

    uint_custom cofactor0 = p_cofactor;
    uint_custom cofactor1 = p_cofactor;
    uint_custom cofactor2 = p_cofactor;

    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        if (es[i] > 0) uint_custom_mul3_64(&cofactor2, &cofactor2, primes[i]);
        else if (es[i] < 0) uint_custom_mul3_64(&cofactor1, &cofactor1, primes[i]);
        else uint_custom_mul3_64(&cofactor0, &cofactor0, primes[i]);
    }

    proj add, K, tmp1, tmp2;
    xMUL(&tmp1, A, points[0], &cofactor1);
    xMUL(&tmp2, A, points[1], &cofactor2);
    ADD(&add, &tmp1, &tmp2, A);
    xMUL(&add, A, &add, &cofactor0);

    proj *eval_points[4] = {&tmp1, &add, points[2], points[3]};

    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        if (~batch[i/8] & (1 << (i % 8))) continue;

        uint_custom cof = uint_custom_1;

        for (size_t j = i+1; j < NUM_PRIMES; ++j)
            if (batch[j/8] & (1 << (j % 8)))
                uint_custom_mul3_64(&cof, &cof, primes[j]);

        if (uint_custom_len(&cof) > (cost_ratio_inv_mul >> !is_affine(A)))
            affinize(A, &K);

        xMUL(&K, A, eval_points[1], &cof);

        if(primes[i] == 4) {
            xISOG_4(A, eval_points, &K, number);
            continue;
        }
        
        if(number == 3) xISOG_3pts(A, eval_points, &K, primes[i]);
        else xISOG_1pt(A, eval_points, &K, primes[i]);
        
    }

    assert(!is_infinity(A));

}

// 0/1 -> 0x00/0xFF
static inline uint8_t ct_mask_u8(int b) {
    return (uint8_t)(0 - (uint8_t)(b != 0));
}

// dst = m ? src : dst
static inline void ct_mem_cmov(void *dst, const void *src, size_t len, uint8_t m) {
    uint8_t *d = (uint8_t *)dst;
    const uint8_t *s = (const uint8_t *)src;
    for (size_t i = 0; i < len; i++) {
        uint8_t x = (uint8_t)(d[i] ^ s[i]);
        x &= m;
        d[i] ^= x;
    }
}

static inline void proj_cmov(proj *dst, const proj *src, uint8_t m) {
    ct_mem_cmov(&dst->x, &src->x, sizeof(dst->x), m);
    ct_mem_cmov(&dst->z, &src->z, sizeof(dst->z), m);
}

static inline void print_uint_custom(uint_custom x) {
    for (int i = 0; i < LIMBS; i++) {
        printf("0x%016lx, ", x.c[i]);
    }
    printf("\n");
}

void isogeny_constant(uint8_t number,
                      proj *A_out,
                      proj **points,
                      int8_t es[NUM_PRIMES],
                      const unsigned primes[NUM_PRIMES])
{
    unsigned char sel_real[(NUM_PRIMES + 7)/8] = {0};  // es[i]!=0
    unsigned char sel_dummy[(NUM_PRIMES + 7)/8] = {0}; // es[i]==0
    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        if (es[i] != 0) sel_real[i>>3] |= (1u << (i & 7));
        else            sel_dummy[i>>3] |= (1u << (i & 7));
    }

    // suf[i] = ∏_{j>=i} primes[j]
    uint_custom suf[NUM_PRIMES + 1];
    suf[NUM_PRIMES] = uint_custom_1;
    for (size_t idx = NUM_PRIMES; idx-- > 0; ) {
        suf[idx] = suf[idx+1];
        uint_custom_mul3_64(&suf[idx], &suf[idx], primes[idx]);
    }

    // ------------------------------------------------------------------
    // add_sel/add_unsel
    //  - c0 = ∏_{es==0} primes * p_cofactor
    //  - c1 = ∏_{es< 0} primes * p_cofactor
    //  - c2 = ∏_{es> 0} primes * p_cofactor
    // ------------------------------------------------------------------
    proj A_sel = *A_out;
    proj t1, t2, add_sel;

    uint_custom c0 = p_cofactor, c1 = p_cofactor, c2 = p_cofactor;
    for (size_t j = 0; j < NUM_PRIMES; ++j) {
        if      (es[j] < 0) uint_custom_mul3_64(&c1, &c1, primes[j]);
        else if (es[j] > 0) uint_custom_mul3_64(&c2, &c2, primes[j]);
        else                uint_custom_mul3_64(&c0, &c0, primes[j]);
    }

    xMUL(&t1, &A_sel, points[0], &c1);
    xMUL(&t2, &A_sel, points[1], &c2);
    ADD(&add_sel, &t1, &t2, &A_sel);
    xMUL(&add_sel, &A_sel, &add_sel, &c0);

    proj A_unsel = *A_out;
    proj add_unsel;
    xMUL(&add_unsel, &A_unsel, points[1], &c1);
    xMUL(&add_unsel, &A_unsel, &add_unsel, &c2);

    proj R_sel, S_sel, R_unsel, S_unsel; 
    if (number == 3){
        R_sel = *points[2]; S_sel = *points[3];
        R_unsel = *points[2]; S_unsel = *points[3];
    } 

    proj K2;

    for (size_t i = 0; i < NUM_PRIMES; ++i) {

        int M = (sel_real[i>>3] >> (i & 7)) & 1;
        uint8_t m  = ct_mask_u8(M);        // 0xFF if real, else 0x00
        uint8_t m0 = (uint8_t)~m;

        proj A = A_unsel, K = add_unsel;
        proj R = (proj){fp2_1,fp2_1}, S = (proj){fp2_1,fp2_1};
        if (number == 3) { R = R_unsel; S = S_unsel; }

        // Select inputs by cmov (real path)
        proj_cmov(&K, &add_sel, m); proj_cmov(&A, &A_sel, m);
        if (number == 3) { proj_cmov(&R, &R_sel, m); proj_cmov(&S, &S_sel, m); }

        // Real and dummy cofactors
        uint_custom cof_real = suf_real[i+1];
        uint_custom cof_dummy = suf_dummy[i+1];
        uint_custom cof_sel = cof_dummy; ct_mem_cmov(&cof_sel, &cof_real, sizeof(cof_sel), m);
        uint_custom cof_uns = cof_real;  ct_mem_cmov(&cof_uns, &cof_dummy, sizeof(cof_uns), m);

        // Real K2
        xMUL(&K2, &A, &K, &cof_sel);

        // Dummy-balanced xMUL on opposite path (discarded)
        proj A_op = A_sel, K_op = add_sel;
        proj_cmov(&A_op, &A_unsel, m); proj_cmov(&K_op, &add_unsel, m);

        proj *PP[4] = {&t1, &K, &R, &S};

        if (primes[i] == 4) xISOG_4(&A, PP, &K2, number);
        else {
            if      (number == 3) xISOG_3pts(&A, PP, &K2, primes[i]);
            else                  xISOG_1pt (&A, PP, &K2, primes[i]);
        }

        proj_cmov(&A_sel, &A, m); proj_cmov(&add_sel, &K, m);
        proj_cmov(&A_unsel, &A, m0); proj_cmov(&add_unsel, &K, m0);

        if (number == 3){
            proj_cmov(&R_sel, &R, m); proj_cmov(&S_sel, &S, m);
            proj_cmov(&R_unsel,&R, m0);    proj_cmov(&S_unsel,&S, m0);
        } 

    }

    *A_out = A_sel;
    if (number == 3) *points[2] = R_sel, *points[3] = S_sel;

    assert(!is_infinity(A_out));
}

/* includes public-key validation. */
/* totally not constant-time. */
__attribute__((visibility("default")))
bool keygen(public_key *out, proj **points, public_key const *in, private_key *priv, bool Alice)
{   
    int8_t es[NUM_PRIMES];
    memcpy(es, priv->e, sizeof(es));
    proj A;
    fp2_enc(&A.x, &in->A);
    A.z = fp2_1;

    // check Alice or Bob
    const unsigned *primes = Alice ? A_primes : B_primes;

    // isogeny computation
    // isogeny(4, &A, points, es, primes);
    // isogeny_reduced(3, &A, points, es, primes);
    isogeny_constant(3, &A, points, es, primes);

    affinize(&A, NULL);
    assert(is_affine(&A));

    fp2_dec(&out->A, &A.x);

    copy_point(&out->xR, points[2]);
    copy_point(&out->xS, points[3]);

    return true;
}

__attribute__((visibility("default")))
bool shared(fp2 *out, public_key *in, private_key const *priv, bool Alice){

    if (!validate_basic(in)){
        printf("invalid");
        return false;
    }

    int8_t es[NUM_PRIMES];
    memcpy(es, priv->e, sizeof(es));

    proj A;
    fp2_enc(&A.x, &in->A);
    A.z = fp2_1;

    proj *points[2] = {&in->xR, &in->xS};

    const unsigned *primes = Alice ? A_primes : B_primes;

    // isogeny computation for R
    // isogeny(2, &A, points, es, primes);
    // isogeny_reduced(1, &A, points, es, primes);
    isogeny_constant(1, &A, points, es, primes);

    j_inv(&A.x, &A.z, out);

    return true;
}

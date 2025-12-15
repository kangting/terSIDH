
#include <string.h>
#include <assert.h>
#include <stdlib.h>

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

typedef struct {
    uint64_t prime;
    unsigned exp;
} prime_power;

typedef struct {
    uint64_t m;             /* small modulus (prime power) */
    uint64_t inv_n_mod_m;   /* (modulus/m)^(-1) mod m */
    uint_custom n;          /* modulus/m */
} crt_piece;

#define MAX_PRIME_FACTORS (NUM_PRIMES * 2)

static prime_power factors_A[MAX_PRIME_FACTORS];
static prime_power factors_B[MAX_PRIME_FACTORS];
static size_t factors_A_len = 0, factors_B_len = 0;
static uint_custom order_A, order_B;
static uint_custom cofactor_A, cofactor_B;
static bool torsion_data_ready = false;
static crt_piece crt_A[MAX_PRIME_FACTORS], crt_B[MAX_PRIME_FACTORS];

static proj gen_PA, gen_QA, gen_PB, gen_QB, gen_PQA, gen_PQB;
static bool generators_ready = false;
static uint_custom active_cofactor;
static bool last_generators_random = false;

/* forward */
static void uint_custom_divmod_64_local(uint_custom *q, uint64_t *r, const uint_custom *x, uint64_t y);

static bool divides_small(uint_custom const *x, uint64_t p) {
    uint_custom q;
    uint64_t rem = 0;
    uint_custom_divmod_64_local(&q, &rem, x, p);
    return rem == 0;
}

static bool coprime_to_factors(uint_custom const *x, const prime_power *factors, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        if (divides_small(x, factors[i].prime))
            return false;
    }
    return true;
}

static size_t collect_prime_powers(prime_power *out, const unsigned *primes) {
    size_t count = 0;
    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        unsigned n = primes[i];
        for (unsigned p = 2; n > 1; ++p) {
            if (n % p) continue;
            unsigned exp = 0;
            while (n % p == 0) {
                n /= p;
                ++exp;
            }
            size_t idx = 0;
            for (; idx < count; ++idx) {
                if (out[idx].prime == p) {
                    out[idx].exp += exp;
                    break;
                }
            }
            if (idx == count) {
                assert(count < MAX_PRIME_FACTORS);
                out[count++] = (prime_power){p, exp};
            }
        }
    }
    return count;
}

static uint64_t prime_power_value(const prime_power *pp) {
    uint64_t v = 1;
    for (unsigned e = 0; e < pp->exp; ++e)
        v *= pp->prime;
    return v;
}

static int uint_custom_cmp(const uint_custom *x, const uint_custom *y) {
    for (size_t i = LIMBS; i-- > 0; ) {
        if (x->c[i] < y->c[i]) return -1;
        if (x->c[i] > y->c[i]) return 1;
    }
    return 0;
}

static void uint_custom_add_mod(uint_custom *out, const uint_custom *a, const uint_custom *b, const uint_custom *mod) {
    uint_custom sum;
    uint_custom_add3(&sum, a, b);
    if (uint_custom_cmp(&sum, mod) >= 0)
        uint_custom_sub3(&sum, &sum, mod);
    *out = sum;
}

static uint64_t inv_mod_u64(uint64_t a, uint64_t mod) {
    int64_t t0 = 0, t1 = 1;
    int64_t r0 = (int64_t)mod, r1 = (int64_t)(a % mod);
    while (r1 != 0) {
        int64_t q = r0 / r1;
        int64_t tmp = r0 - q * r1; r0 = r1; r1 = tmp;
        tmp = t0 - q * t1; t0 = t1; t1 = tmp;
    }
    if (r0 != 1) return 0;
    if (t0 < 0) t0 += mod;
    return (uint64_t)t0;
}

static void fill_crt_precomp(crt_piece *pcs, const prime_power *factors, size_t len,
                             const uint_custom *modulus) {
    for (size_t i = 0; i < len; ++i) {
        uint64_t m = prime_power_value(&factors[i]);
        pcs[i].m = m;

        uint_custom q;
        uint64_t rem = 0;
        uint_custom_divmod_64_local(&q, &rem, modulus, m); /* q = modulus / m */
        assert(rem == 0);
        pcs[i].n = q;

        uint64_t n_mod;
        uint_custom tmp;
        uint_custom_divmod_64_local(&tmp, &n_mod, &q, m);
        pcs[i].inv_n_mod_m = inv_mod_u64(n_mod % m, m);
        assert(pcs[i].inv_n_mod_m);
    }
}

static bool uint_custom_inv_crt_fast64(uint_custom *inv, uint64_t k_small,
                                       const uint_custom *modulus,
                                       const crt_piece *pcs, size_t pcs_len) {
    uint_custom acc = uint_custom_0;
    for (size_t i = 0; i < pcs_len; ++i) {
        uint64_t m = pcs[i].m;
        uint64_t k_mod = k_small % m;
        if (k_mod == 0)
            return false;
        uint64_t inv_k = inv_mod_u64(k_mod, m);
        if (!inv_k)
            return false;
        uint64_t coeff = (uint64_t)(((__uint128_t)inv_k * pcs[i].inv_n_mod_m) % m);

        uint_custom term;
        uint_custom_mul3_64(&term, &pcs[i].n, coeff);
        uint_custom_add_mod(&acc, &acc, &term, modulus);
    }
    *inv = acc;
    return true;
}

static void compute_orders_and_cofactors(void) {
    uint_custom_set(&order_A, 1);
    uint_custom_set(&order_B, 1);

    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        uint_custom_mul3_64(&order_A, &order_A, A_primes[i]);
        uint_custom_mul3_64(&order_B, &order_B, B_primes[i]);
    }

    for (size_t i = 0; i < LIMBS; ++i) {
        cofactor_A.c[i] = active_cofactor.c[i];
        cofactor_B.c[i] = active_cofactor.c[i];
    }
    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        uint_custom_mul3_64(&cofactor_A, &cofactor_A, B_primes[i]);
        uint_custom_mul3_64(&cofactor_B, &cofactor_B, A_primes[i]);
    }

    factors_A_len = collect_prime_powers(factors_A, A_primes);
    factors_B_len = collect_prime_powers(factors_B, B_primes);
    fill_crt_precomp(crt_A, factors_A, factors_A_len, &order_A);
    fill_crt_precomp(crt_B, factors_B, factors_B_len, &order_B);

    torsion_data_ready = true;
}

static bool random_point_on_curve(proj *P, const proj *A) {
    fp2 x, rhs, tmp;
    for (size_t i = 0; i < 256; ++i) {
        fp2_random(&x);

        fp2_sqr2(&rhs, &x);            /* x^2 */
        fp2_mul3(&tmp, &A->x, &rhs);   /* A*x^2 */
        fp2_mul2(&rhs, &x);            /* x^3 */
        fp2_add2(&rhs, &tmp);          /* x^3 + A*x^2 */
        fp2_add2(&rhs, &x);            /* x^3 + A*x^2 + x */

        fp2 rhs_check = rhs;
        if (!fp2_issquare(&rhs_check))
            continue;

        fp2_copy(&P->x, &x);
        P->z = fp2_1;
        return true;
    }
    return false;
}

static void uint_custom_divmod_64_local(uint_custom *q, uint64_t *r, const uint_custom *x, uint64_t y) {
    uint64_t rem = 0;
    for (size_t idx = 0; idx < LIMBS; ++idx) q->c[idx] = 0;

    /* Bit-at-a-time long division to stay within portable 64-bit ops */
    for (size_t limb = LIMBS; limb-- > 0; ) {
        uint64_t word = x->c[limb];
        for (unsigned bit = 0; bit < 64; ++bit) {
            rem = (rem << 1) | (word >> 63);
            word <<= 1;
            if (rem >= y) {
                rem -= y;
                q->c[limb] |= (uint64_t)1 << (63 - bit);
            }
        }
    }
    *r = rem;
}

static bool divide_by_small(uint_custom *x, uint64_t p) {
    uint_custom q;
    uint64_t rem = 0;
    uint_custom_divmod_64_local(&q, &rem, x, p);
    if (rem) return false;
    *x = q;
    return true;
}

static bool has_exact_order(const proj *P, const proj *A, const uint_custom *order,
                            const prime_power *factors, size_t factors_len) {
    proj tmp;

    /* ensure P*order == O */
    xMUL(&tmp, A, P, order);
    if (!is_infinity(&tmp))
        return false;

    /* check removing any prime factor (with multiplicity) kills the exactness */
    for (size_t i = 0; i < factors_len; ++i) {
        uint_custom k = *order;
        for (unsigned e = 0; e < factors[i].exp; ++e) {
            if (!divide_by_small(&k, factors[i].prime))
                return false;
            xMUL(&tmp, A, P, &k);
            if (is_infinity(&tmp))
                return false;
        }
    }
    return true;
}

static bool sample_point_of_order(proj *out, const proj *A, const uint_custom *cofactor,
                                  const uint_custom *order, const prime_power *factors,
                                  size_t factors_len) {
    const clock_t budget = (clock_t)(5 * CLOCKS_PER_SEC); /* ~5s per attempt */
    const clock_t t0 = clock();

    proj candidate;
    size_t attempts = 0;
    for (size_t i = 0; i < 1024; ++i) {
        if (clock() - t0 > budget)
            break;

        ++attempts;

        if (!random_point_on_curve(&candidate, A))
            continue;

        xMUL(&candidate, A, &candidate, cofactor);
        if (is_infinity(&candidate))
            continue;

        bool order_ok = has_exact_order(&candidate, A, order, factors, factors_len);

        if (!order_ok)
            continue;

        *out = candidate;
        return true;
    }
    return false;
}

static bool sample_mask(uint_custom *mask, uint_custom *mask_inv, bool Alice) {
    if (!torsion_data_ready)
        compute_orders_and_cofactors();

    const uint_custom *modulus = Alice ? &order_B : &order_A;
    const prime_power *factors = Alice ? factors_B : factors_A;
    size_t factors_len = Alice ? factors_B_len : factors_A_len;

    for (int tries = 0; tries < 64; ++tries) {
        uint64_t r;
        randombytes(&r, sizeof(r));
        r |= 1; /* avoid zero and even value that would conflict with 2^2 factor */

        bool ok = true;
        for (size_t i = 0; i < factors_len; ++i) {
            if ((r % prime_power_value(&factors[i])) == 0) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;

        uint_custom inv;
        if (!uint_custom_inv_crt_fast64(&inv, r, modulus,
                                        Alice ? crt_B : crt_A, factors_len))
            continue;

        uint_custom_set(mask, r);
        *mask_inv = inv;
        return true;
    }
    return false;
}

static void compute_x_difference(proj *out, const proj *P, const proj *Q, const proj *A) {
    proj P_aff = *P, Q_aff = *Q;
    affinize(&P_aff, &Q_aff);

    fp2 yP, yQ, x3, yQ_neg;
    ec_recover_y(&yP, &P_aff.x, A);
    ec_recover_y(&yQ, &Q_aff.x, A);
    fp2_neg2(&yQ_neg, &yQ);
    affine_add(&x3, &P_aff.x, &yP, &Q_aff.x, &yQ_neg, &A->x);

    fp2_copy(&out->x, &x3);
    out->z = fp2_1;
}

static bool ensure_generators(void) {
    if (generators_ready)
        return true;

    if (!torsion_data_ready)
        compute_orders_and_cofactors();

    const clock_t budget_total = (clock_t)(20 * CLOCKS_PER_SEC); /* stop after ~20s overall */
    const clock_t t_start = clock();

    proj A;
    fp2_enc(&A.x, &base.A);
    A.z = fp2_1;

    if (!sample_point_of_order(&gen_PA, &A, &cofactor_A, &order_A, factors_A, factors_A_len))
        return false;
    if (clock() - t_start > budget_total) return false;
    if (!sample_point_of_order(&gen_QA, &A, &cofactor_A, &order_A, factors_A, factors_A_len))
        return false;
    if (clock() - t_start > budget_total) return false;
    if (!sample_point_of_order(&gen_PB, &A, &cofactor_B, &order_B, factors_B, factors_B_len))
        return false;
    if (clock() - t_start > budget_total) return false;
    if (!sample_point_of_order(&gen_QB, &A, &cofactor_B, &order_B, factors_B, factors_B_len))
        return false;
    if (clock() - t_start > budget_total) return false;

    compute_x_difference(&gen_PQA, &gen_PA, &gen_QA, &A);
    compute_x_difference(&gen_PQB, &gen_PB, &gen_QB, &A);

    generators_ready = true;
    return true;
}

__attribute__((visibility("default")))
bool get_public_generators(proj *PA, proj *QA, proj *PB, proj *QB, proj *PQA, proj *PQB)
{
    if (!ensure_generators())
        return false;

    if (PA) copy_point(PA, &gen_PA);
    if (QA) copy_point(QA, &gen_QA);
    if (PB) copy_point(PB, &gen_PB);
    if (QB) copy_point(QB, &gen_QB);
    if (PQA) copy_point(PQA, &gen_PQA);
    if (PQB) copy_point(PQB, &gen_PQB);

    return true;
}

__attribute__((visibility("default")))
bool generators_were_random(void)
{
    return last_generators_random;
}

__attribute__((visibility("default")))
void setup(proj **points, bool Alice)
{
    proj A;
    fp2_enc(&A.x, &base.A);
    A.z = fp2_1;

    proj base_PA, base_QA, base_PB, base_QB;
    bool use_random = true;

    /* Start from baked generators and apply invertible scalars to each basis element. */
    active_cofactor = p_cofactor;
    torsion_data_ready = false;
    generators_ready = false;
    compute_orders_and_cofactors(); /* recompute with active_cofactor */

    fp2_enc(&base_PA.x, &PA); base_PA.z = fp2_1;
    fp2_enc(&base_QA.x, &QA); base_QA.z = fp2_1;
    fp2_enc(&base_PB.x, &PB); base_PB.z = fp2_1;
    fp2_enc(&base_QB.x, &QB); base_QB.z = fp2_1;

    uint_custom sPA = uint_custom_1, sQA = uint_custom_1, sPB = uint_custom_1, sQB = uint_custom_1;

    for (int tries = 0; tries < 32; ++tries) {
        uint_custom k; uint_custom_random(&k, &order_A);
        if (!uint_custom_eq(&k, &uint_custom_0) && coprime_to_factors(&k, factors_A, factors_A_len)) { sPA = k; break; }
    }
    for (int tries = 0; tries < 32; ++tries) {
        uint_custom k; uint_custom_random(&k, &order_A);
        if (!uint_custom_eq(&k, &uint_custom_0) && coprime_to_factors(&k, factors_A, factors_A_len)) { sQA = k; break; }
    }
    for (int tries = 0; tries < 32; ++tries) {
        uint_custom k; uint_custom_random(&k, &order_B);
        if (!uint_custom_eq(&k, &uint_custom_0) && coprime_to_factors(&k, factors_B, factors_B_len)) { sPB = k; break; }
    }
    for (int tries = 0; tries < 32; ++tries) {
        uint_custom k; uint_custom_random(&k, &order_B);
        if (!uint_custom_eq(&k, &uint_custom_0) && coprime_to_factors(&k, factors_B, factors_B_len)) { sQB = k; break; }
    }

    xMUL(&base_PA, &A, &base_PA, &sPA);
    xMUL(&base_QA, &A, &base_QA, &sQA);
    xMUL(&base_PB, &A, &base_PB, &sPB);
    xMUL(&base_QB, &A, &base_QB, &sQB);

    last_generators_random = use_random;

    if (Alice) {
        copy_point(points[0], &base_PA);
        copy_point(points[1], &base_QA);
        copy_point(points[2], &base_PB);
        copy_point(points[3], &base_QB);
    } else {
        copy_point(points[0], &base_PB);
        copy_point(points[1], &base_QB);
        copy_point(points[2], &base_PA);
        copy_point(points[3], &base_QA);
    }

    /* one-line indicator for the caller/user */
    (void)last_generators_random; /* always scaled now */
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

        uint_custom cofactor_zero = active_cofactor;
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

    uint_custom cofactor0 = active_cofactor;
    uint_custom cofactor1 = active_cofactor;
    uint_custom cofactor2 = active_cofactor;

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
    for (size_t i = 0; i < NUM_PRIMES; ++i) {
        if (es[i] != 0) sel_real[i>>3] |= (1u << (i & 7));
    }

    // Randomize processing order of primes (Fisher-Yates shuffle)
    size_t order[NUM_PRIMES];
    for (size_t i = 0; i < NUM_PRIMES; ++i) order[i] = i;
    for (size_t k = NUM_PRIMES; k > 1; --k) {
        uint32_t r;
        randombytes(&r, sizeof(r));
        size_t j = (size_t)(r % k);
        size_t tmp = order[k-1]; order[k-1] = order[j]; order[j] = tmp;
    }
    // Track which primes have been used already (to exclude from future cofactors)
    uint8_t used[NUM_PRIMES] = {0};

    // ------------------------------------------------------------------
    // add_sel/add_unsel
    //  - c0 = ∏_{es==0} primes * cofactor
    //  - c1 = ∏_{es< 0} primes * cofactor
    //  - c2 = ∏_{es> 0} primes * cofactor
    // ------------------------------------------------------------------
    // Compute public max bitlength for prep xMULs: C_total = cofactor * ∏ primes
    uint_custom C_total = active_cofactor;
    for (size_t j = 0; j < NUM_PRIMES; ++j) {
        uint_custom_mul3_64(&C_total, &C_total, primes[j]);
    }
    size_t max_bits_prep = uint_custom_len(&C_total);
    proj A_sel = *A_out;
    proj t1, t2, add_sel;

    uint_custom c0 = active_cofactor, c1 = active_cofactor, c2 = active_cofactor;
    for (size_t j = 0; j < NUM_PRIMES; ++j) {
        if      (es[j] < 0) uint_custom_mul3_64(&c1, &c1, primes[j]);
        else if (es[j] > 0) uint_custom_mul3_64(&c2, &c2, primes[j]);
        else                uint_custom_mul3_64(&c0, &c0, primes[j]);
    }

    xMUL_fixedbits(&t1, &A_sel, points[0], &c1, max_bits_prep);
    xMUL_fixedbits(&t2, &A_sel, points[1], &c2, max_bits_prep);
    ADD(&add_sel, &t1, &t2, &A_sel);
    xMUL_fixedbits(&add_sel, &A_sel, &add_sel, &c0, max_bits_prep);

    proj A_unsel = *A_out;
    proj add_unsel;
    xMUL_fixedbits(&add_unsel, &A_unsel, points[1], &c1, max_bits_prep);
    xMUL_fixedbits(&add_unsel, &A_unsel, &add_unsel, &c2, max_bits_prep);

    proj R_sel, S_sel, R_unsel, S_unsel; 
    if (number == 3){
        R_sel = *points[2]; S_sel = *points[3];
        R_unsel = *points[2]; S_unsel = *points[3];
    } 
    proj K2;

    for (size_t it = 0; it < NUM_PRIMES; ++it) {
        size_t i = order[it];

        int M = (sel_real[i>>3] >> (i & 7)) & 1;
        uint8_t m  = ct_mask_u8(M);        // 0xFF if real, else 0x00
        uint8_t m0 = (uint8_t)~m;

        proj A = A_unsel, K = add_unsel;
        proj R = (proj){fp2_1,fp2_1}, S = (proj){fp2_1,fp2_1};
        if (number == 3) { R = R_unsel; S = S_unsel; }

        // Select inputs by cmov (real path)
        proj_cmov(&K, &add_sel, m); proj_cmov(&A, &A_sel, m);
        if (number == 3) { proj_cmov(&R, &R_sel, m); proj_cmov(&S, &S_sel, m); }

        // Build cofactors on-the-fly by multiplying remaining primes in each subgroup
        uint_custom cof_sel = uint_custom_1;
        uint_custom cof_uns = uint_custom_1;
        for (size_t j = 0; j < NUM_PRIMES; ++j) {
            if (j == i || used[j]) continue; // skip current and already used
            if ((es[j] != 0) == (M != 0)) {
                uint_custom_mul3_64(&cof_sel, &cof_sel, primes[j]);
            } else {
                uint_custom_mul3_64(&cof_uns, &cof_uns, primes[j]);
            }
        }

        // Real K2
        xMUL(&K2, &A, &K, &cof_sel);

        // Dummy-balanced xMUL on opposite path (discarded)
        proj A_op = A_sel, K_op = add_sel, K2_dummy;
        proj_cmov(&A_op, &A_unsel, m); proj_cmov(&K_op, &add_unsel, m);
        xMUL(&K2_dummy, &A_op, &K_op, &cof_uns); (void)K2_dummy;

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

        used[i] = 1; // mark this prime as processed
    }

    *A_out = A_sel;
    if (number == 3) *points[2] = R_sel, *points[3] = S_sel;

    assert(!is_infinity(A_out));
}

/* includes public-key validation. */
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

    uint_custom mask, mask_inv;
    if (!sample_mask(&mask, &mask_inv, Alice))
        return false;

    xMUL(points[2], &A, points[2], &mask);
    xMUL(points[3], &A, points[3], &mask_inv);

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

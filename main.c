
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include "fpx.h"
#include "tersidh.h"
#include "uint_custom.h"
#include "mont.h"

#define N 1

typedef struct {
    const char *label;
    void (*fill)(private_key *k);
} scenario;

static void key_fill_random(private_key *k){ tersidh_private(k); }
static void key_fill_all_pos(private_key *k){ k->e[0]=0; k->e[1]=1; k->e[2]=-1; for(int i=3;i<NUM_PRIMES;++i) k->e[i]= 1; }
static void key_fill_all_neg(private_key *k){ k->e[0]=0; k->e[1]=1; k->e[2]=-1; for(int i=3;i<NUM_PRIMES;++i) k->e[i]=-1; }
static void key_fill_all_zero(private_key *k){ k->e[0]=0; k->e[1]=1; k->e[2]=-1; for(int i=3;i<NUM_PRIMES;++i) k->e[i]= 0; }
static void key_fill_alt_01(private_key *k){ k->e[0]=0; k->e[1]=1; k->e[2]=-1; for(int i=3;i<NUM_PRIMES;++i) k->e[i]= (i%2==0)?1:0; }
static void key_fill_alt_02(private_key *k){ k->e[0]=0; k->e[1]=1; k->e[2]=-1; for(int i=3;i<NUM_PRIMES;++i) k->e[i]= (i%2==0)?-1:0; }
static void key_fill_alt_12(private_key *k){ k->e[0]=0; k->e[1]=1; k->e[2]=-1; for(int i=3;i<NUM_PRIMES;++i) k->e[i]= (i%2==0)?1:-1; }

static void key_count(private_key *k, int *neg, int *zero, int *pos){
    int n=0,z=0,p=0;
    for (int i=0;i<NUM_PRIMES;i++){ if(k->e[i]<0) n++; else if(k->e[i]>0) p++; else z++; }
    if(neg) *neg=n; 
    if(zero) *zero=z; 
    if(pos) *pos=p;
}

static void key_preview(private_key *k, char *out, size_t outlen){
    // print first 27 trits mapped: -1->2, 0->0, 1->1 (fits 30 with ...)
    const size_t PREV_DIGITS = 27;
    size_t L = (NUM_PRIMES<PREV_DIGITS)?NUM_PRIMES:PREV_DIGITS;
    if(outlen==0) return; 
    out[0]='\0';
    size_t pos=0;
    for (size_t i=0;i<L && pos+1<outlen;i++){
        int v = k->e[i]==-1?2:(k->e[i]==1?1:0);
        out[pos++] = (char)('0'+v);
    }
    if(L<NUM_PRIMES && pos+3<outlen){ out[pos++]='.'; out[pos++]='.'; out[pos++]='.'; }
    out[pos]='\0';
}

static void print_table_header(void){
    printf("\n");
    printf("+----+------------------------+----------+-----------+--------------------------------+\n");
    printf("| #  | Scenario               | neg/0/1  | ms        | key preview                    |\n");
    printf("+----+------------------------+----------+-----------+--------------------------------+\n");
}

static void print_table_row(int idx, const char *label, private_key *k, double ms){
    int n=0,z=0,p=0; key_count(k,&n,&z,&p);
    char prev[40]; key_preview(k, prev, sizeof(prev));
    printf("| %2d | %-22s | %2d/%2d/%2d | %9.3f | %-30s |\n",
           idx, label, n, z, p, ms, prev);
}
static void print_table_footer(void){
    printf("+----+------------------------+----------+----------+---------------------------------+\n");
}

void uint_print(uint_custom *x)
{
    for (size_t i = 8*LIMBS-1; i < 8*LIMBS; --i)
        printf("%02hhx", i[(unsigned char *) x->c]);
}

void priv_print(private_key *k)
{
    size_t count = sizeof(k->e) / sizeof(k->e[0]); // size of private key
    for (size_t i = 0; i < count; ++i) {
        // -1, 0, 1 -> 2, 0, 1
        int converted_value = (k->e[i] == -1) ? 2 : (k->e[i] == 1) ? 1 : 0;

        printf("%d", converted_value);

    }
    printf("\n");
}

int main(void)
{
    clock_t t0, t1;
    bool ret; (void) ret;
    
    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    proj base_alice_points[4], base_bob_points[4];
    proj work_alice_points[4], work_bob_points[4];
    fp2 shared_alice, shared_bob;
    proj *points_alice[4] = {&work_alice_points[0], &work_alice_points[1], &work_alice_points[2], &work_alice_points[3]};
    proj *points_bob[4] = {&work_bob_points[0], &work_bob_points[1], &work_bob_points[2], &work_bob_points[3]};

    double alice_keygen_sum = 0;
    double bob_keygen_sum = 0;
    double alice_shared_sum = 0;
    double bob_shared_sum = 0;

    // Pretty table per private key scenario
    scenario cases[] = {
        {"random",    key_fill_random},
        {"all +1",    key_fill_all_pos},
        {"all -1",    key_fill_all_neg},
        {"all 0",     key_fill_all_zero},
        {"alt 0/1",   key_fill_alt_01},
        {"alt 0/-1",  key_fill_alt_02},
        {"alt 1/-1",  key_fill_alt_12},
    };
    int S = (int)(sizeof(cases)/sizeof(cases[0]));

    /* Prepare base generator sets once, then copy per iteration to avoid expensive setup */
    proj *setup_targets_alice[4] = {&base_alice_points[0], &base_alice_points[1], &base_alice_points[2], &base_alice_points[3]};
    proj *setup_targets_bob[4]   = {&base_bob_points[0],   &base_bob_points[1],   &base_bob_points[2],   &base_bob_points[3]};
    setup(setup_targets_alice, true);
    setup(setup_targets_bob, false);

    print_table_header();
    for (int s = 0; s < S; ++s){
        double sum_ms = 0.0;
        for (int iter = 0; iter < N; ++iter) {
            /* refresh working copies */
            for (int i = 0; i < 4; ++i) copy_point(points_alice[i], &base_alice_points[i]);
            cases[s].fill(&priv_alice);
            t0 = clock();
            ret = keygen(&pub_alice, points_alice, &base, &priv_alice, true);
            assert(ret);
            t1 = clock();
            sum_ms += 1000.0 * (t1 - t0) / CLOCKS_PER_SEC;
        }
        double ms = sum_ms / (double)N;
        print_table_row(s+1, cases[s].label, &priv_alice, ms);
    }
    print_table_footer();
    printf("\n");

    // normal case (legacy output)

    for (int iter = 0; iter < N; iter++) {

        for (int i = 0; i < 4; ++i) {
            copy_point(points_alice[i], &base_alice_points[i]);
            copy_point(points_bob[i],   &base_bob_points[i]);
        }

        tersidh_private(&priv_alice);
        tersidh_private(&priv_bob);

        t0 = clock();
        ret = keygen(&pub_alice, points_alice, &base, &priv_alice, true);
        assert(ret);
        t1 = clock();
        alice_keygen_sum += 1000.0 * (t1 - t0) / CLOCKS_PER_SEC;

        t0 = clock();
        ret = keygen(&pub_bob, points_bob, &base, &priv_bob, false);
        assert(ret);
        t1 = clock();
        bob_keygen_sum += 1000.0 * (t1 - t0) / CLOCKS_PER_SEC;

        t0 = clock();
        ret = shared(&shared_alice, &pub_bob, &priv_alice, true);
        assert(ret);
        t1 = clock();
        alice_shared_sum += 1000.0 * (t1 - t0) / CLOCKS_PER_SEC;

        t0 = clock();
        ret = shared(&shared_bob, &pub_alice, &priv_bob, false);
        assert(ret);
        t1 = clock();
        bob_shared_sum += 1000.0 * (t1 - t0) / CLOCKS_PER_SEC;

        assert(fp2_eq(&shared_alice, &shared_bob));

    }

    printf("\n");
    printf("The entire process of the normal terSIDH");
    printf("\n\n");
    printf("Average alice keygen time: %7.3lf ms\n", alice_keygen_sum / N);
    printf("Average bob keygen time:   %7.3lf ms\n", bob_keygen_sum / N);
    printf("Average alice shared time: %7.3lf ms\n", alice_shared_sum / N);
    printf("Average bob shared time:   %7.3lf ms\n", bob_shared_sum / N);
    printf("\n");

    // print_fp2(shared_alice);
    // print_fp2(shared_bob);
    
}

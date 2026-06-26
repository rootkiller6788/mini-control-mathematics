/**
 * @file ex02_markov_chain.c
 * @brief Example: Discrete-time Markov chain analysis — weather model.
 *
 * Demonstrates: markov_dtmc_create(), markov_dtmc_stationary(),
 * markov_dtmc_nstep(), markov_dtmc_mfpt(), gillespie_simulate().
 *
 * Scenario: Three-state weather model (Sunny, Cloudy, Rainy).
 * Analyses long-run behaviour, n-step transition probabilities,
 * and mean first-passage times.
 *
 * Knowledge: L2 core concepts (Markov property, stationary distribution,
 * Chapman-Kolmogorov), L6 engineering problems (weather/risk modeling).
 */
#include "stochastic_process.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Example 2: Markov Chain Weather Model ===\n\n");

    /* Three-state weather: 0=Sunny, 1=Cloudy, 2=Rainy */
    /* Transition matrix (daily):
       From Sunny:  60% stay sunny, 30% cloudy, 10% rainy
       From Cloudy: 40% sunny, 40% cloudy, 20% rainy
       From Rainy:  25% sunny, 45% cloudy, 30% rainy */
    prob_t P[9] = {
        0.60, 0.30, 0.10,
        0.40, 0.40, 0.20,
        0.25, 0.45, 0.30
    };
    prob_t pi0[3] = {1.0, 0.0, 0.0}; /* Start sunny */

    markov_chain_dtmc_t *weather = markov_dtmc_create(P, pi0, 3);
    if (!weather) { printf("Failed to create DTMC\n"); return 1; }

    /* Stationary distribution */
    prob_t pi_stat[3];
    markov_dtmc_stationary(weather, pi_stat);
    printf("Stationary distribution (long-run weather):\n");
    printf("  Sunny:  %.1f%%\n", pi_stat[0] * 100.0);
    printf("  Cloudy: %.1f%%\n", pi_stat[1] * 100.0);
    printf("  Rainy:  %.1f%%\n\n", pi_stat[2] * 100.0);

    /* Detailed balance check */
    double db = markov_dtmc_detailed_balance(weather, pi_stat);
    printf("Detailed balance deviation: %.2e\n", db);
    printf("Chain is %s\n\n", (db < 1e-10) ? "REVERSIBLE" : "NOT reversible");

    /* N-step transitions */
    printf("N-step transition probabilities P^n(0, j):\n");
    printf("  n      Sunny    Cloudy   Rainy\n");
    int steps[] = {1, 2, 5, 10, 50};
    int i;
    for (i = 0; i < 5; i++) {
        prob_t Pn[9];
        markov_dtmc_nstep(weather, Pn, (uint32_t)steps[i]);
        printf("  %-5d  %.3f   %.3f   %.3f\n", steps[i],
               Pn[0*3+0], Pn[0*3+1], Pn[0*3+2]);
    }
    printf("  (stationary) %.3f   %.3f   %.3f\n\n",
           pi_stat[0], pi_stat[1], pi_stat[2]);

    /* Mean first-passage times */
    printf("Mean first-passage time (days):\n");
    const char *names[] = {"Sunny", "Cloudy", "Rainy"};
    int from, to;
    for (from = 0; from < 3; from++) {
        for (to = 0; to < 3; to++) {
            if (from == to) continue;
            double mfpt = markov_dtmc_mfpt(weather, (size_t)from, (size_t)to);
            printf("  %s -> %s: %.1f days\n", names[from], names[to], mfpt);
        }
    }

    /* Simulate 30-day trajectory */
    printf("\nSimulated 30-day weather trajectory:\n  Day: ");
    size_t states[31];
    markov_dtmc_simulate(weather, 30, states, 12345);
    const char sym[] = {'S', 'C', 'R'};
    for (i = 0; i <= 30; i++) {
        if (i > 0 && i % 10 == 0) printf("\n       ");
        printf("%c ", sym[states[i]]);
    }
    printf("\n  (S=Sunny, C=Cloudy, R=Rainy)\n");

    markov_dtmc_free(weather);

    /* Bonus: Gillespie simulation of a simple birth-death process */
    printf("\n--- Bonus: Gillespie SSA of birth-death process ---\n");
    prob_real_t Q_bd[9] = {
        -0.5, 0.5, 0.0,
         0.3,-0.8, 0.5,
         0.0, 0.4,-0.4
    };
    prob_real_t times[100];
    size_t bd_states[100];
    size_t nj = gillespie_simulate(Q_bd, 3, 0, 20.0, times, bd_states, 100, 42);
    printf("  Simulated %zu jumps in [0, 20]\n", nj);
    printf("  Final state: %zu\n", bd_states[nj-1]);

    return 0;
}

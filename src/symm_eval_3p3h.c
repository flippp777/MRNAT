#include "../include/insertionSort_perm.h"
#include "../include/numberOfStates.h"
#include "../include/permute.h"
#include "../include/printb.h"
#include "../include/symm_eval.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void symm_eval_3p3h(int states[], unsigned long long statesb[], int smax,
                    int **phi, int phib[], int kmax, int Mref) {
  int mu; // Index for SD mu
  int nu; // Index for SD nu
  int a;  // Index for particle a
  int r;  // Index for paritcle r
  int b;  // Index for hole b
  int s;  // Index for hole s

  // Values for comparisons
  int p3h = 0; // Number of 3p3h case
  int p3h_contr = 0;
  int p3h_pos = 0;
  int p3h_neg = 0;
  int p3h_steps = 0;
  int p2h = 0; // Number of 2p2h case
  int p2h_contr = 0;
  int p2h_pos = 0;
  int p2h_neg = 0;
  int p2h_steps = 0;
  int p1h = 0; // Number of 1p1h case
  int p1h_contr = 0;
  int p1h_pos = 0;
  int p1h_neg = 0;
  int p1h_steps = 0;
  int p0h = 0; // Number of 0p0h case
  int p0h_contr = 0;
  int p0h_pos = 0;
  int p0h_neg = 0;
  int p0h_steps = 0;

  for (mu = 0; mu < Mref; mu++) {
    for (nu = 0; nu < Mref; nu++) {
      int ph = phib[mu] ^ phib[nu];
      int ph_num; // Number #p#h excitations between SD's
      ph_num = numberOfStates(ph);

      // 3p3h difference
      if (ph_num == 6) {
        p3h += 1;

        int ph_states[6];
        int ph_statesb[6];
        int ph_states_idx = 0;
        int ph_count = 0;
        for (int i = 0; i < smax; i++) {
          p3h_steps += 1;
          if ((ph & 1) == 1) {
            ph_states[ph_states_idx] = states[i];
            ph_statesb[ph_states_idx] = statesb[i];
            ph_states_idx++;
            ph_count++;
          }
          if (ph_count == 6) {
            break;
          }
          ph >>= 1;
        }

        int p_states[3];
        int p_states_idx = 0;
        int h_states[3];
        int h_states_idx = 0;
        for (int i = 0; i < 6; i++) {
          p3h_steps += 1;
          if ((ph_statesb[i] & phib[nu]) == ph_statesb[i]) {
            h_states[h_states_idx] = ph_states[i];
            h_states_idx++;
          } else {
            p_states[p_states_idx] = ph_states[i];
            p_states_idx++;
          }
        }

        int perm_count = 6;
        int count = 0;

        int **h_perms = malloc(perm_count * sizeof(unsigned long *));
        for (int i = 0; i < perm_count; i++) {
          h_perms[i] = malloc(3 * sizeof(unsigned long));
        }
        permute(h_states, 0, 3, h_perms, &count);

        int **p_perms = malloc(perm_count * sizeof(unsigned long *));
        for (int i = 0; i < perm_count; i++) {
          p_perms[i] = malloc(3 * sizeof(unsigned long));
        }
        count = 0;
        permute(p_states, 0, 3, p_perms, &count);

        // phi_nu with its 3p-3h "excitations"
        for (int hi = 0; hi < perm_count; hi++) {
          p3h_steps++;
          for (int pi = 0; pi < perm_count; pi++) {
            p3h_steps++;
            int phi_nu[kmax];
            for (int k = 0; k < kmax; k++) {
              p3h_steps++;
              int replaced = 0;
              for (int idx = 0; idx < 3; idx++) {
                p3h_steps++;
                if (phi[nu][k] == h_perms[hi][idx]) {
                  phi_nu[k] = p_perms[pi][idx];
                  replaced = 1;
                  break;
                }
              }
              if (!replaced)
                phi_nu[k] = phi[nu][k];
            }

            // Track permutations and sorting for sign
            int perm = 0;
            insertionSort_perm_steps(phi_nu, kmax, &perm, &p3h_steps);
            int sign = pow(-1, perm);

            p3h_contr += 1;
            if (sign > 0) {
              p3h_pos += 1;
            } else {
              p3h_neg += 1;
            }
          }
        }
      }

      // 2p2h difference
      if (ph_num == 4) {
        p2h += 1;

        int ph_states[4];
        int ph_statesb[4];
        int ph_states_idx = 0;
        int ph_count = 0;
        for (int i = 0; i < smax; i++) {
          p2h_steps += 1;
          if ((ph & 1) == 1) {
            ph_states[ph_states_idx] = states[i];
            ph_statesb[ph_states_idx] = statesb[i];
            ph_states_idx++;
            ph_count++;
          }
          if (ph_count == 4) {
            break;
          }
          ph >>= 1;
        }

        int p_states[2];
        int p_states_idx = 0;
        int h_states[2];
        int h_states_idx = 0;
        for (int i = 0; i < 4; i++) {
          p2h_steps += 1;
          if ((ph_statesb[i] & phib[nu]) == ph_statesb[i]) {
            h_states[h_states_idx] = ph_states[i];
            h_states_idx++;
          } else {
            p_states[p_states_idx] = ph_states[i];
            p_states_idx++;
          }
        }

        int perm_count = 4;
        int count = 0;

        int **h_perms = malloc(perm_count * sizeof(unsigned long *));
        for (int i = 0; i < perm_count; i++) {
          h_perms[i] = malloc(2 * sizeof(unsigned long));
        }
        permute(h_states, 0, 2, h_perms, &count);

        int **p_perms = malloc(perm_count * sizeof(unsigned long *));
        for (int i = 0; i < perm_count; i++) {
          p_perms[i] = malloc(2 * sizeof(unsigned long));
        }
        count = 0;
        permute(p_states, 0, 2, p_perms, &count);

        // phi_nu with its 2p-2h "excitations"
        for (int hi = 0; hi < perm_count; hi++) {
          p3h_steps++;
          for (int pi = 0; pi < perm_count; pi++) {
            p3h_steps++;
            int phi_nu[kmax];
            for (int k = 0; k < kmax; k++) {
              p3h_steps++;
              int replaced = 0;
              for (int idx = 0; idx < 2; idx++) {
                p3h_steps++;
                if (phi[nu][k] == h_perms[hi][idx]) {
                  phi_nu[k] = p_perms[pi][idx];
                  replaced = 1;
                  break;
                }
              }
              if (!replaced)
                phi_nu[k] = phi[nu][k];
            }

            // Track permutations and sorting for sign
            int perm = 0;
            insertionSort_perm_steps(phi_nu, kmax, &perm, &p3h_steps);
            int sign = pow(-1, perm);

            p3h_contr += 1;
            if (sign > 0) {
              p3h_pos += 1;
            } else {
              p3h_neg += 1;
            }
          }
        }
      }

      // 1p1h difference
      if (ph_num == 2) {
        p1h += 1;

        int ph_states[2];
        int ph_statesb[2];
        int ph_states_idx = 0;
        for (int i = 0; i < smax; i++) {
          p1h_steps += 1;
          if ((ph & 1) == 1) {
            ph_states[ph_states_idx] = states[i];
            ph_statesb[ph_states_idx] = statesb[i];
            ph_states_idx++;
          }
          if (ph_states_idx == 2) {
            break;
          }
          ph >>= 1;
        }

        int p_state;
        int h_state;
        for (int i = 0; i < 2; i++) {
          p1h_steps += 1;
          if ((ph_statesb[i] & phib[nu]) == ph_statesb[i]) {
            h_state = ph_states[i];
          } else {
            p_state = ph_states[i];
          }
        }

        // the c_i / c^dagger_j order implemented
        // is c^d_r c_s c^d_a c_b
        // case: b=a, r=p, s=h
        //  so all cases where a h-state can be
        //  annihilated and then a p-state created
        int phi_nu[kmax];
        for (int k = 0; k < kmax; k++) {
          p1h_steps += 1;
          if (phi[nu][k] == h_state) {
            phi_nu[k] = p_state;
          } else {
            phi_nu[k] = phi[nu][k];
          }
        }
        // Track permutations and sorting for sign
        int perm = 0;
        insertionSort_perm_steps(phi_nu, kmax, &perm, &p1h_steps);
        float sign = pow(-1, perm);
        int contr_num = numberOfStates(phib[mu]);
        if (sign > 0) {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_pos++;
          }
        } else {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_neg++;
          }
        }
        // The same happens for
        // case: s=r, a=p, b=h
        //  It has the same sign as before.
        // (The factor will be different in
        // the final implementation)
        // and now the h- / p-summation states are
        // annihilated / created after the ph-states
        contr_num = numberOfStates(phib[mu]);
        if (sign > 0) {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_pos++;
          }
        } else {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_neg++;
          }
        }
        // case: b=r, a=p, s=h
        // Now the summation state is destroyed and
        // the p-state is created in its place, then
        // the h-states is destroyed and the other
        // summation state put in its place.
        // It is the same as before but with the
        // opposite sign.
        contr_num = numberOfStates(phib[nu] & phib[mu]);
        if (sign < 0) {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_pos++;
          }
        } else {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_neg++;
          }
        }
        // case: s=a, r=p, b=h
        // Here the summation states start
        // with a p-creation so the contrubutions
        // are determined by the number of p-states
        // the sign should be opposite to the first
        // two cases
        contr_num = smax - numberOfStates(phib[nu] | phib[mu]);
        if (sign < 0) {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_pos++;
          }
        } else {
          for (int i = 0; i < contr_num; i++) {
            p1h_steps++;
            p1h_contr++;
            p1h_neg++;
          }
        }
      }
      // 0p0h difference
      if (ph_num == 0) {
        p0h += 1;

        // The signs are now all positive
        unsigned long long phi_nu_b = phib[nu];
        unsigned long long phi_nu_b_tmp = phib[nu];
        b = 0;
        while (phi_nu_b_tmp) {
          p0h_steps += 1;
          for (s = 0; s < smax; s++) {
            p0h_steps += 1;

            // b=h-state
            if ((phi_nu_b_tmp & 1) == 1) {
              // s=h-state
              //  b=a & s=r
              if ((phi_nu_b & statesb[s]) == statesb[s]) {
                p0h_contr++;
                p0h_pos++;
                // b <-> s
                // if (states[b] < statesb[s]) {
                //   p0h_contr++;
                //   p0h_pos++;
                // }
              }
              // s=p-state
              //  b=r & s=a
              else {
                p0h_contr++;
                p0h_pos++;
              }
            }
          }
          b++;
          phi_nu_b_tmp >>= 1;
        }
      }
    }
  }
  int cases = p3h + p2h + p1h + p0h;
  int contributions = p3h_contr + p2h_contr + p1h_contr + p0h_contr;
  int steps = p3h_steps + p2h_steps + p1h_steps + p0h_steps;

  printf("### Symmetry Evaluation 2p2h ###\n");
  printf("\n");

  printf("Number of Cases:          %d  =>  100%%\n", cases);
  printf("Number of Contributions:  %d  =>  100%%\n", contributions);
  printf("Number of Steps:          %d  =>  100%%\n", steps);
  printf("\n");

  printf("  Number of 0p0h-Cases:         %d  =>  %f%%\n", p0h,
         (float)p0h / cases * 100);
  printf("  Number of 0p0h-Contributions: %d  =>  %f%%\n", p0h_contr,
         (float)p0h_contr / contributions * 100);
  printf("  Number of 0p0h-Steps:         %d  =>  %f%%\n", p0h_steps,
         (float)p0h_steps / steps * 100);
  printf("    Number of +1 0p0h-Contributions: %d\n", p0h_pos);
  printf("    Number of -1 0p0h-Contributions: %d\n", p0h_neg);
  printf("\n");

  printf("  Number of 1p1h-Cases:         %d  =>  %f%%\n", p1h,
         (float)p1h / cases * 100);
  printf("  Number of 1p1h-Contributions: %d  =>  %f%%\n", p1h_contr,
         (float)p1h_contr / contributions * 100);
  printf("  Number of 1p1h-Steps:         %d  =>  %f%%\n", p1h_steps,
         (float)p1h_steps / steps * 100);
  printf("    Number of +1 1p1h-Contributions: %d\n", p1h_pos);
  printf("    Number of -1 1p1h-Contributions: %d\n", p1h_neg);
  printf("\n");

  printf("  Number of 2p2h-Cases:         %d  =>  %f%%\n", p2h,
         (float)p2h / cases * 100);
  printf("  Number of 2p2h-Contributions: %d  =>  %f%%\n", p2h_contr,
         (float)p2h_contr / contributions * 100);
  printf("  Number of 2p2h-Steps:         %d  =>  %f%%\n", p2h_steps,
         (float)p2h_steps / steps * 100);
  printf("    Number of +1 2p2h-Contributions: %d\n", p2h_pos);
  printf("    Number of -1 2p2h-Contributions: %d\n", p2h_neg);
  printf("\n");

  printf("  Number of 3p3h-Cases:         %d  =>  %f%%\n", p3h,
         (float)p3h / cases * 100);
  printf("  Number of 3p3h-Contributions: %d  =>  %f%%\n", p3h_contr,
         (float)p3h_contr / contributions * 100);
  printf("  Number of 2p2h-Steps:         %d  =>  %f%%\n", p3h_steps,
         (float)p3h_steps / steps * 100);
  printf("    Number of +1 3p3h-Contributions: %d\n", p3h_pos);
  printf("    Number of -1 3p3h-Contributions: %d\n", p3h_neg);
  printf("\n");
}

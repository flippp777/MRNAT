#include "../include/bruteForce_eval.h"
#include "../include/insertionSort_perm.h"
#include "../include/numberOfStates.h"
#include <math.h>
#include <stdio.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void bruteForce_eval_2p2h(int states[], unsigned long long statesb[], int smax,
                          int **phi, int phib[], int kmax, int Mref) {
  int mu; // Index for SD mu
  int nu; // Index for SD nu
  int a;  // Index for particle a
  int r;  // Index for paritcle r
  int b;  // Index for hole b
  int s;  // Index for hole s

  // Values for comparisons
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
      int ph; // Number #p#h excitations between SD's
      ph = numberOfStates(phib[mu] ^ phib[nu]);

      // 2p2h difference
      if (ph == 4) {
        p2h += 1;

        for (b = 0; b < smax; b++) {
          p2h_steps++;
          unsigned long long tmp_b = phib[nu];
          if ((tmp_b & statesb[b]) != statesb[b])
            continue;
          tmp_b ^= statesb[b];

          for (a = 0; a < smax; a++) {
            p2h_steps++;
            unsigned long long tmp_a = tmp_b;
            if ((tmp_a & statesb[a]) != 0)
              continue;
            tmp_a ^= statesb[a];

            for (s = 0; s < smax; s++) {
              p2h_steps++;
              unsigned long long tmp_s = tmp_a;
              if ((tmp_s & statesb[s]) != statesb[s])
                continue;
              tmp_s ^= statesb[s];

              for (r = 0; r < smax; r++) {
                p2h_steps++;
                unsigned long long tmp_r = tmp_s;
                if ((tmp_r & statesb[r]) != 0)
                  continue;
                tmp_r ^= statesb[r];

                if (tmp_r == phib[mu]) {
                  // phi_nu with its 2p-2h "excitations"
                  int phi_nu[kmax];
                  for (int i = 0; i < kmax; i++) {
                    p2h_steps += 1;
                    if (phi[nu][i] == states[b]) {
                      phi_nu[i] = states[a];
                    } else {
                      phi_nu[i] = phi[nu][i];
                    }
                  }
                  for (int i = 0; i < kmax; i++) {
                    p2h_steps += 1;
                    if (phi_nu[i] == states[s]) {
                      phi_nu[i] = states[r];
                      break;
                    }
                  }
                  // Track permutations and sorting for sign
                  int perm = 0;
                  insertionSort_perm_steps(phi_nu, kmax, &perm, &p2h_steps);

                  p2h_contr += 1;
                  if (pow(-1, perm) > 0) {
                    p2h_pos += 1;
                  } else {
                    p2h_neg += 1;
                  }
                }
              } // r
            } // s
          } // a
        } // b
      }
      // 1p1h difference
      if (ph == 2) {
        p1h += 1;

        for (b = 0; b < smax; b++) {
          p1h_steps++;
          unsigned long long tmp_b = phib[nu];
          if ((tmp_b & statesb[b]) != statesb[b])
            continue;
          tmp_b ^= statesb[b];

          for (a = 0; a < smax; a++) {
            p1h_steps++;
            unsigned long long tmp_a = tmp_b;
            if ((tmp_a & statesb[a]) != 0)
              continue;
            tmp_a ^= statesb[a];

            for (s = 0; s < smax; s++) {
              p1h_steps++;
              unsigned long long tmp_s = tmp_a;
              if ((tmp_s & statesb[s]) != statesb[s])
                continue;
              tmp_s ^= statesb[s];

              for (r = 0; r < smax; r++) {
                p1h_steps++;
                unsigned long long tmp_r = tmp_s;
                if ((tmp_r & statesb[r]) != 0)
                  continue;
                tmp_r ^= statesb[r];

                if (tmp_r == phib[mu]) {
                  // phi_nu with its 2p-2h "excitations"
                  int phi_nu[kmax];
                  for (int i = 0; i < kmax; i++) {
                    p1h_steps += 1;
                    if (phi[nu][i] == states[b]) {
                      phi_nu[i] = states[a];
                    } else {
                      phi_nu[i] = phi[nu][i];
                    }
                  }
                  for (int i = 0; i < kmax; i++) {
                    p1h_steps += 1;
                    if (phi_nu[i] == states[s]) {
                      phi_nu[i] = states[r];
                      break;
                    }
                  }
                  // Track permutations and sorting for sign
                  int perm = 0;
                  insertionSort_perm_steps(phi_nu, kmax, &perm, &p1h_steps);

                  p1h_contr += 1;
                  if (pow(-1, perm) > 0) {
                    p1h_pos += 1;
                  } else {
                    p1h_neg += 1;
                  }
                }
              } // r
            } // s
          } // a
        } // b
      }
      // 0p0h difference
      if (ph == 0) {
        p0h += 1;

        for (b = 0; b < smax; b++) {
          p0h_steps++;
          unsigned long long tmp_b = phib[nu];
          if ((tmp_b & statesb[b]) != statesb[b])
            continue;
          tmp_b ^= statesb[b];

          for (a = 0; a < smax; a++) {
            p0h_steps++;
            unsigned long long tmp_a = tmp_b;
            if ((tmp_a & statesb[a]) != 0)
              continue;
            tmp_a ^= statesb[a];

            for (s = 0; s < smax; s++) {
              p0h_steps++;
              unsigned long long tmp_s = tmp_a;
              if ((tmp_s & statesb[s]) != statesb[s])
                continue;
              tmp_s ^= statesb[s];

              for (r = 0; r < smax; r++) {
                p0h_steps++;
                unsigned long long tmp_r = tmp_s;
                if ((tmp_r & statesb[r]) != 0)
                  continue;
                tmp_r ^= statesb[r];

                if (tmp_r == phib[mu]) {
                  // phi_nu with its 2p-2h "excitations"
                  int phi_nu[kmax];
                  for (int i = 0; i < kmax; i++) {
                    p0h_steps += 1;
                    if (phi[nu][i] == states[b]) {
                      phi_nu[i] = states[a];
                    } else {
                      phi_nu[i] = phi[nu][i];
                    }
                  }
                  for (int i = 0; i < kmax; i++) {
                    p0h_steps += 1;
                    if (phi_nu[i] == states[s]) {
                      phi_nu[i] = states[r];
                      break;
                    }
                  }
                  // Track permutations and sorting for sign
                  int perm = 0;
                  insertionSort_perm_steps(phi_nu, kmax, &perm, &p0h_steps);

                  p0h_contr += 1;
                  if (pow(-1, perm) > 0) {
                    p0h_pos += 1;
                  } else {
                    p0h_neg += 1;
                  }
                }
              } // r
            } // s
          } // a
        } // b
      }
    }
  }
  int cases = p2h + p1h + p0h;
  int contributions = p2h_contr + p1h_contr + p0h_contr;
  int steps = p2h_steps + p1h_steps + p0h_steps;

  printf("### Brute Force Evaluation 2p2h ###\n");
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
}

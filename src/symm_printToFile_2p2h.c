#include "../include/count_bits_right.h"
#include "../include/get_npnh.h"
#include "../include/hash.h"
#include "../include/permute.h"
#include "../include/printToFile.h"
#include "../include/printb.h"
#include "../include/quickperm.h"
#include "../include/symm_printToFile.h"
#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void
symm_printToFile_2p2h(FILE* fptr,
                      mp_bitcnt_t smax,
                      int** phi,
                      mpz_t phib[],
                      size_t kmax,
                      size_t Mref)
{
  size_t mu; // Index for SD mu
  size_t nu; // Index for SD nu
  mp_bitcnt_t s;  // Index for summed hole

  mpz_t phi_ph;
  mpz_t phi_nub;
  mpz_t phi_nub_exc;
  mpz_inits(phi_ph, phi_nub, phi_nub_exc, NULL);

  size_t n = 2; // # particles / holes
  // All permutations of particles and holes
  size_t perm_count = 2; // 2!

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
    // Print current SD mu
    fprintf(fptr, "########## phi_mu = ");
    fprintb(fptr, phib[mu], smax);
    fprintf(fptr, " ##########");
    fprintf(fptr, "\n");
    fprintf(fptr, "\n");

    for (nu = 0; nu < Mref; nu++) {
      fprintf(fptr, "########## phi_nu = ");
      fprintb(fptr, phib[nu], smax);
      fprintf(fptr, " ##########");
      fprintf(fptr, "\n");
      fprintf(fptr, "\n");

      mpz_xor(phi_ph, phib[mu], phib[nu]);
      mp_bitcnt_t ph_num = mpz_popcount(phi_ph); // Number #p#h excitations between SD's
      mpz_set(phi_nub, phib[nu]);

      // 2p2h difference
      if (ph_num == 4) {
        p2h++;

        mp_bitcnt_t particles[n];
        mp_bitcnt_t holes[n];
        get_npnh(particles, holes, 4, phi_ph, phi_nub);
        p2h_steps += 4;

        mp_bitcnt_t* h_perms;
        quickperm(holes, n, &h_perms, &perm_count);
        p2h_steps += 2;

        mp_bitcnt_t* p_perms;
        quickperm(particles, n, &p_perms, &perm_count);
        p2h_steps += 2;

        // 4 contributions
        // phi_nu with its 2p-2h "excitations"
        for (size_t hi = 0; hi < perm_count; hi++) {
          p2h_steps++;
          for (size_t pi = 0; pi < perm_count; pi++) {
            p2h_steps++;

            mpz_set(phi_nub_exc, phi_nub);
            int perm = 0;
            for (size_t i = 0; i < n; i++) {
              p2h_steps++;
              mpz_clrbit(phi_nub_exc, h_perms[hi * n + i]);
              perm += count_bits_right(phi_nub_exc, h_perms[hi * n + i]);

              mpz_setbit(phi_nub_exc, p_perms[pi * n + i]);
              perm += count_bits_right(phi_nub_exc, p_perms[pi * n + i]);
            }
            p2h_contr++;
            int sign;
            if (perm % 2 == 0) {
              sign = 1;
              p2h_pos++;
            } else {
              sign = -1;
              p2h_neg++;
            }
            printToFile_2p2h(fptr,
                             phib,
                             smax,
                             phi,
                             kmax,
                             mu,
                             nu,
                             h_perms[hi * n + 0],
                             p_perms[pi * n + 0],
                             h_perms[hi * n + 1],
                             p_perms[pi * n + 1],
                             phi_nub_exc,
                             perm,
                             sign);
          }
        }
        free(h_perms);
        free(p_perms);
      }

      // 1p1h difference
      if (ph_num == 2) {
        p1h++;

        mp_bitcnt_t p_state;
        mp_bitcnt_t h_state;
        get_npnh(&p_state, &h_state, ph_num, phi_ph, phi_nub);
        p1h_steps += 2;

        mp_bitcnt_t* particles[2] = { &s, &p_state };
        mp_bitcnt_t* holes[2] = { &s, &h_state };

        mp_bitcnt_t** h_perms;
        quickperm_ptr(holes, n, &h_perms, &perm_count);
        p1h_steps += 2;

        mp_bitcnt_t** p_perms;
        quickperm_ptr(particles, n, &p_perms, &perm_count);
        p1h_steps += 2;

        for (s = 0; s < smax; s++) {
          p1h_steps += 2;

          size_t set_size = n;
          uint64_t h_keys[set_size];
          uint8_t h_used[set_size];
          HashSet H;
          hashset_init(&H, h_keys, h_used, set_size);

          for (size_t hi = 0; hi < perm_count; hi++) {
            p1h_steps += 2;

            uint64_t h_hash = hash_func_ptr(h_perms, hi, n);
            if (hashset_check(&H, h_hash))
              continue;

            uint64_t p_keys[set_size];
            uint8_t p_used[set_size];
            HashSet P;
            hashset_init(&P, p_keys, p_used, set_size);

            for (size_t pi = 0; pi < perm_count; pi++) {
              p1h_steps += 2;

              uint64_t p_hash = hash_func_ptr(p_perms, pi, n);
              if (hashset_check(&P, p_hash))
                continue;

              mpz_set(phi_nub_exc, phi_nub);
              int perm = 0;
              int valid = 1;
              for (size_t i = 0; i < n; i++) {
                if (mpz_tstbit(phi_nub_exc, *h_perms[hi * n + i]) != 1) {
                  valid = 0;
                  break;
                }
                mpz_clrbit(phi_nub_exc, *h_perms[hi * n + i]);
                perm += count_bits_right(phi_nub_exc, *h_perms[hi * n + i]);

                if (mpz_tstbit(phi_nub_exc, *p_perms[pi * n + i]) != 0) {
                  valid = 0;
                  break;
                }
                mpz_setbit(phi_nub_exc, *p_perms[pi * n + i]);
                perm += count_bits_right(phi_nub_exc, *p_perms[pi * n + i]);
              }
              if (valid) {
                p1h_contr++;
                int sign;
                if (perm % 2 == 0) {
                  p1h_pos++;
                  sign = 1;
                } else {
                  p1h_neg++;
                  sign = -1;
                }
                printToFile_2p2h(fptr,
                                 phib,
                                 smax,
                                 phi,
                                 kmax,
                                 mu,
                                 nu,
                                 *h_perms[hi * n + 0],
                                 *p_perms[pi * n + 0],
                                 *h_perms[hi * n + 1],
                                 *p_perms[pi * n + 1],
                                 phi_nub_exc,
                                 perm,
                                 sign);
              }
            }
          }
        }

        free(h_perms);
        free(p_perms);
      }

      // 0p0h difference
      if (ph_num == 0) {
        p0h++;

        // The signs are now all positive
        int b = 0;
        for (size_t k = 0; k < kmax; k++) {
          p0h_steps++;
          // b=h-state
          b = mpz_scan1(phi_nub, 0);
          for (s = 0; s < smax; s++) {
            p0h_steps++;
            // s=h-state
            //  b=a & s=r
            if (mpz_tstbit(phi_nub, s) == 1) {
              p0h_contr++;
              p0h_pos++;
              mpz_clrbit(phi_nub_exc, b);
              mpz_setbit(phi_nub_exc, b);
              mpz_clrbit(phi_nub_exc, s);
              mpz_clrbit(phi_nub_exc, s);

              printToFile_2p2h(fptr,
                               phib,
                               smax,
                               phi,
                               kmax,
                               mu,
                               nu,
                               b,
                               b,
                               s,
                               s,
                               phi_nub_exc,
                               0,
                               1);
            }
            // s=p-state
            //  b=r & s=a
            else {
              p0h_contr++;
              p0h_pos++;
              mpz_clrbit(phi_nub_exc, b);
              mpz_setbit(phi_nub_exc, s);
              mpz_clrbit(phi_nub_exc, s);
              mpz_clrbit(phi_nub_exc, b);

              printToFile_2p2h(fptr,
                               phib,
                               smax,
                               phi,
                               kmax,
                               mu,
                               nu,
                               b,
                               s,
                               s,
                               b,
                               phi_nub_exc,
                               0,
                               1);
            }
          }
          b++;
        }
      }
    }
  }
  mpz_clears(phi_ph, phi_nub, phi_nub_exc, NULL);

  int cases = p2h + p1h + p0h;
  int contributions = p2h_contr + p1h_contr + p0h_contr;
  int steps = p2h_steps + p1h_steps + p0h_steps;

  printf("### Symmetry Force Evaluation 2p2h ###\n");
  printf("\n");

  printf("Number of Cases:          %d  =>  100%%\n", cases);
  printf("Number of Contributions:  %d  =>  100%%\n", contributions);
  printf("Number of Steps:          %d  =>  100%%\n", steps);
  printf("\n");

  printf("  Number of 0p0h-Cases:         %d  =>  %f%%\n",
         p0h,
         (float)p0h / cases * 100);
  printf("  Number of 0p0h-Contributions: %d  =>  %f%%\n",
         p0h_contr,
         (float)p0h_contr / contributions * 100);
  printf("  Number of 0p0h-Steps:         %d  =>  %f%%\n",
         p0h_steps,
         (float)p0h_steps / steps * 100);
  printf("    Number of +1 0p0h-Contributions: %d\n", p0h_pos);
  printf("    Number of -1 0p0h-Contributions: %d\n", p0h_neg);
  printf("\n");

  printf("  Number of 1p1h-Cases:         %d  =>  %f%%\n",
         p1h,
         (float)p1h / cases * 100);
  printf("  Number of 1p1h-Contributions: %d  =>  %f%%\n",
         p1h_contr,
         (float)p1h_contr / contributions * 100);
  printf("  Number of 1p1h-Steps:         %d  =>  %f%%\n",
         p1h_steps,
         (float)p1h_steps / steps * 100);
  printf("    Number of +1 1p1h-Contributions: %d\n", p1h_pos);
  printf("    Number of -1 1p1h-Contributions: %d\n", p1h_neg);
  printf("\n");

  printf("  Number of 2p2h-Cases:         %d  =>  %f%%\n",
         p2h,
         (float)p2h / cases * 100);
  printf("  Number of 2p2h-Contributions: %d  =>  %f%%\n",
         p2h_contr,
         (float)p2h_contr / contributions * 100);
  printf("  Number of 2p2h-Steps:         %d  =>  %f%%\n",
         p2h_steps,
         (float)p2h_steps / steps * 100);
  printf("    Number of +1 2p2h-Contributions: %d\n", p2h_pos);
  printf("    Number of -1 2p2h-Contributions: %d\n", p2h_neg);
  printf("\n");
}

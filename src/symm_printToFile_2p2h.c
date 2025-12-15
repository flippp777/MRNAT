#include "../include/count_bits_right.h"
#include "../include/get_npnh.h"
#include "../include/permute.h"
#include "../include/printToFile.h"
#include "../include/printb.h"
#include "../include/symm_printToFile.h"
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void
symm_printToFile_2p2h(FILE* fptr,
                      int states[],
                      int smax,
                      int** phi,
                      mpz_t phib[],
                      int kmax,
                      int Mref)
{
  int mu; // Index for SD mu
  int nu; // Index for SD nu
  int a;  // Index for particle a
  int r;  // Index for paritcle r
  int b;  // Index for hole b
  int s;  // Index for hole s
  mpz_t phi_ph;
  mpz_t phi_nub;
  mpz_t phi_nub_exc;
  mpz_inits(phi_ph, phi_nub, phi_nub_exc, NULL);

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
      int ph_num = mpz_popcount(phi_ph); // Number #p#h excitations between SD's
      int n = ph_num / 2;
      mpz_set(phi_nub, phib[nu]);

      // 2p2h difference
      if (ph_num == 4) {

        int particles[2];
        int holes[2];
        get_npnh(particles, holes, 4, smax, phi_ph, phi_nub);

        // All permutations of particles and holes
        int perm_count = 2; // 2!

        int** h_perms;
        permute(h_perms, holes, 0, n, perm_count);

        int** p_perms;
        permute(p_perms, particles, 0, n, perm_count);

        // 4 contributions
        // phi_nu with its 2p-2h "excitations"
        for (int hi = 0; hi < perm_count; hi++) {
          mpz_set(phi_nub_exc, phi_nub);
          for (int pi = 0; pi < perm_count; pi++) {
            int perm = 0;
            for (int i = 0; i < n; i++) {
              mpz_clrbit(phi_nub_exc, h_perms[hi][i]);
              perm += count_bits_right(phi_nub_exc, h_perms[hi][i]);

              mpz_setbit(phi_nub_exc, p_perms[pi][i]);
              perm += count_bits_right(phi_nub_exc, p_perms[pi][i]);
            }
            int sign;
            if (perm % 2 == 0) {
              sign = 1;
            } else {
              sign = -1;
            }
            printToFile(fptr,
                        phib,
                        smax,
                        phi,
                        kmax,
                        mu,
                        nu,
                        h_perms[hi][1],
                        p_perms[pi][1],
                        h_perms[hi][2],
                        p_perms[pi][2],
                        phi_nub_exc,
                        sign);
          }
        }

        // 1p1h difference
        if (ph_num == 2) {

          int p_state;
          int h_state;
          get_npnh(&p_state, &h_state, ph_num, smax, phi_ph, phi_nub);

          int s;
          int* particles[2] = { &s, &p_state };
          int* holes[2] = { &s, &h_state };

          int perm_count = 2; // 2!

          int*** h_perms;
          permute_ptr(h_perms, holes, 0, n, perm_count);

          int*** p_perms;
          permute_ptr(p_perms, particles, 0, n, perm_count);

          for (s = 0; s < smax; s++) {
            for (int hi = 0; hi < perm_count; hi++) {
              mpz_set(phi_nub_exc, phi_nub);

              for (int pi = 0; pi < perm_count; pi++) {
                int perm = 0;
                int valid = 1;
                for (int i = 0; i < n; i++) {
                  if (mpz_tstbit(phi_nub, *h_perms[hi][i]) != 1) {
                    valid = 0;
                    break;
                  }
                  mpz_clrbit(phi_nub_exc, *h_perms[hi][i]);
                  perm += count_bits_right(phi_nub_exc, *h_perms[hi][i]);

                  if (mpz_tstbit(phi_nub_exc, *p_perms[pi][i]) != 0) {
                    valid = 0;
                    break;
                  }
                  mpz_setbit(phi_nub_exc, *p_perms[pi][i]);
                  perm += count_bits_right(phi_nub_exc, *p_perms[pi][i]);
                }
                if (valid) {
                  int sign;
                  if (perm % 2 == 0) {
                    sign = 1;
                  } else {
                    sign = -1;
                  }
                  printToFile(fptr,
                              phib,
                              smax,
                              phi,
                              kmax,
                              mu,
                              nu,
                              *h_perms[hi][1],
                              *p_perms[pi][1],
                              *h_perms[hi][2],
                              *p_perms[pi][2],
                              phi_nub_exc,
                              sign);
                }
              }
            }
          }
        }

        // 0p0h difference
        if (ph_num == 0) {

          // The signs are now all positive
          int b = 0;
          for (int k = 0; k < kmax; k++) {
            // b=h-state
            b = mpz_scan1(phi_nub, 0);
            for (s = 0; s < smax; s++) {
              // s=h-state
              //  b=a & s=r
              if (mpz_tstbit(phi_nub, s) == 1) {
                mpz_clrbit(phi_nub_exc, b);
                mpz_setbit(phi_nub_exc, b);
                mpz_clrbit(phi_nub_exc, s);
                mpz_clrbit(phi_nub_exc, s);

                printToFile(fptr,
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
                            1);
              }
              // s=p-state
              //  b=r & s=a
              else {
                mpz_clrbit(phi_nub_exc, b);
                mpz_setbit(phi_nub_exc, s);
                mpz_clrbit(phi_nub_exc, s);
                mpz_clrbit(phi_nub_exc, b);

                printToFile(fptr,
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
                            1);
              }
            }
            b++;
          }
        }
      }
    }
  }
  mpz_clears(phi_ph, phi_nub, phi_nub_exc, NULL);
}

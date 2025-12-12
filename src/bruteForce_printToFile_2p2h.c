#include "../include/bruteForce_printToFile.h"
#include "../include/printToFile.h"
#include "../include/printb.h"
#include <gmp.h>
#include <stdio.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void
bruteForce_printToFile_2p2h(FILE* fptr,
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

  mpz_t tmp_b, tmp_a, tmp_s, tmp_r;
  mpz_inits(tmp_b, tmp_a, tmp_s, tmp_r, NULL);

  for (mu = 0; mu < Mref; mu++) {
    // Print current SD mu
    fprintf(fptr, "########## phi_mu = ");
    fprintb(fptr, phib[mu], smax);
    fprintf(fptr, " ##########");
    fprintf(fptr, "\n");

    for (nu = 0; nu < Mref; nu++) {
      fprintf(fptr, "########## phi_nu = ");
      fprintb(fptr, phib[nu], smax);
      fprintf(fptr, " ##########");
      fprintf(fptr, "\n");
      fprintf(fptr, "\n");

      for (b = 0; b < smax; b++) {
        mpz_set(tmp_b, phib[nu]);
        if (mpz_tstbit(tmp_b, b) == 0)
          continue;
        mpz_combit(tmp_b, b);

        for (a = 0; a < smax; a++) {
          mpz_set(tmp_a, tmp_b);
          if (mpz_tstbit(tmp_a, a) == 1)
            continue;
          mpz_combit(tmp_a, a);

          for (s = 0; s < smax; s++) {
            mpz_set(tmp_s, tmp_a);
            if (mpz_tstbit(tmp_s, s) == 0)
              continue;
            mpz_combit(tmp_s, s);

            for (r = 0; r < smax; r++) {
              mpz_set(tmp_r, tmp_s);
              if (mpz_tstbit(tmp_r, r) == 1)
                continue;
              mpz_combit(tmp_r, r);

              if (mpz_cmp(tmp_r, phib[mu]) == 0) {

                printToFile(fptr,
                            phib,
                            smax,
                            phi,
                            kmax,
                            mu,
                            nu,
                            b,
                            a,
                            s,
                            r,
                            tmp_b,
                            tmp_a,
                            tmp_s,
                            tmp_r);
              }
            }
          }
        }
      }
    }
  }
  mpz_clears(tmp_b, tmp_a, tmp_s, tmp_r, NULL);
}

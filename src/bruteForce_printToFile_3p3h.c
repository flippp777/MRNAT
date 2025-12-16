#include "../include/bruteForce_printToFile.h"
#include "../include/count_bits_right.h"
#include "../include/printToFile.h"
#include "../include/printb.h"
#include <gmp.h>
#include <stdio.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void
bruteForce_printToFile_3p3h(FILE* fptr,
                            int smax,
                            int** phi,
                            mpz_t phib[],
                            int kmax,
                            int Mref)
{
  int mu; // Index for SD mu
  int nu; // Index for SD nu
  int a;  // Index for particle a
  int c;  // Index for particle c
  int r;  // Index for paritcle r
  int b;  // Index for hole b
  int d;  // Index for hole d
  int s;  // Index for hole s

  mpz_t phi_ph, tmp_b, tmp_a, tmp_d, tmp_c, tmp_s, tmp_r;
  mpz_inits(phi_ph, tmp_b, tmp_a, tmp_d, tmp_c, tmp_s, tmp_r, NULL);

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

      mpz_xor(phi_ph, phib[mu], phib[nu]);
      int ph_num = mpz_popcount(phi_ph); // Number #p#h excitations between SD's

      if (ph_num == 6) {
        p3h++;
        for (b = 0; b < smax; b++) {
          p3h_steps++;
          mpz_set(tmp_b, phib[nu]);
          if (mpz_tstbit(tmp_b, b) == 0)
            continue;
          mpz_combit(tmp_b, b);

          for (a = 0; a < smax; a++) {
            p3h_steps++;
            mpz_set(tmp_a, tmp_b);
            if (mpz_tstbit(tmp_a, a) == 1)
              continue;
            mpz_combit(tmp_a, a);

            for (d = 0; d < smax; d++) {
              p3h_steps++;
              mpz_set(tmp_d, tmp_a);
              if (mpz_tstbit(tmp_d, d) == 0)
                continue;
              mpz_combit(tmp_d, d);

              for (c = 0; c < smax; c++) {
                p3h_steps++;
                mpz_set(tmp_c, tmp_d);
                if (mpz_tstbit(tmp_c, c) == 1)
                  continue;
                mpz_combit(tmp_c, c);

                for (s = 0; s < smax; s++) {
                  p3h_steps++;
                  mpz_set(tmp_s, tmp_c);
                  if (mpz_tstbit(tmp_s, s) == 0)
                    continue;
                  mpz_combit(tmp_s, s);

                  for (r = 0; r < smax; r++) {
                    p3h_steps++;
                    mpz_set(tmp_r, tmp_s);
                    if (mpz_tstbit(tmp_r, r) == 1)
                      continue;
                    mpz_combit(tmp_r, r);

                    if (mpz_cmp(tmp_r, phib[mu]) == 0) {
                      p3h_contr++;

                      int perm = count_bits_right(tmp_b, b) +
                                 count_bits_right(tmp_a, a) +
                                 count_bits_right(tmp_d, d) +
                                 count_bits_right(tmp_c, c) +
                                 count_bits_right(tmp_s, s) +
                                 count_bits_right(tmp_r, r);
                      int sign;
                      if (perm % 2 == 0) {
                        p3h_pos++;
                        sign = 1;
                      } else {
                        p3h_neg++;
                        sign = -1;
                      }

                      printToFile_3p3h(fptr,
                                       phib,
                                       smax,
                                       phi,
                                       kmax,
                                       mu,
                                       nu,
                                       b,
                                       a,
                                       d,
                                       c,
                                       s,
                                       r,
                                       tmp_r,
                                       perm,
                                       sign);
                    } // if
                  } // r
                } // s
              } // c
            } // d
          } // a
        } // b
      } // if ph_num
      if (ph_num == 4) {
        p2h++;
        for (b = 0; b < smax; b++) {
          p2h_steps++;
          mpz_set(tmp_b, phib[nu]);
          if (mpz_tstbit(tmp_b, b) == 0)
            continue;
          mpz_combit(tmp_b, b);

          for (a = 0; a < smax; a++) {
            p2h_steps++;
            mpz_set(tmp_a, tmp_b);
            if (mpz_tstbit(tmp_a, a) == 1)
              continue;
            mpz_combit(tmp_a, a);

            for (d = 0; d < smax; d++) {
              p3h_steps++;
              mpz_set(tmp_d, tmp_a);
              if (mpz_tstbit(tmp_d, d) == 0)
                continue;
              mpz_combit(tmp_d, d);

              for (c = 0; c < smax; c++) {
                p3h_steps++;
                mpz_set(tmp_c, tmp_d);
                if (mpz_tstbit(tmp_c, c) == 1)
                  continue;
                mpz_combit(tmp_c, c);

                for (s = 0; s < smax; s++) {
                  p2h_steps++;
                  mpz_set(tmp_s, tmp_c);
                  if (mpz_tstbit(tmp_s, s) == 0)
                    continue;
                  mpz_combit(tmp_s, s);

                  for (r = 0; r < smax; r++) {
                    p2h_steps++;
                    mpz_set(tmp_r, tmp_s);
                    if (mpz_tstbit(tmp_r, r) == 1)
                      continue;
                    mpz_combit(tmp_r, r);

                    if (mpz_cmp(tmp_r, phib[mu]) == 0) {
                      p2h_contr++;

                      int perm = count_bits_right(tmp_b, b) +
                                 count_bits_right(tmp_a, a) +
                                 count_bits_right(tmp_d, d) +
                                 count_bits_right(tmp_c, c) +
                                 count_bits_right(tmp_s, s) +
                                 count_bits_right(tmp_r, r);
                      int sign;
                      if (perm % 2 == 0) {
                        p2h_pos++;
                        sign = 1;
                      } else {
                        p2h_neg++;
                        sign = -1;
                      }

                      printToFile_3p3h(fptr,
                                       phib,
                                       smax,
                                       phi,
                                       kmax,
                                       mu,
                                       nu,
                                       b,
                                       a,
                                       d,
                                       c,
                                       s,
                                       r,
                                       tmp_r,
                                       perm,
                                       sign);
                    } // if
                  } // r
                } // s
              } // c
            } // d
          } // a
        } // b
      } // if ph_num
      if (ph_num == 2) {
        p1h++;
        for (b = 0; b < smax; b++) {
          p1h_steps++;
          mpz_set(tmp_b, phib[nu]);
          if (mpz_tstbit(tmp_b, b) == 0)
            continue;
          mpz_combit(tmp_b, b);

          for (a = 0; a < smax; a++) {
            p1h_steps++;
            mpz_set(tmp_a, tmp_b);
            if (mpz_tstbit(tmp_a, a) == 1)
              continue;
            mpz_combit(tmp_a, a);

            for (d = 0; d < smax; d++) {
              p3h_steps++;
              mpz_set(tmp_d, tmp_a);
              if (mpz_tstbit(tmp_d, d) == 0)
                continue;
              mpz_combit(tmp_d, d);

              for (c = 0; c < smax; c++) {
                p3h_steps++;
                mpz_set(tmp_c, tmp_d);
                if (mpz_tstbit(tmp_c, c) == 1)
                  continue;
                mpz_combit(tmp_c, c);

                for (s = 0; s < smax; s++) {
                  p1h_steps++;
                  mpz_set(tmp_s, tmp_c);
                  if (mpz_tstbit(tmp_s, s) == 0)
                    continue;
                  mpz_combit(tmp_s, s);

                  for (r = 0; r < smax; r++) {
                    p1h_steps++;
                    mpz_set(tmp_r, tmp_s);
                    if (mpz_tstbit(tmp_r, r) == 1)
                      continue;
                    mpz_combit(tmp_r, r);

                    if (mpz_cmp(tmp_r, phib[mu]) == 0) {
                      p1h_contr++;

                      int perm = count_bits_right(tmp_b, b) +
                                 count_bits_right(tmp_a, a) +
                                 count_bits_right(tmp_d, d) +
                                 count_bits_right(tmp_c, c) +
                                 count_bits_right(tmp_s, s) +
                                 count_bits_right(tmp_r, r);
                      int sign;
                      if (perm % 2 == 0) {
                        p1h_pos++;
                        sign = 1;
                      } else {
                        p1h_neg++;
                        sign = -1;
                      }

                      printToFile_3p3h(fptr,
                                       phib,
                                       smax,
                                       phi,
                                       kmax,
                                       mu,
                                       nu,
                                       b,
                                       a,
                                       d,
                                       c,
                                       s,
                                       r,
                                       tmp_r,
                                       perm,
                                       sign);
                    } // if
                  } // r
                } // s
              } // c
            } // d
          } // a
        } // b
      } // if ph_num
      if (ph_num == 0) {
        p0h++;
        for (b = 0; b < smax; b++) {
          p0h_steps++;
          mpz_set(tmp_b, phib[nu]);
          if (mpz_tstbit(tmp_b, b) == 0)
            continue;
          mpz_combit(tmp_b, b);

          for (a = 0; a < smax; a++) {
            p0h_steps++;
            mpz_set(tmp_a, tmp_b);
            if (mpz_tstbit(tmp_a, a) == 1)
              continue;
            mpz_combit(tmp_a, a);

            for (d = 0; d < smax; d++) {
              p3h_steps++;
              mpz_set(tmp_d, tmp_a);
              if (mpz_tstbit(tmp_d, d) == 0)
                continue;
              mpz_combit(tmp_d, d);

              for (c = 0; c < smax; c++) {
                p3h_steps++;
                mpz_set(tmp_c, tmp_d);
                if (mpz_tstbit(tmp_c, c) == 1)
                  continue;
                mpz_combit(tmp_c, c);

                for (s = 0; s < smax; s++) {
                  p0h_steps++;
                  mpz_set(tmp_s, tmp_c);
                  if (mpz_tstbit(tmp_s, s) == 0)
                    continue;
                  mpz_combit(tmp_s, s);

                  for (r = 0; r < smax; r++) {
                    p0h_steps++;
                    mpz_set(tmp_r, tmp_s);
                    if (mpz_tstbit(tmp_r, r) == 1)
                      continue;
                    mpz_combit(tmp_r, r);

                    if (mpz_cmp(tmp_r, phib[mu]) == 0) {
                      p0h_contr++;

                      int perm = count_bits_right(tmp_b, b) +
                                 count_bits_right(tmp_a, a) +
                                 count_bits_right(tmp_d, d) +
                                 count_bits_right(tmp_c, c) +
                                 count_bits_right(tmp_s, s) +
                                 count_bits_right(tmp_r, r);
                      int sign;
                      if (perm % 2 == 0) {
                        p0h_pos++;
                        sign = 1;
                      } else {
                        p0h_neg++;
                        sign = -1;
                      }

                      printToFile_3p3h(fptr,
                                       phib,
                                       smax,
                                       phi,
                                       kmax,
                                       mu,
                                       nu,
                                       b,
                                       a,
                                       d,
                                       c,
                                       s,
                                       r,
                                       tmp_r,
                                       perm,
                                       sign);
                    } // if
                  } // r
                } // s
              } // c
            } // d
          } // a
        } // b
      }
    }
  }
  mpz_clears(tmp_b, tmp_a, tmp_d, tmp_c, tmp_s, tmp_r, NULL);

  int cases = p3h + p2h + p1h + p0h;
  int contributions = p3h_contr + p2h_contr + p1h_contr + p0h_contr;
  int steps = p3h_steps + p2h_steps + p1h_steps + p0h_steps;

  printf("### Brute Force Evaluation 3p3h ###\n");
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

  printf("  Number of 3p3h-Cases:         %d  =>  %f%%\n",
         p3h,
         (float)p3h / cases * 100);
  printf("  Number of 3p3h-Contributions: %d  =>  %f%%\n",
         p3h_contr,
         (float)p3h_contr / contributions * 100);
  printf("  Number of 3p3h-Steps:         %d  =>  %f%%\n",
         p3h_steps,
         (float)p3h_steps / steps * 100);
  printf("    Number of +1 3p3h-Contributions: %d\n", p3h_pos);
  printf("    Number of -1 3p3h-Contributions: %d\n", p3h_neg);
  printf("\n");
}

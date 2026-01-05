#include "../include/bruteForce_eval.h"
#include "../include/count_bits_right.h"
#include "../include/printToFile.h"
#include "../include/printb.h"
#include <gmp.h>
#include <stdio.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void
bruteForce_eval_5p5h(int smax, mpz_t phib[], int kmax, int Mref)
{
  int mu; // Index for SD mu
  int nu; // Index for SD nu
  int p1; // Index for particle
  int h1; // Index for hole
  int p2; // Index for particle
  int h2; // Index for hole
  int p3; // Index for paritcle
  int h3; // Index for hole
  int p4; // Index for paritcle
  int h4; // Index for hole
  int p5; // Index for paritcle
  int h5; // Index for hole

  mpz_t phi_ph, tmp_h1, tmp_p1, tmp_h2, tmp_p2, tmp_h3, tmp_p3, tmp_h4, tmp_p4,
    tmp_h5, tmp_p5;
  mpz_inits(phi_ph,
            tmp_h1,
            tmp_p1,
            tmp_h2,
            tmp_p2,
            tmp_h3,
            tmp_p3,
            tmp_h4,
            tmp_p4,
            tmp_h5,
            tmp_p5,
            NULL);

  // Values for comparisons
  int p5h = 0; // Number of 5p5h case
  int p5h_contr = 0;
  int p5h_pos = 0;
  int p5h_neg = 0;
  unsigned long p5h_steps = 0;
  int p4h = 0; // Number of 4p4h case
  int p4h_contr = 0;
  int p4h_pos = 0;
  int p4h_neg = 0;
  unsigned long p4h_steps = 0;
  int p3h = 0; // Number of 3p3h case
  int p3h_contr = 0;
  int p3h_pos = 0;
  int p3h_neg = 0;
  unsigned long p3h_steps = 0;
  int p2h = 0; // Number of 2p2h case
  int p2h_contr = 0;
  int p2h_pos = 0;
  int p2h_neg = 0;
  unsigned long p2h_steps = 0;
  int p1h = 0; // Number of 1p1h case
  int p1h_contr = 0;
  int p1h_pos = 0;
  int p1h_neg = 0;
  int p1h_steps = 0;
  int p0h = 0; // Number of 0p0h case
  int p0h_contr = 0;
  int p0h_pos = 0;
  int p0h_neg = 0;
  unsigned long p0h_steps = 0;

  for (mu = 0; mu < Mref; mu++) {

    for (nu = 0; nu < Mref; nu++) {

      mpz_xor(phi_ph, phib[mu], phib[nu]);
      int ph_num = mpz_popcount(phi_ph); // Number #p#h excitations between SD's

      if (ph_num == 10) {
        p5h++;
        for (h1 = 0; h1 < smax; h1++) {
          p5h_steps++;
          mpz_set(tmp_h1, phib[nu]);
          if (mpz_tstbit(tmp_h1, h1) == 0)
            continue;
          mpz_combit(tmp_h1, h1);

          for (p1 = 0; p1 < smax; p1++) {
            p5h_steps++;
            mpz_set(tmp_p1, tmp_h1);
            if (mpz_tstbit(tmp_p1, p1) == 1)
              continue;
            mpz_combit(tmp_p1, p1);

            for (h2 = 0; h2 < smax; h2++) {
              p5h_steps++;
              mpz_set(tmp_h2, tmp_p1);
              if (mpz_tstbit(tmp_h2, h2) == 0)
                continue;
              mpz_combit(tmp_h2, h2);

              for (p2 = 0; p2 < smax; p2++) {
                p5h_steps++;
                mpz_set(tmp_p2, tmp_h2);
                if (mpz_tstbit(tmp_p2, p2) == 1)
                  continue;
                mpz_combit(tmp_p2, p2);

                for (h3 = 0; h3 < smax; h3++) {
                  p5h_steps++;
                  mpz_set(tmp_h3, tmp_p2);
                  if (mpz_tstbit(tmp_h3, h3) == 0)
                    continue;
                  mpz_combit(tmp_h3, h3);

                  for (p3 = 0; p3 < smax; p3++) {
                    p5h_steps++;
                    mpz_set(tmp_p3, tmp_h3);
                    if (mpz_tstbit(tmp_p3, p3) == 1)
                      continue;
                    mpz_combit(tmp_p3, p3);

                    for (h4 = 0; h4 < smax; h4++) {
                      p5h_steps++;
                      mpz_set(tmp_h4, tmp_p3);
                      if (mpz_tstbit(tmp_h4, h4) == 0)
                        continue;
                      mpz_combit(tmp_h4, h4);

                      for (p4 = 0; p4 < smax; p4++) {
                        p5h_steps++;
                        mpz_set(tmp_p4, tmp_h4);
                        if (mpz_tstbit(tmp_p4, p4) == 1)
                          continue;
                        mpz_combit(tmp_p4, p4);

                        for (h5 = 0; h5 < smax; h5++) {
                          p5h_steps++;
                          mpz_set(tmp_h5, tmp_p4);
                          if (mpz_tstbit(tmp_h5, h5) == 0)
                            continue;
                          mpz_combit(tmp_h5, h5);

                          for (p5 = 0; p5 < smax; p5++) {
                            p5h_steps++;
                            mpz_set(tmp_p5, tmp_h5);
                            if (mpz_tstbit(tmp_p5, p5) == 1)
                              continue;
                            mpz_combit(tmp_p5, p5);

                            if (mpz_cmp(tmp_p5, phib[mu]) == 0) {
                              p5h_contr++;

                              int perm = count_bits_right(tmp_h1, h1) +
                                         count_bits_right(tmp_p1, p1) +
                                         count_bits_right(tmp_h2, h2) +
                                         count_bits_right(tmp_p2, p2) +
                                         count_bits_right(tmp_h3, h3) +
                                         count_bits_right(tmp_p3, p3) +
                                         count_bits_right(tmp_h4, h4) +
                                         count_bits_right(tmp_p4, p4) +
                                         count_bits_right(tmp_h5, h5) +
                                         count_bits_right(tmp_p5, p5);
                              if (perm % 2 == 0) {
                                p5h_pos++;
                              } else {
                                p5h_neg++;
                              }
                            } // if
                          } // r
                        }
                      }
                    }
                  }
                } // s
              } // c
            } // d
          } // a
        } // b
      } // if ph_num
      if (ph_num == 8) {
        p4h++;
        for (h1 = 0; h1 < smax; h1++) {
          p4h_steps++;
          mpz_set(tmp_h1, phib[nu]);
          if (mpz_tstbit(tmp_h1, h1) == 0)
            continue;
          mpz_combit(tmp_h1, h1);

          for (p1 = 0; p1 < smax; p1++) {
            p4h_steps++;
            mpz_set(tmp_p1, tmp_h1);
            if (mpz_tstbit(tmp_p1, p1) == 1)
              continue;
            mpz_combit(tmp_p1, p1);

            for (h2 = 0; h2 < smax; h2++) {
              p4h_steps++;
              mpz_set(tmp_h2, tmp_p1);
              if (mpz_tstbit(tmp_h2, h2) == 0)
                continue;
              mpz_combit(tmp_h2, h2);

              for (p2 = 0; p2 < smax; p2++) {
                p4h_steps++;
                mpz_set(tmp_p2, tmp_h2);
                if (mpz_tstbit(tmp_p2, p2) == 1)
                  continue;
                mpz_combit(tmp_p2, p2);

                for (h3 = 0; h3 < smax; h3++) {
                  p4h_steps++;
                  mpz_set(tmp_h3, tmp_p2);
                  if (mpz_tstbit(tmp_h3, h3) == 0)
                    continue;
                  mpz_combit(tmp_h3, h3);

                  for (p3 = 0; p3 < smax; p3++) {
                    p4h_steps++;
                    mpz_set(tmp_p3, tmp_h3);
                    if (mpz_tstbit(tmp_p3, p3) == 1)
                      continue;
                    mpz_combit(tmp_p3, p3);

                    for (h4 = 0; h4 < smax; h4++) {
                      p4h_steps++;
                      mpz_set(tmp_h4, tmp_p3);
                      if (mpz_tstbit(tmp_h4, h4) == 0)
                        continue;
                      mpz_combit(tmp_h4, h4);

                      for (p4 = 0; p4 < smax; p4++) {
                        p4h_steps++;
                        mpz_set(tmp_p4, tmp_h4);
                        if (mpz_tstbit(tmp_p4, p4) == 1)
                          continue;
                        mpz_combit(tmp_p4, p4);

                        for (h5 = 0; h5 < smax; h5++) {
                          p4h_steps++;
                          mpz_set(tmp_h5, tmp_p4);
                          if (mpz_tstbit(tmp_h5, h5) == 0)
                            continue;
                          mpz_combit(tmp_h5, h5);

                          for (p5 = 0; p5 < smax; p5++) {
                            p4h_steps++;
                            mpz_set(tmp_p5, tmp_h5);
                            if (mpz_tstbit(tmp_p5, p5) == 1)
                              continue;
                            mpz_combit(tmp_p5, p5);

                            if (mpz_cmp(tmp_p5, phib[mu]) == 0) {
                              p4h_contr++;

                              int perm = count_bits_right(tmp_h1, h1) +
                                         count_bits_right(tmp_p1, p1) +
                                         count_bits_right(tmp_h2, h2) +
                                         count_bits_right(tmp_p2, p2) +
                                         count_bits_right(tmp_h3, h3) +
                                         count_bits_right(tmp_p3, p3) +
                                         count_bits_right(tmp_h4, h4) +
                                         count_bits_right(tmp_p4, p4) +
                                         count_bits_right(tmp_h5, h5) +
                                         count_bits_right(tmp_p5, p5);
                              if (perm % 2 == 0) {
                                p4h_pos++;
                              } else {
                                p4h_neg++;
                              }
                            } // if
                          } // r
                        }
                      }
                    }
                  }
                } // s
              } // c
            } // d
          } // a
        } // b
      } // if ph_num
      if (ph_num == 6) {
        p3h++;
        for (h1 = 0; h1 < smax; h1++) {
          p3h_steps++;
          mpz_set(tmp_h1, phib[nu]);
          if (mpz_tstbit(tmp_h1, h1) == 0)
            continue;
          mpz_combit(tmp_h1, h1);

          for (p1 = 0; p1 < smax; p1++) {
            p3h_steps++;
            mpz_set(tmp_p1, tmp_h1);
            if (mpz_tstbit(tmp_p1, p1) == 1)
              continue;
            mpz_combit(tmp_p1, p1);

            for (h2 = 0; h2 < smax; h2++) {
              p3h_steps++;
              mpz_set(tmp_h2, tmp_p1);
              if (mpz_tstbit(tmp_h2, h2) == 0)
                continue;
              mpz_combit(tmp_h2, h2);

              for (p2 = 0; p2 < smax; p2++) {
                p3h_steps++;
                mpz_set(tmp_p2, tmp_h2);
                if (mpz_tstbit(tmp_p2, p2) == 1)
                  continue;
                mpz_combit(tmp_p2, p2);

                for (h3 = 0; h3 < smax; h3++) {
                  p3h_steps++;
                  mpz_set(tmp_h3, tmp_p2);
                  if (mpz_tstbit(tmp_h3, h3) == 0)
                    continue;
                  mpz_combit(tmp_h3, h3);

                  for (p3 = 0; p3 < smax; p3++) {
                    p3h_steps++;
                    mpz_set(tmp_p3, tmp_h3);
                    if (mpz_tstbit(tmp_p3, p3) == 1)
                      continue;
                    mpz_combit(tmp_p3, p3);

                    for (h4 = 0; h4 < smax; h4++) {
                      p3h_steps++;
                      mpz_set(tmp_h4, tmp_p3);
                      if (mpz_tstbit(tmp_h4, h4) == 0)
                        continue;
                      mpz_combit(tmp_h4, h4);

                      for (p4 = 0; p4 < smax; p4++) {
                        p3h_steps++;
                        mpz_set(tmp_p4, tmp_h4);
                        if (mpz_tstbit(tmp_p4, p4) == 1)
                          continue;
                        mpz_combit(tmp_p4, p4);

                        for (h5 = 0; h5 < smax; h5++) {
                          p3h_steps++;
                          mpz_set(tmp_h5, tmp_p4);
                          if (mpz_tstbit(tmp_h5, h5) == 0)
                            continue;
                          mpz_combit(tmp_h5, h5);

                          for (p5 = 0; p5 < smax; p5++) {
                            p3h_steps++;
                            mpz_set(tmp_p5, tmp_h5);
                            if (mpz_tstbit(tmp_p5, p5) == 1)
                              continue;
                            mpz_combit(tmp_p5, p5);

                            if (mpz_cmp(tmp_p5, phib[mu]) == 0) {
                              p3h_contr++;

                              int perm = count_bits_right(tmp_h1, h1) +
                                         count_bits_right(tmp_p1, p1) +
                                         count_bits_right(tmp_h2, h2) +
                                         count_bits_right(tmp_p2, p2) +
                                         count_bits_right(tmp_h3, h3) +
                                         count_bits_right(tmp_p3, p3) +
                                         count_bits_right(tmp_h4, h4) +
                                         count_bits_right(tmp_p4, p4) +
                                         count_bits_right(tmp_h5, h5) +
                                         count_bits_right(tmp_p5, p5);
                              if (perm % 2 == 0) {
                                p3h_pos++;
                              } else {
                                p3h_neg++;
                              }
                            } // if
                          } // r
                        }
                      }
                    }
                  }
                } // s
              } // c
            } // d
          } // a
        } // b
      } // if ph_num
      if (ph_num == 4) {
        p2h++;
        for (h1 = 0; h1 < smax; h1++) {
          p2h_steps++;
          mpz_set(tmp_h1, phib[nu]);
          if (mpz_tstbit(tmp_h1, h1) == 0)
            continue;
          mpz_combit(tmp_h1, h1);

          for (p1 = 0; p1 < smax; p1++) {
            p2h_steps++;
            mpz_set(tmp_p1, tmp_h1);
            if (mpz_tstbit(tmp_p1, p1) == 1)
              continue;
            mpz_combit(tmp_p1, p1);

            for (h2 = 0; h2 < smax; h2++) {
              p2h_steps++;
              mpz_set(tmp_h2, tmp_p1);
              if (mpz_tstbit(tmp_h2, h2) == 0)
                continue;
              mpz_combit(tmp_h2, h2);

              for (p2 = 0; p2 < smax; p2++) {
                p2h_steps++;
                mpz_set(tmp_p2, tmp_h2);
                if (mpz_tstbit(tmp_p2, p2) == 1)
                  continue;
                mpz_combit(tmp_p2, p2);

                for (h3 = 0; h3 < smax; h3++) {
                  p2h_steps++;
                  mpz_set(tmp_h3, tmp_p2);
                  if (mpz_tstbit(tmp_h3, h3) == 0)
                    continue;
                  mpz_combit(tmp_h3, h3);

                  for (p3 = 0; p3 < smax; p3++) {
                    p2h_steps++;
                    mpz_set(tmp_p3, tmp_h3);
                    if (mpz_tstbit(tmp_p3, p3) == 1)
                      continue;
                    mpz_combit(tmp_p3, p3);

                    for (h4 = 0; h4 < smax; h4++) {
                      p2h_steps++;
                      mpz_set(tmp_h4, tmp_p3);
                      if (mpz_tstbit(tmp_h4, h4) == 0)
                        continue;
                      mpz_combit(tmp_h4, h4);

                      for (p4 = 0; p4 < smax; p4++) {
                        p2h_steps++;
                        mpz_set(tmp_p4, tmp_h4);
                        if (mpz_tstbit(tmp_p4, p4) == 1)
                          continue;
                        mpz_combit(tmp_p4, p4);

                        for (h5 = 0; h5 < smax; h5++) {
                          p2h_steps++;
                          mpz_set(tmp_h5, tmp_p4);
                          if (mpz_tstbit(tmp_h5, h5) == 0)
                            continue;
                          mpz_combit(tmp_h5, h5);

                          for (p5 = 0; p5 < smax; p5++) {
                            p2h_steps++;
                            mpz_set(tmp_p5, tmp_h5);
                            if (mpz_tstbit(tmp_p5, p5) == 1)
                              continue;
                            mpz_combit(tmp_p5, p5);

                            if (mpz_cmp(tmp_p5, phib[mu]) == 0) {
                              p2h_contr++;

                              int perm = count_bits_right(tmp_h1, h1) +
                                         count_bits_right(tmp_p1, p1) +
                                         count_bits_right(tmp_h2, h2) +
                                         count_bits_right(tmp_p2, p2) +
                                         count_bits_right(tmp_h3, h3) +
                                         count_bits_right(tmp_p3, p3) +
                                         count_bits_right(tmp_h4, h4) +
                                         count_bits_right(tmp_p4, p4) +
                                         count_bits_right(tmp_h5, h5) +
                                         count_bits_right(tmp_p5, p5);
                              if (perm % 2 == 0) {
                                p2h_pos++;
                              } else {
                                p2h_neg++;
                              }
                            } // if
                          } // r
                        }
                      }
                    }
                  }
                } // s
              } // c
            } // d
          } // a
        } // b
      } // if ph_num
      if (ph_num == 2) {
        p1h++;
        for (h1 = 0; h1 < smax; h1++) {
          p1h_steps++;
          mpz_set(tmp_h1, phib[nu]);
          if (mpz_tstbit(tmp_h1, h1) == 0)
            continue;
          mpz_combit(tmp_h1, h1);

          for (p1 = 0; p1 < smax; p1++) {
            p1h_steps++;
            mpz_set(tmp_p1, tmp_h1);
            if (mpz_tstbit(tmp_p1, p1) == 1)
              continue;
            mpz_combit(tmp_p1, p1);

            for (h2 = 0; h2 < smax; h2++) {
              p1h_steps++;
              mpz_set(tmp_h2, tmp_p1);
              if (mpz_tstbit(tmp_h2, h2) == 0)
                continue;
              mpz_combit(tmp_h2, h2);

              for (p2 = 0; p2 < smax; p2++) {
                p1h_steps++;
                mpz_set(tmp_p2, tmp_h2);
                if (mpz_tstbit(tmp_p2, p2) == 1)
                  continue;
                mpz_combit(tmp_p2, p2);

                for (h3 = 0; h3 < smax; h3++) {
                  p1h_steps++;
                  mpz_set(tmp_h3, tmp_p2);
                  if (mpz_tstbit(tmp_h3, h3) == 0)
                    continue;
                  mpz_combit(tmp_h3, h3);

                  for (p3 = 0; p3 < smax; p3++) {
                    p1h_steps++;
                    mpz_set(tmp_p3, tmp_h3);
                    if (mpz_tstbit(tmp_p3, p3) == 1)
                      continue;
                    mpz_combit(tmp_p3, p3);

                    for (h4 = 0; h4 < smax; h4++) {
                      p1h_steps++;
                      mpz_set(tmp_h4, tmp_p3);
                      if (mpz_tstbit(tmp_h4, h4) == 0)
                        continue;
                      mpz_combit(tmp_h4, h4);

                      for (p4 = 0; p4 < smax; p4++) {
                        p1h_steps++;
                        mpz_set(tmp_p4, tmp_h4);
                        if (mpz_tstbit(tmp_p4, p4) == 1)
                          continue;
                        mpz_combit(tmp_p4, p4);

                        for (h5 = 0; h5 < smax; h5++) {
                          p1h_steps++;
                          mpz_set(tmp_h5, tmp_p4);
                          if (mpz_tstbit(tmp_h5, h5) == 0)
                            continue;
                          mpz_combit(tmp_h5, h5);

                          for (p5 = 0; p5 < smax; p5++) {
                            p1h_steps++;
                            mpz_set(tmp_p5, tmp_h5);
                            if (mpz_tstbit(tmp_p5, p5) == 1)
                              continue;
                            mpz_combit(tmp_p5, p5);

                            if (mpz_cmp(tmp_p5, phib[mu]) == 0) {
                              p1h_contr++;

                              int perm = count_bits_right(tmp_h1, h1) +
                                         count_bits_right(tmp_p1, p1) +
                                         count_bits_right(tmp_h2, h2) +
                                         count_bits_right(tmp_p2, p2) +
                                         count_bits_right(tmp_h3, h3) +
                                         count_bits_right(tmp_p3, p3) +
                                         count_bits_right(tmp_h4, h4) +
                                         count_bits_right(tmp_p4, p4) +
                                         count_bits_right(tmp_h5, h5) +
                                         count_bits_right(tmp_p5, p5);
                              if (perm % 2 == 0) {
                                p1h_pos++;
                              } else {
                                p1h_neg++;
                              }
                            } // if
                          } // r
                        }
                      }
                    }
                  }
                } // s
              } // c
            } // d
          } // a
        } // b
      }
      if (ph_num == 0) {
        p0h++;
        for (h1 = 0; h1 < smax; h1++) {
          p0h_steps++;
          mpz_set(tmp_h1, phib[nu]);
          if (mpz_tstbit(tmp_h1, h1) == 0)
            continue;
          mpz_combit(tmp_h1, h1);

          for (p1 = 0; p1 < smax; p1++) {
            p0h_steps++;
            mpz_set(tmp_p1, tmp_h1);
            if (mpz_tstbit(tmp_p1, p1) == 1)
              continue;
            mpz_combit(tmp_p1, p1);

            for (h2 = 0; h2 < smax; h2++) {
              p0h_steps++;
              mpz_set(tmp_h2, tmp_p1);
              if (mpz_tstbit(tmp_h2, h2) == 0)
                continue;
              mpz_combit(tmp_h2, h2);

              for (p2 = 0; p2 < smax; p2++) {
                p0h_steps++;
                mpz_set(tmp_p2, tmp_h2);
                if (mpz_tstbit(tmp_p2, p2) == 1)
                  continue;
                mpz_combit(tmp_p2, p2);

                for (h3 = 0; h3 < smax; h3++) {
                  p0h_steps++;
                  mpz_set(tmp_h3, tmp_p2);
                  if (mpz_tstbit(tmp_h3, h3) == 0)
                    continue;
                  mpz_combit(tmp_h3, h3);

                  for (p3 = 0; p3 < smax; p3++) {
                    p0h_steps++;
                    mpz_set(tmp_p3, tmp_h3);
                    if (mpz_tstbit(tmp_p3, p3) == 1)
                      continue;
                    mpz_combit(tmp_p3, p3);

                    for (h4 = 0; h4 < smax; h4++) {
                      p0h_steps++;
                      mpz_set(tmp_h4, tmp_p3);
                      if (mpz_tstbit(tmp_h4, h4) == 0)
                        continue;
                      mpz_combit(tmp_h4, h4);

                      for (p4 = 0; p4 < smax; p4++) {
                        p0h_steps++;
                        mpz_set(tmp_p4, tmp_h4);
                        if (mpz_tstbit(tmp_p4, p4) == 1)
                          continue;
                        mpz_combit(tmp_p4, p4);

                        for (h5 = 0; h5 < smax; h5++) {
                          p0h_steps++;
                          mpz_set(tmp_h5, tmp_p4);
                          if (mpz_tstbit(tmp_h5, h5) == 0)
                            continue;
                          mpz_combit(tmp_h5, h5);

                          for (p5 = 0; p5 < smax; p5++) {
                            p0h_steps++;
                            mpz_set(tmp_p5, tmp_h5);
                            if (mpz_tstbit(tmp_p5, p5) == 1)
                              continue;
                            mpz_combit(tmp_p5, p5);

                            if (mpz_cmp(tmp_p5, phib[mu]) == 0) {
                              p0h_contr++;

                              int perm = count_bits_right(tmp_h1, h1) +
                                         count_bits_right(tmp_p1, p1) +
                                         count_bits_right(tmp_h2, h2) +
                                         count_bits_right(tmp_p2, p2) +
                                         count_bits_right(tmp_h3, h3) +
                                         count_bits_right(tmp_p3, p3) +
                                         count_bits_right(tmp_h4, h4) +
                                         count_bits_right(tmp_p4, p4) +
                                         count_bits_right(tmp_h5, h5) +
                                         count_bits_right(tmp_p5, p5);
                              if (perm % 2 == 0) {
                                p0h_pos++;
                              } else {
                                p0h_neg++;
                              }
                            } // if
                          } // r
                        }
                      }
                    }
                  }
                } // s
              } // c
            } // d
          } // a
        } // b
      }
    }
  }
  mpz_clears(tmp_h1,
             tmp_p1,
             tmp_h2,
             tmp_p2,
             tmp_h3,
             tmp_p3,
             tmp_h4,
             tmp_p4,
             tmp_h5,
             tmp_p5,
             NULL);

  int cases = p5h + p4h + p3h + p2h + p1h + p0h;
  int contributions =
    p5h_contr + p4h_contr + p3h_contr + p2h_contr + p1h_contr + p0h_contr;
  unsigned long steps =
    p5h_steps + p4h_steps + p3h_steps + p2h_steps + p1h_steps + p0h_steps;

  printf("### Brute Force Evaluation 5p5h ###\n");
  printf("\n");

  printf("Number of Cases:          %d  =>  100%%\n", cases);
  printf("Number of Contributions:  %d  =>  100%%\n", contributions);
  printf("Number of Steps:          %lu  =>  100%%\n", steps);
  printf("\n");

  printf("  Number of 0p0h-Cases:         %d  =>  %f%%\n",
         p0h,
         (float)p0h / cases * 100);
  printf("  Number of 0p0h-Contributions: %d  =>  %f%%\n",
         p0h_contr,
         (float)p0h_contr / contributions * 100);
  printf("  Number of 0p0h-Steps:         %lu  =>  %f%%\n",
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
  printf("  Number of 2p2h-Steps:         %lu  =>  %f%%\n",
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
  printf("  Number of 3p3h-Steps:         %lu  =>  %f%%\n",
         p3h_steps,
         (float)p3h_steps / steps * 100);
  printf("    Number of +1 3p3h-Contributions: %d\n", p3h_pos);
  printf("    Number of -1 3p3h-Contributions: %d\n", p3h_neg);
  printf("\n");

  printf("  Number of 4p4h-Cases:         %d  =>  %f%%\n",
         p4h,
         (float)p4h / cases * 100);
  printf("  Number of 4p4h-Contributions: %d  =>  %f%%\n",
         p4h_contr,
         (float)p4h_contr / contributions * 100);
  printf("  Number of 4p4h-Steps:         %lu  =>  %f%%\n",
         p4h_steps,
         (float)p4h_steps / steps * 100);
  printf("    Number of +1 4p4h-Contributions: %d\n", p4h_pos);
  printf("    Number of -1 4p4h-Contributions: %d\n", p4h_neg);
  printf("\n");

  printf("  Number of 5p5h-Cases:         %d  =>  %f%%\n",
         p5h,
         (float)p5h / cases * 100);
  printf("  Number of 5p5h-Contributions: %d  =>  %f%%\n",
         p5h_contr,
         (float)p5h_contr / contributions * 100);
  printf("  Number of 5p5h-Steps:         %lu  =>  %f%%\n",
         p5h_steps,
         (float)p5h_steps / steps * 100);
  printf("    Number of +1 5p5h-Contributions: %d\n", p5h_pos);
  printf("    Number of -1 5p5h-Contributions: %d\n", p5h_neg);
  printf("\n");
}

#include "../include/bruteForce_eval.h"
#include "../include/bruteForce_printToFile.h"
#include "../include/printb.h"
#include "../include/symm_eval.h"
#include "../include/symm_printToFile.h"
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
toNum (mpz_t result, mpz_t statesb[], int phi[], int kmax)
{
  mpz_t phi_num;
  mpz_init (phi_num);
  for (int i = 0; i < kmax; i++)
    {
      mpz_add (phi_num, phi_num, statesb[phi[i]]);
    }
  mpz_set (result, phi_num);
}

void
toStatesB (mpz_t result, int state)
{
  mpz_set_ui (result, 0);
  mpz_setbit (result, state);
}

int
main ()
{

  int smax = 32; // Number of SP states
  int states[smax];
  for (int i = 0; i < smax; i++)
    {
      states[i] = i;
    }
  mpz_t statesB[smax];
  for (int i = 0; i < smax; i++)
    {
      mpz_init (statesB[i]);
      toStatesB (statesB[i], states[i]);
    }

  int kmax = 4; // Number of Particles
  int Mref = 4; // Size of Mref
  int phi0[] = { 0, 1, 2, 3 };
  int phi2[] = { 0, 1, 2, 4 };
  int phi4[] = { 0, 1, 4, 5 };
  int phi6[] = { 0, 4, 5, 6 };
  int **phi = malloc (Mref * sizeof (int *));
  for (int i = 0; i < Mref; i++)
    {
      phi[i] = malloc (kmax * sizeof (int));
    }
  memcpy (phi[0], phi0, kmax * sizeof (int));
  memcpy (phi[1], phi2, kmax * sizeof (int));
  memcpy (phi[2], phi4, kmax * sizeof (int));
  memcpy (phi[3], phi6, kmax * sizeof (int));

  mpz_t phi0b, phi2b, phi4b, phi6b;
  mpz_inits (phi0b, phi2b, phi4b, phi6b);

  toNum (phi0b, statesB, phi0, kmax);
  toNum (phi2b, statesB, phi2, kmax);
  toNum (phi4b, statesB, phi4, kmax);
  toNum (phi6b, statesB, phi6, kmax);

  mpz_t phib[4];
  for (int i = 0; i < 4; i++)
    {
      mpz_init (phib[i]);
    }

  mpz_set (phib[0], phi0b);
  mpz_set (phib[1], phi2b);
  mpz_set (phib[2], phi4b);
  mpz_set (phib[3], phi6b);

  // FILE* f1 = fopen("bruteForceBitEval_smax33_kmax4.txt", "w");
  // bruteForce_printToFile(f1, states, statesb, smax, phi, phib, kmax, Mref);
  // fclose(f1);
  // FILE* f2 = fopen("symmBitEval_smax8_kmax4.txt", "w");
  // symm_printToFile_2p2h(f2, states, statesb, smax, phi, phib, kmax, Mref);
  // fclose(f2);
  //
  bruteForce_eval_2p2h (states, statesB, smax, phi, phib, kmax, Mref);
  symm_eval_2p2h (states, statesB, smax, phi, phib, kmax, Mref);

  // FILE* f3 = fopen("bruteForceBitEval_smax8_kmax4_3p3h.txt", "w");
  // bruteForce_printToFile_3p3h(f3, states, statesb, smax, phi, phib, kmax,
  // Mref); fclose(f3);
  //
  // bruteForce_eval_3p3h(states, statesb, smax, phi, phib, kmax, Mref);
  // symm_eval_3p3h(states, statesb, smax, phi, phib, kmax, Mref);

  free (phi);
}

#include "../include/bruteForce_printToFile.h"
#include "../include/symm_printToFile.h"
#include "../include/bruteForce_eval.h"
#include "../include/printb.h"
#include "../include/symm_eval.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int toNum(unsigned long long statesb[], int phi[], int kmax) {
  int phi_num = 0;
  for (int i = 0; i < kmax; i++) {
    phi_num += statesb[phi[i]];
  }
  return phi_num;
}

int main() {

  int smax = 64; // Number of SP states
  int states[smax];
  for (int i = 0; i < smax; i++) {
    states[i] = i;
  }
  unsigned long long statesb[smax];
  for (int i = 0; i < smax; i++) {
    statesb[i] = pow(2, i);
  }

  int kmax = 4; // Number of Particles
  int Mref = 4;  // Size of Mref
  int phi0[] = {0, 1, 2, 3};
  int phi2[] = {0, 1, 2, 4};
  int phi4[] = {0, 1, 4, 5};
  int phi6[] = {0, 4, 5, 6};
  int **phi = malloc(Mref * sizeof(int*));
  for (int i = 0; i < Mref; i++) {
    phi[i] = malloc(kmax * sizeof(int));
  }
  memcpy(phi[0], phi0, kmax * sizeof(int));
  memcpy(phi[1], phi2, kmax * sizeof(int));
  memcpy(phi[2], phi4, kmax * sizeof(int));
  memcpy(phi[3], phi6, kmax * sizeof(int));

  int phi0b = toNum(statesb, phi0, kmax);
  int phi2b = toNum(statesb, phi2, kmax);
  int phi4b = toNum(statesb, phi4, kmax);
  int phi6b = toNum(statesb, phi6, kmax);

  int phib[] = {phi0b, phi2b, phi4b, phi6b};

  // FILE* f1 = fopen("bruteForceBitEval_smax33_kmax4.txt", "w");
  // bruteForce_printToFile(f1, states, statesb, smax, phi, phib, kmax, Mref);
  // fclose(f1);
  // FILE* f2 = fopen("symmBitEval_smax8_kmax4.txt", "w");
  // symm_printToFile_2p2h(f2, states, statesb, smax, phi, phib, kmax, Mref);
  // fclose(f2);
  // 
  bruteForce_eval_2p2h(states, statesb, smax, phi, phib, kmax, Mref);
  symm_eval_2p2h(states, statesb, smax, phi, phib, kmax, Mref);

  // FILE* f3 = fopen("bruteForceBitEval_smax8_kmax4_3p3h.txt", "w");
  // bruteForce_printToFile_3p3h(f3, states, statesb, smax, phi, phib, kmax, Mref);
  // fclose(f3);
  //
  // bruteForce_eval_3p3h(states, statesb, smax, phi, phib, kmax, Mref);
  // symm_eval_3p3h(states, statesb, smax, phi, phib, kmax, Mref);

  free(phi);

}

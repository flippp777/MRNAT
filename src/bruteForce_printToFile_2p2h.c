#include "../include/bruteForce_printToFile.h"
#include "../include/insertionSort_perm.h"
#include "../include/printb.h"
#include <math.h>
#include <stdio.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void
bruteForce_printToFile_2p2h(FILE* fptr,
                            int states[],
                            unsigned long long statesb[],
                            int smax,
                            int** phi,
                            int phib[],
                            int kmax,
                            int Mref)
{
  int mu; // Index for SD mu
  int nu; // Index for SD nu
  int a;  // Index for particle a
  int r;  // Index for paritcle r
  int b;  // Index for hole b
  int s;  // Index for hole s
  for (mu = 0; mu < Mref; mu++) {
    // Print current SD mu
    fprintf(fptr, "########## phi_mu = ");
    fprintb(fptr, phib[mu], smax);
    fprintf(fptr, " ##########");
    fprintf(fptr, "\n");
    fprintf(fptr, "\n");
    for (nu = 0; nu < Mref; nu++) {

      for (b = 0; b < smax; b++) {
        unsigned long long tmp_b = phib[nu];
        if ((tmp_b & statesb[b]) != statesb[b])
          continue;
        tmp_b ^= statesb[b];

        for (a = 0; a < smax; a++) {
          unsigned long long tmp_a = tmp_b;
          if ((tmp_a & statesb[a]) != 0)
            continue;
          tmp_a ^= statesb[a];

          for (s = 0; s < smax; s++) {
            unsigned long long tmp_s = tmp_a;
            if ((tmp_s & statesb[s]) != statesb[s])
              continue;
            tmp_s ^= statesb[s];

            for (r = 0; r < smax; r++) {
              unsigned long long tmp_r = tmp_s;
              if ((tmp_r & statesb[r]) != 0)
                continue;
              tmp_r ^= statesb[r];

              if (tmp_r == phib[mu]) {

                // phi_nu with its 2p-2h "excitations"
                int phi_nu[kmax];
                for (int i = 0; i < kmax; i++) {
                  if (phi[nu][i] == states[b]) {
                    phi_nu[i] = states[a];
                  } else {
                    phi_nu[i] = phi[nu][i];
                  }
                }
                for (int i = 0; i < kmax; i++) {
                  if (phi_nu[i] == states[s]) {
                    phi_nu[i] = states[r];
                    break;
                  }
                }

                fprintf(fptr, "phi_mu:           ");
                fprintb(fptr, phib[mu], smax);
                fprintf(fptr, " = ");
                for (int i = 0; i < kmax; i++) {
                  fprintf(fptr, "%d", phi[mu][i]);
                }
                fprintf(fptr, "\n");

                fprintf(fptr, "phi_nu:           ");
                fprintb(fptr, phib[nu], smax);
                fprintf(fptr, " = ");
                for (int i = 0; i < kmax; i++) {
                  fprintf(fptr, "%d", phi[nu][i]);
                }
                fprintf(fptr, "\n");

                fprintf(fptr, "Hole b:           ");
                fprintb(fptr, statesb[b], smax);
                fprintf(fptr, " = %d", states[b]);
                fprintf(fptr, "\n");

                fprintf(fptr, "Particle a:       ");
                fprintb(fptr, statesb[a], smax);
                fprintf(fptr, " = %d", states[a]);
                fprintf(fptr, "\n");

                fprintf(fptr, "Hole s:           ");
                fprintb(fptr, statesb[s], smax);
                fprintf(fptr, " = %d", states[s]);
                fprintf(fptr, "\n");

                fprintf(fptr, "Particle r:       ");
                fprintb(fptr, statesb[r], smax);
                fprintf(fptr, " = %d", states[r]);
                fprintf(fptr, "\n");
                fprintf(fptr, "\n");

                fprintf(fptr, "phi_nu^{ar}_{bs}: ");
                fprintb(fptr, tmp_r, smax);
                fprintf(fptr, "\n");
                fprintf(fptr, "\n");
                fprintf(fptr, "phi_mu =           {");
                for (int i = 0; i < kmax; i++) {
                  fprintf(fptr, "%d", phi[mu][i]);
                }
                fprintf(fptr, "}\n");
                fprintf(fptr, "phi_nu^{ar}_{bs} = {");
                for (int i = 0; i < kmax; i++) {
                  fprintf(fptr, "%d", phi_nu[i]);
                }
                fprintf(fptr, "}\n");

                // Track permutations and sorting for sign
                int perm = 0;
                insertionSort_perm(phi_nu, kmax, &perm);
                fprintf(fptr, "Permutations: %d", perm);
                fprintf(fptr, "\n");

                fprintf(fptr, "<phi_mu|phi_nu^{ar}_{bs}> = %f", pow(-1, perm));
                fprintf(fptr, "\n");
                fprintf(fptr, "\n");
              }
            }
          }
        }
      }
    }
  }
}

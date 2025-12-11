#include "../include/symm_printToFile.h"
#include "../include/insertionSort_perm.h"
#include "../include/numberOfStates.h"
#include "../include/printb.h"
#include <math.h>
#include <stdio.h>

// Look for all 2p-2h excitations for <phi_mu|phi_nu^{ar}_{bs}
void symm_printToFile_2p2h(FILE *fptr, int states[], unsigned long long statesb[],
                      int smax, int **phi, int phib[], int kmax, int Mref) {
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
      int ph = phib[mu] ^ phib[nu];
      int ph_num; // Number #p#h excitations between SD's
      ph_num = numberOfStates(ph);
      int phi_nu_b = phib[nu];

      // 2p2h difference
      if (ph_num == 4) {

        int ph_states[4];
        int ph_statesb[4];
        int ph_states_idx = 0;
        int ph_count = 0;
        for (int i = 0; i < smax; i++) {
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
        int p_statesb[2];
        int p_states_idx = 0;
        int h_states[2];
        int h_statesb[2];
        int h_states_idx = 0;
        for (int i = 0; i < 4; i++) {
          if ((ph_statesb[i] & phib[nu]) == ph_statesb[i]) {
            h_states[h_states_idx] = ph_states[i];
            h_statesb[h_states_idx] = ph_statesb[i];
            h_states_idx++;
          } else {
            p_states[p_states_idx] = ph_states[i];
            p_statesb[p_states_idx] = ph_statesb[i];
            p_states_idx++;
          }
        }

        // 4 contributions
        // phi_nu with its 2p-2h "excitations"
        for (int h = 0; h < 2; h++) {
          int phi_nu[kmax];
          for (int p = 0; p < 2; p++) {
            for (int k = 0; k < kmax; k++) {
              if (phi[nu][k] == h_states[h]) {
                phi_nu[k] = p_states[p];
                phi_nu_b = phi_nu_b ^ statesb[h];
                phi_nu_b = phi_nu_b ^ statesb[p];
              } else if (phi[nu][k] == h_states[(h + 1) % 2]) {
                phi_nu[k] = p_states[(p + 1) % 2];
                phi_nu_b = phi_nu_b ^ statesb[(h + 1) % 2];
                phi_nu_b = phi_nu_b ^ statesb[(p + 1) % 2];
              } else {
                phi_nu[k] = phi[nu][k];
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
            fprintb(fptr, h_statesb[h], smax);
            fprintf(fptr, " = %d", h_states[h]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle a:       ");
            fprintb(fptr, p_statesb[p], smax);
            fprintf(fptr, " = %d", p_states[p]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Hole s:           ");
            fprintb(fptr, h_statesb[(h + 1) % 2], smax);
            fprintf(fptr, " = %d", h_states[(h + 1) % 2]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle r:       ");
            fprintb(fptr, p_statesb[(p + 1) % 2], smax);
            fprintf(fptr, " = %d", p_states[(p + 1) % 2]);
            fprintf(fptr, "\n");
            fprintf(fptr, "\n");

            fprintf(fptr, "phi_nu^{ar}_{bs}: ");
            fprintb(fptr, phi_nu_b, smax);
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

      // 1p1h difference
      if (ph_num == 2) {

        int ph_states[2];
        int ph_statesb[2];
        int ph_states_idx = 0;
        for (int i = 0; i < smax; i++) {
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
        int p_stateb;
        int h_state;
        int h_stateb;
        for (int i = 0; i < 2; i++) {
          if ((ph_statesb[i] & phib[nu]) == ph_statesb[i]) {
            h_state = ph_states[i];
            h_stateb = ph_statesb[i];
          } else {
            p_state = ph_states[i];
            p_stateb = ph_statesb[i];
          }
        }

        // the c_i / c^dagger_j order implemented
        // is c^d_r c_s c^d_a c_b
        // case: b=a, r=p, s=h
        //  so all cases where a h-state can be
        //  annihilated and then a p-state created
        int phi_nu[kmax];
        for (int k = 0; k < kmax; k++) {
          if (phi[nu][k] == h_state) {
            phi_nu[k] = p_state;
            phi_nu_b = phi_nu_b ^ h_stateb;
            phi_nu_b = phi_nu_b ^ p_stateb;
          } else {
            phi_nu[k] = phi[nu][k];
          }
        }
        int phib_nu = phib[nu];

        int idx = 0;
        while (phib_nu) {
          if ((phib_nu & 1) == 1) {
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
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle a:       ");
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Hole s:           ");
            fprintb(fptr, h_stateb, smax);
            fprintf(fptr, " = %d", h_state);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle r:       ");
            fprintb(fptr, p_stateb, smax);
            fprintf(fptr, " = %d", p_state);
            fprintf(fptr, "\n");
            fprintf(fptr, "\n");

            fprintf(fptr, "phi_nu^{ar}_{bs}: ");
            fprintb(fptr, phi_nu_b, smax);
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
          idx++;
          phib_nu >>= 1;
        }
        // The same happens for
        // case: s=r, a=p, b=h
        //  It has the same sign as before.
        // (The factor will be different in
        // the final implementation)
        // and now the h- / p-summation states are
        // annihilated / created after the ph-states
        int phib_mu = phib[mu];

        idx = 0;
        while (phib_mu) {
          if ((phib_mu & 1) == 1) {
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
            fprintb(fptr, h_stateb, smax);
            fprintf(fptr, " = %d", h_state);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle a:       ");
            fprintb(fptr, p_stateb, smax);
            fprintf(fptr, " = %d", p_state);
            fprintf(fptr, "\n");

            fprintf(fptr, "Hole s:           ");
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle r:       ");
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");
            fprintf(fptr, "\n");

            fprintf(fptr, "phi_nu^{ar}_{bs}: ");
            fprintb(fptr, phi_nu_b, smax);
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
          idx++;
          phib_mu >>= 1;
        }
        // case: b=r, a=p, s=h
        // Now the summation state is destroyed and
        // the p-state is created in its place, then
        // the h-states is destroyed and the other
        // summation state put in its place.
        // It is the same as before but with the
        // opposite sign.
        int equal_holes = phib[mu] & phib[nu];

        idx = 0;
        while (equal_holes) {
          if ((equal_holes & 1) == 1) {
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
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle a:       ");
            fprintb(fptr, p_stateb, smax);
            fprintf(fptr, " = %d", p_state);
            fprintf(fptr, "\n");

            fprintf(fptr, "Hole s:           ");
            fprintb(fptr, h_stateb, smax);
            fprintf(fptr, " = %d", h_state);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle r:       ");
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");
            fprintf(fptr, "\n");

            fprintf(fptr, "phi_nu^{ar}_{bs}: ");
            fprintb(fptr, phi_nu_b, smax);
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
            int perm = 1;
            insertionSort_perm(phi_nu, kmax, &perm);
            fprintf(fptr, "Permutations: %d", perm);
            fprintf(fptr, "\n");

            fprintf(fptr, "<phi_mu|phi_nu^{ar}_{bs}> = %f", pow(-1, perm));
            fprintf(fptr, "\n");
            fprintf(fptr, "\n");
          }
          idx++;
          equal_holes >>= 1;
        }
        // case: s=a, r=p, b=h
        // Here the summation states start
        // with a p-creation so the contrubutions
        // are determined by the number of p-states
        // the sign should be opposite to the first
        // two cases
        int equal_particles = ~(phib[mu] | phib[nu]);

        idx = 0;
        for (int idx = 0; idx < smax; idx++) {
          if ((equal_particles & 1) == 1) {
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
            fprintb(fptr, h_stateb, smax);
            fprintf(fptr, " = %d", h_state);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle a:       ");
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Hole s:           ");
            fprintb(fptr, statesb[idx], smax);
            fprintf(fptr, " = %d", states[idx]);
            fprintf(fptr, "\n");

            fprintf(fptr, "Particle r:       ");
            fprintb(fptr, p_stateb, smax);
            fprintf(fptr, " = %d", p_state);
            fprintf(fptr, "\n");
            fprintf(fptr, "\n");

            fprintf(fptr, "phi_nu^{ar}_{bs}: ");
            fprintb(fptr, phi_nu_b, smax);
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
          equal_particles >>= 1;
        }
      }
      // 0p0h difference
      if (ph_num == 0) {

        // The signs are now all positive
        unsigned long long phi_nu_b = phib[nu];
        unsigned long long phi_nu_b_tmp = phib[nu];
        b = 0;
        int *phi_nu = phi[nu];
        while (phi_nu_b_tmp) {
          for (s = 0; s < smax; s++) {

            // b=h-state
            if ((phi_nu_b_tmp & 1) == 1) {
              // s=h-state
              //  b=a & s=r
              if ((phi_nu_b & statesb[s]) == statesb[s]) {
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
                fprintb(fptr, statesb[b], smax);
                fprintf(fptr, " = %d", states[b]);
                fprintf(fptr, "\n");

                fprintf(fptr, "Hole s:           ");
                fprintb(fptr, statesb[s], smax);
                fprintf(fptr, " = %d", states[s]);
                fprintf(fptr, "\n");

                fprintf(fptr, "Particle r:       ");
                fprintb(fptr, statesb[s], smax);
                fprintf(fptr, " = %d", states[s]);
                fprintf(fptr, "\n");
                fprintf(fptr, "\n");

                fprintf(fptr, "phi_nu^{ar}_{bs}: ");
                fprintb(fptr, phi_nu_b, smax);
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

                // b <-> s
              //   if (statesb[b] < statesb[s]) {
              //     fprintf(fptr, "phi_mu:           ");
              //     fprintb(fptr, phib[mu], smax);
              //     fprintf(fptr, " = ");
              //     for (int i = 0; i < kmax; i++) {
              //       fprintf(fptr, "%d", phi[mu][i]);
              //     }
              //     fprintf(fptr, "\n");
              //
              //     fprintf(fptr, "phi_nu:           ");
              //     fprintb(fptr, phib[nu], smax);
              //     fprintf(fptr, " = ");
              //     for (int i = 0; i < kmax; i++) {
              //       fprintf(fptr, "%d", phi[nu][i]);
              //     }
              //     fprintf(fptr, "\n");
              //
              //     fprintf(fptr, "Hole b:           ");
              //     fprintb(fptr, statesb[s], smax);
              //     fprintf(fptr, " = %d", states[s]);
              //     fprintf(fptr, "\n");
              //
              //     fprintf(fptr, "Particle a:       ");
              //     fprintb(fptr, statesb[s], smax);
              //     fprintf(fptr, " = %d", states[s]);
              //     fprintf(fptr, "\n");
              //
              //     fprintf(fptr, "Hole s:           ");
              //     fprintb(fptr, statesb[b], smax);
              //     fprintf(fptr, " = %d", states[b]);
              //     fprintf(fptr, "\n");
              //
              //     fprintf(fptr, "Particle r:       ");
              //     fprintb(fptr, statesb[b], smax);
              //     fprintf(fptr, " = %d", states[b]);
              //     fprintf(fptr, "\n");
              //     fprintf(fptr, "\n");
              //
              //     fprintf(fptr, "phi_nu^{ar}_{bs}: ");
              //     fprintb(fptr, phi_nu_b, smax);
              //     fprintf(fptr, "\n");
              //     fprintf(fptr, "\n");
              //     fprintf(fptr, "phi_mu =           {");
              //     for (int i = 0; i < kmax; i++) {
              //       fprintf(fptr, "%d", phi[mu][i]);
              //     }
              //     fprintf(fptr, "}\n");
              //     fprintf(fptr, "phi_nu^{ar}_{bs} = {");
              //     for (int i = 0; i < kmax; i++) {
              //       fprintf(fptr, "%d", phi_nu[i]);
              //     }
              //     fprintf(fptr, "}\n");
              //
              //     // Track permutations and sorting for sign
              //     perm = 0;
              //     insertionSort_perm(phi_nu, kmax, &perm);
              //     fprintf(fptr, "Permutations: %d", perm);
              //     fprintf(fptr, "\n");
              //
              //     fprintf(fptr, "<phi_mu|phi_nu^{ar}_{bs}> = %f",
              //             pow(-1, perm));
              //     fprintf(fptr, "\n");
              //     fprintf(fptr, "\n");
              //   }
              }
              // s=p-state
              //  b=r & s=a
              else {
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
                fprintb(fptr, statesb[s], smax);
                fprintf(fptr, " = %d", states[s]);
                fprintf(fptr, "\n");

                fprintf(fptr, "Hole s:           ");
                fprintb(fptr, statesb[s], smax);
                fprintf(fptr, " = %d", states[s]);
                fprintf(fptr, "\n");

                fprintf(fptr, "Particle r:       ");
                fprintb(fptr, statesb[b], smax);
                fprintf(fptr, " = %d", states[b]);
                fprintf(fptr, "\n");
                fprintf(fptr, "\n");

                fprintf(fptr, "phi_nu^{ar}_{bs}: ");
                fprintb(fptr, phi_nu_b, smax);
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
          b++;
          phi_nu_b_tmp >>= 1;
        }
      }
    }
  }
}

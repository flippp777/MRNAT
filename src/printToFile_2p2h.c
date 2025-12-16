#include "../include/printToFile.h"
#include "../include/excite.h"
#include "../include/printb.h"
#include <gmp.h>
#include <stdio.h>

void
printToFile_2p2h(FILE* fptr,
            mpz_t phib[],
            int smax,
            int** phi,
            int kmax,
            int mu,
            int nu,
            int b,
            int a,
            int s,
            int r,
            mpz_t phi_nub_exc,
            int perm,
            int sign)
{
  // phi_nu with its 2p-2h "excitations"
  int phi_nu[kmax];
  excite_2p2h(phi_nu, phi, nu, kmax, b, a, s, r);

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

  mpz_t stateB;
  mpz_init(stateB);

  fprintf(fptr, "Hole b:           ");
  mpz_setbit(stateB, b);
  fprintb(fptr, stateB, smax);
  mpz_clrbit(stateB, b);
  fprintf(fptr, " = %d", b);
  fprintf(fptr, "\n");

  fprintf(fptr, "Particle a:       ");
  mpz_setbit(stateB, a);
  fprintb(fptr, stateB, smax);
  mpz_clrbit(stateB, a);
  fprintf(fptr, " = %d", a);
  fprintf(fptr, "\n");

  fprintf(fptr, "Hole s:           ");
  mpz_setbit(stateB, s);
  fprintb(fptr, stateB, smax);
  mpz_clrbit(stateB, s);
  fprintf(fptr, " = %d", s);
  fprintf(fptr, "\n");

  fprintf(fptr, "Particle r:       ");
  mpz_setbit(stateB, r);
  fprintb(fptr, stateB, smax);
  mpz_clrbit(stateB, r);
  fprintf(fptr, " = %d", r);
  fprintf(fptr, "\n");
  fprintf(fptr, "\n");

  fprintf(fptr, "phi_nu^{ar}_{bs}: ");
  fprintb(fptr, phi_nub_exc, smax);
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

  // Permutations for destroying and inserting + sorting

  fprintf(fptr, "Permutations: %d", perm);
  fprintf(fptr, "\n");
  fprintf(fptr, "<phi_mu|phi_nu^{ar}_{bs}> = %d", sign);
  fprintf(fptr, "\n");
  fprintf(fptr, "\n");

  mpz_clear(stateB);
}

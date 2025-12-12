#include "../include/printb.h"
#include <gmp.h>
#include <stdio.h>

void
printb(int n, int smax)
{
  if (n == 0) {
    for (int i = 0; i < smax; i++)
      printf("0");
    return;
  }

  int binary[smax];
  for (int i = 0; i < smax; i++)
    binary[i] = 0;
  int index = 0;

  while (n > 0) {
    binary[index] = n % 2;
    n = n / 2;
    index++;
  }

  for (int i = smax - 1; i >= 0; i--) {
    printf("%d", binary[i]);
  }
}

void
fprintb(FILE* fptr, mpz_t n, int smax)
{
  if (mpz_cmp_ui(n, 0) == 0) {
    for (int i = 0; i < smax; i++)
      fprintf(fptr, "0");
    return;
  }

  mpz_t binary[smax];
  for (int i = 0; i < smax; i++)
    mpz_init(binary[i]);

  int index = 0;
  while (mpz_cmp_ui(n, 0)) {
    mpz_mod_ui(binary[index], n, 2);
    mpz_tdiv_q_ui(n, n, 2);
    index++;
  }

  for (int i = smax - 1; i >= 0; i--) {
    fprintf(fptr, "%d", mpz_get_ui(binary[i]));
  }

  for (int i = 0; i < smax; i++)
    mpz_clear(binary[i]);

}

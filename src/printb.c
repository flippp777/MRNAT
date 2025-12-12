#include "../include/printb.h"
#include <gmp.h>
#include <stdio.h>

void
printb (mpz_t num, int smax)
{
  for (int i = smax - 1; i >= 0; i--)
    {
      printf ("%d", mpz_tstbit (num, i));
    }
}

void
fprintb (FILE *fptr, mpz_t num, int smax)
{
  for (int i = smax - 1; i >= 0; i--)
    {
      fprintf (fptr, "%d", mpz_tstbit (num, i));
    }
}

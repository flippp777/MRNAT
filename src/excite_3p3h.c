#include "../include/excite.h"

void
excite_3p3h (int phi_nu[], int **phi, int nu, int kmax, int b,
             int a, int d, int c, int s, int r)
{
  for (int i = 0; i < kmax; i++)
    {
      if (phi[nu][i] == b)
        {
          phi_nu[i] = a;
        }
      else
        {
          phi_nu[i] = phi[nu][i];
        }
    }
  for (int i = 0; i < kmax; i++)
    {
      if (phi_nu[i] == d)
        {
          phi_nu[i] = c;
          break;
        }
    }
  for (int i = 0; i < kmax; i++)
    {
      if (phi_nu[i] == s)
        {
          phi_nu[i] = r;
          break;
        }
    }
}

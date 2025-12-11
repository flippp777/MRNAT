#include "../include/numberOfStates.h"

int numberOfStates(int phib) {
  int num = 0;
  while (phib) {
    num += phib & 1;
    phib >>= 1;
  }
  return num;
}

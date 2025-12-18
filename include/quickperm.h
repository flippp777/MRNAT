#ifndef QUICKPERM
#define QUICKPERM

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int
quickperm(const int* input, size_t n, int** out, size_t* count);

#endif // !QUICKPERM

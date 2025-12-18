#ifndef QUICKPERM
#define QUICKPERM

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int
quickperm(int* input, size_t n, int** out, size_t* count);

int
quickperm_ptr(int* const* input, size_t n, int*** out, size_t* count);

#endif // !QUICKPERM

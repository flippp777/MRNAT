#ifndef HASH
#define HASH
#include <stddef.h>
#include <stdint.h>

uint64_t
hash_func(const int* a, size_t block, size_t n);

#endif // !HASH

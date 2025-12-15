#ifndef PERMUTE
#define PERMUTE

void
permute(int** perms, int* arr, int start, int n, int perm_count);

// void
// permute_rek(int** perms, int* arr, int start, int n, int perm_count, int* count);

void
permute_ptr(int*** perms, int** arr, int start, int n, int perm_count);

// void
// permute_ptr_rek(int*** perms, int** arr, int start, int n, int perm_count, int* count);


#endif // !PERMUTE

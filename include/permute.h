#ifndef PERMUTE
#define PERMUTE

int**
permute(int* arr, int start, int n, int perm_count);

// void
// permute_rek(int** perms, int* arr, int start, int n, int perm_count, int* count);

int***
permute_ptr(int** arr, int start, int n, int perm_count);

// void
// permute_ptr_rek(int*** perms, int** arr, int start, int n, int perm_count, int* count);


#endif // !PERMUTE

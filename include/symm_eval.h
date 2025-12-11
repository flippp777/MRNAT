#ifndef SYMM_EVAL
#define SYMM_EVAL

void symm_eval_2p2h(int states[], unsigned long long statesb[], int smax, int **phi, int phib[],
               int kmax, int Mref);

void symm_eval_3p3h(int states[], unsigned long long statesb[], int smax, int **phi, int phib[],
               int kmax, int Mref);

#endif // !SYMM_EVAL

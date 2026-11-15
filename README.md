# MRNAT – Multi-Reference Natural Orbitals (prototype)

> **Note:** This repository holds early prototyping and testing code from the start of my Master's thesis. The final implementation was developed separately in a different code base, so this code is not maintained and was not used for any results.

## What it does

The goal was to evaluate, for a multi-reference state built from several Slater determinants, which reference determinants $|\Phi_\mu\rangle$ and $|\Phi_\nu\rangle$ are connected by $n$-particle–$n$-hole excitations. Evaluating these matrix elements is the basis for building density matrices and natural orbitals.

- Slater determinants are stored as **bit strings** using GMP big integers.
- The excitation level between two determinants follows from `XOR` + `popcount`.
- For up to 5p5h excitations, two approaches are implemented and checked against each other:
  - a **brute-force** loop over all single-particle indices
  - a **symmetry-based** evaluation that enumerates only the particle/hole permutations, with phases from bit counting, and removes duplicates with a small custom **hash set** (splitmix64-based)
- `test/` contains standalone tests of the permutation and hashing routines.

## Build

Requires a C compiler and [GMP](https://gmplib.org/).

```bash
gcc -O2 -o mrnat src/*.c -lgmp -lm
./mrnat
```

The number of single-particle states, the particle number and the reference determinants are set directly in `src/main.c`.

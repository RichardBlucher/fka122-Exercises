# FKA121/FIM540 Computational Physics 2022

## C programs for HOMEWORK H1

Three C functions are provided:

- `potential.c`: Contains functions that calculates the energy, forces, and
  virial using an Neuroevolution potential (NEP) for aluminum.
- `lattice.c`: Contains a function that generates the positions of atoms in a
  fcc lattice.
- `run.c`: You should write your MD program in this file.
  This makes it possible to test 100% of your code.
  Since you can't have two main functions in one program.
- `tools.c`: Add your library from C4 here (and the corresponding header file
  under include).

The instruction on how to use the functions in lattice.c and potential.c is
in the corresponding header files.

Run `make` to compile your program.

Additionally in this folder you will find `nep.txt` which defines the parameters
for the NEP (don't touch this file) and `src/nepc/` where the code for the NEP
is implemented (we don't recommend touching this either).

## Unit testing

You have the possibility to unit test your code. Write your unit tests in
`unit-test/src/` and do `make test && ./run-test`.

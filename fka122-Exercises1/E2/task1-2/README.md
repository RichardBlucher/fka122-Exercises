---
header-includes:
- \usepackage[a4paper, total={6in, 9in}]{geometry}
- \usepackage{fancyhdr}
- \pagestyle{fancy}
- \fancyhead[L]{FKA121/FIM540 Computational physics}
- \fancyhead[R]{C4 Writing a simple library}
colorlinks: true
output: pdf_document
---

# Exercise C4a: Writing a simple library

In this task you are supposed to write a small tools file that contains
useful functions. This file can be used through out this course, for
your convenience. The input and output of the functions are both defined
in the tools.c and tools.h. Complete the functions, then run the test to
make sure that you have implemented the functions correctly. The
makefiles and project structure is something that you can reuse
throughout the course.

Functions that are to be written

* `elementwise_addition`
* `addition_with_constant`
* `elementwise_multiplication`
* `multiplication_with_constant`
* `dot_product`
* `create_2D_array`
* `destroy_2D_array`
* `matrix_multiplcation`
* `vector_norm`
* `normalize_vector`
* `average`
* `standard_deviation`
* `distance_between_vectors`
* `cumulative_integration`
* `write_xyz`
* `fft_freq`

look in the header file for descriptions of the functions.
The `fft_freq` should have the same output as [numpy.fft.fftfreq](https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html#numpy.fft.fftfreq) with 2$\pi$.

You can visualize the atoms written to atoms.extxyz with ase, `ase gui atoms.extxyz` (which might be useful for later).
There is also more advanced programs, such as, [ovito](https://www.ovito.org/python-downloads/) and vmd.


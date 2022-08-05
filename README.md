# Parallelization study of a 2D N-body problem
This code was developed by Abhiroop Ghosh as a part of the CMSE 822 Parallel Computing Fall 2020 course at Michigan State University. It aims to compare the performance of serial execution and multiple parallelization techniques (MPI, OpenMP, etc.) on an N-body problem predicting the positions and velocities of celestial objects, considered as point particles for simplicity.

To run the code simply run ```make``` followed by ```./nbody_all.o <total timesteps> <no. of particles>```. Results are stored in an HDF5 file.

A report consisting of the problem details and findings of the study is given in the docs folder.

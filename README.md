------------
Graph 500 Generator
------------------------
This file contains two components (folders). 
- g500_2d_tuple: MPI enabled kronecker generator. It outputs **tuple list**.
- 2d_tuple_to_2d_csr: Reads output from **g500_2d_tuple**, convert it to **csr** and associated **beg_pos**.

---
Compile
-------------
- Machine should install MPI.
- **make** will compile the file to executable binary

----
Run
------------
**./executable** will the required input formats to run the code.
- **Single node**: 
  - mpirun -n **number-processes** -host localhost /path/to/generator_test_mpi  log_numverts degree row-partitions column-partitions
  - mpirun -n num-partitions(row-partitions x col-partitions) -host localhost /path/to/tuple_to_csr log_numverts degree row-partitions column-partitions **number-processes**
- **Distributed node** (cluster): please refer to **run.bash**, which is assuming "slurm" resource management policy.


--------
Contact: Hang_Liu@uml.edu. 

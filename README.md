------------
Graph 500 Generator
------------------------
This file contains two components (folders). 
- g500_2d_tuple: MPI enabled kronecker generator. It outputs **tuple list**.
- 2d_tuple_to_2d_csr: Reads output from **g500_2d_tuple**, convert it to **csr** and associated **beg_pos**.

---
Compile
-------------
- Machine should install MPI (on colonial0, colonial03 and colonial04, as well as colonialone).
- **make** will compile the file to executable binary

--
Run
------------
**./executable** will the required input formats to run the code.
- **colonial0** (local machine): 
  - mpirun -n **number-processes** -host localhost /path/to/generator_test_mpi  log_numverts degree row-partitions column-partitions
  - mpirun -n num-partitions(row-partitions x col-partitions) -host localhost /path/to/tuple_to_csr log_numverts degree row-partitions column-partitions **number-processes**
- **colonialone** (cluster): please refer to **run.bash** script in this folder.


---
Let's talk about the detail functionality of each line, if you need to

--Hang

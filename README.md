# Wave-propagation
Below is just a small exerpt of the full statement [Statement.pdf](Statement.pdf).

>Wave propagation at the surface of the ocean can be modelled using the so-called “shallow
>water” equations, which are derived by depth-integrating the Navier-Stokes equations.
>The domain of study is a rectangle [0, a] × [0, b] whose dimensions and bathymetric map
>(map of h(x, y)) are given by an input file. Sample input files are provided in the directory [`Maps/`]

In the `Scripts` directory are multiple bash scripts and Slurm scripts to run the program on a cluster. 
There are 2 python scripts, `contour_plot.py` to display the input map files and the output files `eta_XXX.dat`, `u_XXX.dat` and `v_XXX.dat` and `plot.py` to plot the scalabilty of the program using data gathered from running the program using the parameter files in `Report/Parameters`.

A full analysis can be seen in the Report [Report.pdf](Report.pdf).

## Getting Started

This can be run on a local machine or on a cluster using a Slurm script.

### Prerequisites

Make sure that you have gcc installed using the command below:
```
sudo apt-get install gcc
```
OpenMP comes with the gcc compiler.

Next, you'll have to install OpenMPI to do multiprocess computation. There are may possible ways, one of which is to build it from source [here](https://www.open-mpi.org/software/ompi/v4.0/).

This project was compiled and run using `gcc (Ubuntu 7.4.0-1ubuntu1~18.04.1) 7.4.0` and opemmpi `mpirun (Open MPI) 2.1.1
` . If you want to make sure it works, I invite you to install the same one.
You can check your version with 
```
gcc --version
```
and
```
mpirun --version
```

### Compilation

#### Locally
To compile the program and run it locally (on your machine), run the following commands in the terminal.

```
mpicc project2_Delaunoy_Crasset*.c -std=c99 -lm -fopenmp -o waves               
```
and run the program using the template 
```
OMP_NUM_THREADS=NB_OF_THREADS mpirun -np NB_PROCESSORS waves PARAMETER_FILE  MAP_FILE FINITE_DIFFERENCE_SCHEME
```

Replace `NB_OF_THREADS` with the number of threads you wish to run it with, `NB_PROCESSORS` with the number of processes you want to launch, `PARAMETER_FILE` with one of the parameters files in direcotry `Parameters/` or `Report/Parameters` and `MAP_FILE` with a bathymetric map file in directory `Maps/`.
Lastly, choose the `FINITE_DIFFERENCE_SCHEME` that you want, 0 for the explicit scheme and 1 for the implicit scheme.

For example,

Running the implicit scheme on the map `sriLanka` using 2 threads and 1 process.
```
OMP_NUM_THREADS=2 mpirun -np 1 waves Parameters/sriLanka_implicit_result.txt  Maps/sriLanka.dat 1
```

Running the explicit scheme on the map `sriLanka` using 4 threads and 2 process.
```
OMP_NUM_THREADS=4 mpirun -np 2 waves Parameters/sriLanka_explicit_result.txt  Maps/sriLanka.dat 0
```

#### On a cluster with Slurm

If you want to fully use the capabilities of the program, we suggest you run this code on a cluster managed by the Slurm Workload Manager.

Sample scripts are provided in the `Scripts/` directory.

## Authors

* **Arnaud Delaunoy ** - *Implementation of the problem, scalability and stability analysis*

* **Tom Crasset ** - *Implementation of the problem, scalability and stability analysis*

* **Anthony Royer ** - *Creation of statement, provided map files and `makeMovie.m` to create a video of the wave propagation.*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

Thanks to Pr. C. Geuzaine for his course of High Perfomance Scientific Computing at the University of Liège as well as his assistant, Anthony Royer.
Moreover, we would like to thank the CECI(Consortium des Équipements de Calcul Intensif) for the cluster they provided.

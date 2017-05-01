# PMC-VMC


## Description

A Projection Monte Carlo simulator and a Variational Monte Carlo simulator written for the course "Understanding electron correlation with computer simulation" held at SISSA in the A.A. 2016-2017 by Sandro Sorella and Federico Becca.

The PMC simulator considers a particle in a 1D lattice (OBC) with a toy potential, while the VMC simulator considers a 1D Heisenberg model (PBC). Take a look at the [presentation](./Presentation/presentation.pdf) for more details.


## How to compile the code

### Prerequisites
- `g++` (or `icpc`, not tested with `clang`)
- Armadillo, **>= 7.100** ([http://arma.sourceforge.net/](http://arma.sourceforge.net/)). It's important to meet the minimum version requirement, otherwise the program will not compile.
- `qmake`. Note that an installation of QT is not needed, you just need `qmake` (which doesn't depend on QT).

### Compilation
Grab the latest version of **PMC-VMC** from the GitHub repository:

	git clone https://github.com/matteosecli/PMC-VMC
	
Then, if you want e.g. to compile the PMC simulator, type

	cd PMC-VMC/ProjectionMC && qmake && make
	
If everything goes as it should, you now have an executable `ProjectionMC` in your folder. You can run it via

	./ProjectionMC
	
**Notes:**
For the VMC simulator, the `VariationalMC.pro` file also contains linking to LAPACK, BLAS and OpenMP. While not strictly needed here, they can be used to speedup Armadillo operations at a certain point. Since you have Armadillo you should also have those libraries (and therefore there should be no issue), but if you really have problems you can just link only Armadillo by commenting the proper lines as explained in the file. I tested it and it seems to work.

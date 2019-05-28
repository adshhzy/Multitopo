Multitopo - Zhiyang Huang (2019)

------------------------------------

This code implements the algorithm described in

Huang, Zhi Yang, et al. "Topology-controlled Reconstruction of Multi-labelled Domains from Cross-sections
"  published at SIGGRAPH 2017.

The primary function of this program is to provide the topological-control surface reconstruction in multi-labeled domain.

Currently, the code is only tested on Mac OS X 10.13 or above, it should work at other platforms with minor modification.


BUILDING
======================================================================================================


The code has three dependencies: 1)Eigen/Eigen3,  2)Gurobi,   3)Tetgen 1.4

1. http://eigen.tuxfamily.org/index.php?title=Main_Page

2. http://www.gurobi.com/

3. http://wias-berlin.de/software/tetgen/

After download/install the three dependencies (Tetgen has been included in the external folder, rebuild it if necessary), please go to folder src and modify the relevant path according to your installation in the CMakeList.txt, e.g.

Eigen include file path: SET(EIGEN_INCLUDE_DIRS "/opt/local/include/")

Gurobi include file and lib path: SET(GUROBI_INCLUDE_DIRS "/Library/gurobi811/mac64/include/") SET(GUROBI_LIB_DIRS "/Library/gurobi811/mac64/lib/")


Then build the Cmake file and make:
$cmake .
$make

In the src directory, there should be an executable called "multitopo" (or "multitopo.exe" on Windows if it is successfully built).


RUNNING
======================================================================================================

TBD

For further questions about the code and the paper, please contact Zhiyang Huang at adshhzy@gmail.com or zhiyang.huang@wustl.edu (might be invalid after he graduated). You can also contact Prof. Tao Ju at taoju@wustl.edu.




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

In the "src" directory, there should be an executable called "multitopo" (or "multitopo.exe" on Windows if it is successfully built).


RUNNING
======================================================================================================

To run the code from the command line, type:

$./multitopo -i input_file_name -o output_file_path -p protocol_file_name [-d ray_density] [-m propgation_method] [-s compute_save_load_mode] [-j]

where:
1. -i: followed by the path of the input file. input_file_name is a path to the input file. currently, support file format includes ".contour". See https://www.cse.wustl.edu/~taoju/lliu/paper/ctr2suf/program.html for the description of the format. Note that in order to clearly specifying the topological constraints, we strongly suggest you use 0 as background, and 1, 2, ... n for the rest labels to describe your n+1 labeled domain (n material + 1 background). Please avoid using arbitrary numbers as lable for specifying the materials.

2. -o: followed by the path of the output path. output_file_path is a path to the folder for generating output files. In the output folder, three subfolder will be created to store the output data. Including the "cross"

3. -p: followed by the path of the protocol file for specifying the topological constriain. protocol_file_name is a path to the protocol file. The protocol file should be a .txt file. In the first row, it should be the constraint on the number of component of each label, orderly. If you don't want to put constraint on some of the labels, put -1 in the position. Note that the first number should be the constraint of numbers of components on the 0 lable, which usually be the background. Starting from the second row, you should specifying the genus constraint on each component of the label, if you don't want to specify the genus constraint, put -1 on that position. Note that the numbers on each row should match the number of components specified in the first row for the label. If you don't specify the numbers of component in the first row for that label, you are not allowed to specify the genus constraint on that label and you have to put -1 in the correspond row.

-1 -1 1 2 -1 3
-1
-1
0
0 0
-1
1 2 -1


3. -s: optional argument. Followed by a unsigned integer number indicating the number of voxels in each dimension for the implicit surfacing. Only If -s is included in the command line, the program would output the surface ([input file name]_surface.ply). We recomment using 100 for a default value, and you should set this according to your inputs and the precision of the output. Notices that the surfacing algorithm takes quite a long time for surfacing the zero-level set, and it depends on the resolution and the shape of the zero-level set.

4. -o: optional argument. followed by the path of the output path. output_file_path is a path to the folder for generating output files. Default the folder of the input file.


Some examples have been placed at data folder for testing:
1. $./vipss -i ../data/hand_ok/input.xyz -l 0 -s 200
2. $./vipss -i ../data/walrus/input.xyz -l 0.003 -s 100

The program will generate the predicted normal in [input file name]_normal.ply.
If -s is included in the command line, the program will generate the surface as the zero-level set of the solved implicit function ([input file name]_surface.ply).


For further questions about the code and the paper, please contact Zhiyang Huang at adshhzy@gmail.com or zhiyang.huang@wustl.edu (might be invalid after he graduated). You can also contact Prof. Tao Ju at taoju@wustl.edu.




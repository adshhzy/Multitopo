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

$./multitopo -i input_file_name -o output_file_path -p protocol_file_name [-d ray_density] [-m propgation_method] [-s compute_save_load_mode] [-t io_meta] [-j] [-k]

where:
1. -i: followed by the path of the input file. input_file_name is a path to the input file. currently, support file format includes ".contour". See https://www.cse.wustl.edu/~taoju/lliu/paper/ctr2suf/program.html for the description of the format. Note that in order to clearly specifying the topological constraints, You have to use 0 as background, and 1, 2, ... n for the rest labels to describe your n+1 labeled domain (n material + 1 background). Please do not use arbitrary numbers as label for specifying the materials.

2. -o: followed by the path of the output path. output_file_path is a path to the folder for generating output files. In the output folder, three subfolder will be created to store the output data. Including the "cross-section" folder for storing the cell complexity (space arrangement created by the planes of the input), "suf" folder for storing the candidate topologies with each cell, "meta" folder for all infomation needed in the dynamic programming stage.

3. -p: followed by the path of the protocol file for specifying the topological constriain. protocol_file_name is a path to the protocol file. The protocol file should be a .txt file. In the first row, it should be the constraint on the number of component of each label, orderly. If you don't want to put constraint on some of the labels, put -1 in the position. Note that the first number should be the constraint of numbers of components on the 0 lable, which usually be the background. Starting from the second row, you should then specify the genus constraint on each component of the label, if you don't want to specify the genus constraint, put -1 on that position. Note that the numbers on each row should match the number of components specified in the first row for the label. If you don't specify the numbers of component in the first row for that label, you are not allowed to specify the genus constraint on that label and you have to put -1 in the correspond row.

-1 2 1 2 -1 3

-1

-1 -1

1

0 0

-1

1 2 -1

In the above protocol example, you are working on a 6-label domain. The protocol specifies that the label 1, 2, 3, 5 should have 2, 1, 2, 3 components respetively, while no constraints are placed on the lable 0 and 4, as shown at the first row. In the second row, a -1 is there since no constraint is placed on label 0 so you have to put -1 here. In the third row the protocol put no genus constraint on the two components of the lable 1. In the forth row the protocol put a constraint that the component of the lable 2 should be genus 1. In the fifth row the protocol put a constraint that the two components of the lable 3 should all be genus 0. In the sixth row, a -1 is there since no constraint is placed on label 4 so you have to put -1 here. In the seventh row,  the protocol put a constraint that the two of the components of the lable 5 should be genus 1 and 2 respectively, while no genus constraint is placed on the last component.

4. -d: optional argument. Followed by a unsigned integer number indicating the density of the ray using in the ray shooting algorithm while exploring the topologies within each cell. Default 2 and minimum 2. Higher this value, more ray will be used in the algorithm and you are supposed to get more topology variances in each cell.

5. -m: optional argument. Followed by a 0 or 1 (default) to indicate the method for propagating the lables from faces to each vertex in the triangulated plane structures. When we triangulate the plane structure created by all cross-sections, we include those contours/curves in the triangulation so we can get the label of each triangles in the structure. Now we need to specify the label of each vertex in the domain by propagating the labels from the triangles. If it is 0, then the label of the vertex will be maximum label of its neighborhood triangles. If it is 1, then the label of the vertex will be minimum label of its neighborhood triangles.    

6. -s: optional argument. Followed by a 0, 1, 2, 3  (default 0) to indicate whether to save or load the topological variaties explored by the ray shooting algorithm. 0: only run the ray shooting algorithm and output the cadidates to the "suf" folder. 1: run the ray shooting algorithm and output the cadidates to the "suf" folder, also output the meta data to the folder specified by the -t optional parameter. 2: load the meta data from the folder specified by the -t optional parameter, and also output the cadidates to the "suf" folder. 3: only load the meta data from the folder specified by the -t optional parameter, this mode is used for trying the dynamic programming of different protocols.

7. -t: followed by the path of the meta data of the ray shooting algorithm. The program will save or load the meta data that contain the result of the ray shooting algorithm if you choose -s 1, 2 or 3. Please make sure you didn't change the input when you load your precomputed result of the ray shooting algorithm.

8. -j: optional argument. If -j is included, the program will enable the junction point minimization process described in the future work part in the paper.

9. -k: optional argument. If -k is included, the program will enter iteratively dynamic programming mode. When you modify the protocol file and ready to rerun the dynamic programming, input any char except e, then the program will rerun the dynamic programming based on your modified protocol file. If you want to exit the program, input e.

We provide few examples in the data folder for testing.
1. $./multitopo -i ../data/mug.contour -o ./data/mug/ -p ../protocol/mug.txt -m 0
2. $./multitopo -i ../data/liver.contour -o ./data/liver/ -p ../protocol/liver.txt -m 1
3. $./multitopo -i ../data/mousebrain.contour -o ./data/mousebrain/ -p ../protocol/mousebrain.txt -m 0

The program will create three subfolder in the output path to store the output data. Including the "cross-section" folder for storing the cell complexity (space arrangement created by the planes of the input), "suf" folder for storing the candidate topologies with each cell as well as the combined surface, "meta" folder for all infomation needed and computed in the dynamic programming stage. To get the final result, check suf/outcombine.suf and suf/outcombine.obj. The picked candidated of each cell is output at meta/pickTopo_vec.txt

To read the meta file for further utilization, please go to the "reader" folder and check the code (C++) to see how to interpret the meta files.

For further questions about the code and the paper, please contact Zhiyang Huang at adshhzy@gmail.com or zhiyang.huang@wustl.edu (might be invalid after he graduated). You can also contact Prof. Tao Ju at taoju@wustl.edu.




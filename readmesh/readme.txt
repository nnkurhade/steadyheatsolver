read.f90 is basically the main program.
------------------------------------------------
all other .f90 files are dependancies.
you have to run the solution in visual studio
tested to work on VS2022 intel fortran latest version.
-----------------------------------------
does not work on gfortran! no idea why
-----------------------------------------
all elements_%%.txt and nodes_%%.txt files correspond to some mesh, you will notice they are in pairs
so while giving file names as input from these files (mesh input) make sure the %% part matches in both files

mesh files used in test cases of the project are

nodes_ca_0.1.txt -> test case 1,2,3 and 4
elements_ca_0.1.txt -> test case 1,2,3 and 4

nodes_c.txt -> test case 5
elements_c.txt -> test case 5

nodes_crazy.txt -> case 6
elements_crazy.txt -> case 6

you have to assign these with given boundary conditions and gma values while setting up the environment as instructed in the readme file of the parent directory
-----------------------------------------------------------------------------------------------------------------------------------------------------------------
phi.txt are all the temperature values along with centroids. the data structure is
--------------------------------------------
x-centroid | y-centroid | Temperature
--------------------------------------------
line.txt are all the temperature values along the mid-y horizontal line. the data structure is
--------------------------------------------
x-mid-y | y-mid-y | Temperature-mid-y
--------------------------------------------
one can process these files using the matlab and R post processing codes provided in the report folder of the parent dictionary
--------------------------------------------------------------------------------------------------------------------------------
some outputs are pre-recorded for convenience
phi_%% and line_%% correspond to the condition they are solved for.
detailed problem is provided in the report
dd-dirichlet inner dirichlet outer BC
nn-neumann inner  neumann outer  BC
dn-dirichlet inner  neumann outer BC
nd-naumann inner dirichlet outer BC
source-source term with neumann BC
crazy - the complex geometry case



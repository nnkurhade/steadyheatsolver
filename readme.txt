prerequisites: the code has been developed on microsoft Visual Studio 2022 with intel fortran environment. code works definitely here, not sure about other platforms.

VS 2022: https://visualstudio.microsoft.com/vs/
ifort compiler: https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.yfpdrg

----------------------------------------------------
complete solution is inside readmesh folder
it has another readme file. Kindly go through it to develope the thorough understanding of the code
-----------------------------------------------------------
objective: to solve steady state heat diffusion equation over unstructured grids. 
disclaimer: does not solve for mixed boundary conditions, only triangular elements allowed
---------------------------------------------------------------------------------------------
only been tested for gmsh's .msh output file, 
if you want to have any other format, reading that mesh into the required data structures if solely your responsibility
on the bright side, we only need:
--------------------------------------------
node id | nodeXcoordinate | nodeYcoordinate
--------------------------------------------
and,
-----------------------------------
element id | node1 | node2 | node3
-----------------------------------
as an output from our mesh reading codeblock, so this should not be too difficult to manage.
once we have this info at our hand, our code prepares some data structures that are required for the solver, this can take some time for a complicated mesh.
--------------------------------------------------------------------------------------------------------------------------------------------------------------
once you have the .msh file from gmsh, you cannot directly run the code with that.
the code needs two txt files derived from the .msh output file of gmsh
sample gmsh output file:
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
4 4 1 0
1 0 0 0 0 
2 1 0 0 0 
3 1 1 0 0 
4 0 1 0 0 
1 0 0 0 1 0 0 0 2 1 -2 
2 1 0 0 1 1 0 0 2 2 -3 
3 0 1 0 1 1 0 0 2 3 -4 
4 0 0 0 0 1 0 0 2 4 -1 
1 0 0 0 1 1 0 0 4 3 4 1 2 
$EndEntities
$Nodes
9 5 1 5
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
1 1 0 0
1 2 0 0
1 3 0 0
1 4 0 0
2 1 0 1
5
0.5 0.5 0
$EndNodes
$Elements
9 12 1 12
0 1 15 1
1 1 
0 2 15 1
2 2 
0 3 15 1
3 3 
0 4 15 1
4 4 
1 1 1 1
5 1 2 
1 2 1 1
6 2 3 
1 3 1 1
7 3 4 
1 4 1 1
8 4 1 
2 1 2 4
9 1 2 5 
10 4 1 5 
11 2 3 5 
12 3 4 5 
$EndElements
you want two files out of this .msh file, nodes.txt and elements.txt
to prepare these two files simply follow the following steps
---------
nodes.txt
---------
-open the .msh file using any text reader
-copy everything that is enclosed between $Nodes and $EndNodes
-for ex. from the file above you want to copy,
9 5 1 5
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
1 1 0 0
1 2 0 0
1 3 0 0
1 4 0 0
2 1 0 1
5
0.5 0.5 0
-save this copied data as nodes.txt in the same directory that the program is present, or you can specify the path at readtime.
-------------
elements.txt
-------------
-open the .msh file using any text reader
-copy everything that is starts with something that looks like 2 1 2 4 (which is followed by lines having 4 entries that represent element ID and vertices of the triangle)
and ends just before $EndElements
-do not copy all of the information available between $Elements and $EndElements, the program is not designed to handle that, just copy from the data format mentioned above
example
2 1 2 4
9 1 2 5 
10 4 1 5 
11 2 3 5 
12 3 4 5 
-you can go through gmsh documentation for the speficics, but in short these lines represent the element vertex connectivity
-save the copied data as elements.txt in the same directory as the solution is present
if you have followed everything uptill now you now have two files that include mesh information that the code can read.
you are always welcome to modify this codeblock and develop your own reading programs.
line 16:39 (codeblock to read nodes.txt) 
line 114:125 (codeblock to read elements.txt)
-----------------------------------------------
line 47:107, modified jarvis march algorithm, although slow it works, if you want to change the algorithm that identifies the convex hull, make those changes here.
--------------------------------------------------------------------------------------------------------------------------------------------------------------
Note:- when I say can be defined in any way user wants, I mean that in case a quantity has spatial dependancy, there are data structures in code that specify
some important positions that can be used to define that 'quantity' value at that point using that information. On the other hand if a quantity depends on the
'scalar' one is solving for, it can be placed inside the source term loop so that it evolves with the solution, given that evolution is physical
--------------------------------------------------------------------------------------------------------------------------------------------------------------
user inputs
--------
line 20
--------
filename of desired nodes.txt file corresponding to the geometry of user's interest
------------------------------------------------------------------------------------
line 115
--------
filename of desired elements.txt file corresponding to the geometry of user's interest
---------------------------------------------------------------------------------------
line 574
----------
temperature initialization (initial guess) for the domain cell centroids. can be defined in any way user wants.
------------------------------------------------------------------------------------------------
line 594
----------
diffusion coefficient 'gma' has been defined at cell faces, 
you can have a constant 'gma' or any variation of it as long as you have a 'gma' defined at every face.
the code only requires to know what is the value of 'gma' at the faces.
for the sake of simplicity as observed in most cases, 'gma' has been kept independant of 'scalar' we are solving for, (temperature in our case)
in case 'gma' is 'scalar' dependant, you can simply move the code blocks defining 'gma' and 'gma' dependent variables inside the source term loop and it will work.
(take care of the allocation)
------------------------------
line 650
---------
basically initialization of all values of the scalar at face centroids
initialization of the inner bounadary, 
in case of dirichlet solver, this is the inner boundary condition itself,
since it is basically a 'scalar' defined at faces, can be set up in any way user desires.
-----------------------------------------------------------------------------------------
line 654
---------
initialization of the outer boundary,
in case of dirichlet solver, this is the outer boundary condition itself,
since it is basically a 'scalar' defined at faces, can be set up in any way user desires.
-----------------------------------------------------------------------------------------
line 661
---------
neuman solver switch,
assign ".true." for presence of neuman boundary conditions.
assign ".false." if neuman boundary conditions are not present.
---------------------------------------------------------------
line 665
---------
heat flux seen by inner boundary faces, 
it is a 'scalar' defined at faces and can be initialized in any way user wants
positive value indicates flux entering the cell
------------------------------------------------
line 669
---------
heat flux seen by the outer boundary faces,
it is a 'scalar' defined at faces and can be initialized in any way user wants
positive value indicates flux entering the cell
note: in case inner boundary is absent, it does not matter what heat flux value is assigned at inner boundary faces, only outer boundary will be solved for.
------------------------------------------------------------------------------------------------------------------------------------------------------------
line 673
---------
x-coordinate of fixed value point needed to fix the all neuman boundary solution
must be inside the domain of the mesh
--------------------------------------
line 674
---------
y-coordinate of fixed value point needed to fix the all neuman boundary solution
must be inside the domain of the mesh
--------------------------------------
line 683
---------
fixed 'scalar'(or temperature in our case) value that one wants the fixed point to have
----------------------------------------------------------------------------------------
line 690
---------
switch for the situation where all inner boundary faces have dirichlet boundary conditions and all outer boundary faces have neuman boundary conditions
if this is what one is solving for, assign ".true."
otherwise should be assigned ".false."
strictly follow the instructions mentioned in the comments from line 686 to line 689
------------------------------------------------------------------------------------
line 696
---------
switch for the situation where all outer boundary faces have dirichlet boundary conditions and all inner boundary faces have neuman boundary conditions
if this is what one is solving for, assign ".true."
otherwise should be assigned ".false."
strictly follow the instructions mentioned in the comments from line 692 to line 695
------------------------------------------------------------------------------------
line 706
--------
relaxation parameter. can be used to make the code converge, or for faster convergence (underrelaxation and overrelaxation)
use with caution
in case you do not want any relaxation, default value should be assigned 1.0
feature is only enabled for neuman boundary conditions
-----------------------------------------------------------------------------
line 708
--------
if (1) then
if you want the code to go through the solver
if (0) then
skips the solver
----------------
line 721
---------
constant source term. defined at every element centroid, can be defined in any way user wants
----------------------------------------------------------------------------------------------
line 722
---------
phi dependent source term. defined at every element centroid, can be defined in any way user wants,
but should have at least some dependancy on phi
evolves with the solution
--------------------------
line 725
---------
if (1) then
solver uses iterative green gauss cell based gradient calculation method, which in my implementation was found to be very unstable and difficult to converge
and hence not recommended at all
if someone wants to debug, the codeblocks are there, you are welcome!
technically according to the literature, this algorithm should work
if (0) then
the solver skips this method
-----------------------------
line 785
---------
if (1) then
solver uses simplified green gauss cell based gradient calculation method, 
this method was found to work very well in case of coarse mesh, increasing skewness in the mesh decreases the accuracy of the solution
on the other hand I know a test case where this method consistantly gives inaccurate solution.
hence although this works well for some cases it is not a recommended method to calculate gradients
if (0) then
the solver skips this method
----------------------------
line 803
--------
if (1) then
the solver uses least squares method to calculate cell gradients
this method can be computationally expensive but since we are dealing with only 3 faces here it is actually not.
this method is very accurate and skewness independent hence recommended to always use this method.
if (0) then
the solver skips this method
-----------------------------
line 931
--------
inner loop tolerance
set 1e-6 value to your desired tolerance
---------------------
line 935
--------
outer loop tolerance
set 1e-6 value to your desired tolerance
--------------------
line 958
--------
user can specify filename for the file where phi values are written along with the coordinates of the centroids where they are found
-------------------------------------------------------------------------------------------------------------------------------------
line 983
--------
user can specify filename for the file where the mid-y line is extracted
it basically does the same job as phi file does but the values are extracted on a line instead of extracting it from the whole domain
-------------------------------------------------------------------------------------------------------------------------------------
print counter controls are present at line 881 (inner), 717(outer) and 748(grad), ideally you want inner counter off, as it slows the code down
-------------------------------------------------------------------------------------------------------------------------------------
for additional information kindly refer the comments made in the code
----------------------------------------------------------------------
although the code is bit long, it is modular and hence basically a build up of different code blocks, 
the code thus can be broken down into functions and made short
but given that such a function is called mostly once,
 I did not see the point of making it that way which will indeed increase the complexity of the code and result in too many dependancies at build time
-------------------------------------------------------------------------------------------------------------------------------------------------------
some unavoidable dependancies are still present. They are included in the project and should be compiled together to make the code work.
-----------------------------------------------------------------------------------------------------------------------------------------
debugging blocks are available at the end of each code block, users can make use of it easily to have a better understanding of data structures and,
in case something goes wrong; for debugging purposes.
-----------------------------------------------------
this solver, given a physically possible problem, can solve steady scalar diffusion equation for any geometry, with or without source term.
outer and inner boundary are thus merely features of this code that makes assigning boundary values easier for the specific set of problems we are targeting.
--------------------------------------------------------------------------------------------------------------------------------------------------------------
visualization
-------------
two files me6151projectplotting.m (for matlab) and visualization_rstudio (for R) are available in the project folder. You have to modify their paths of "phi" and 
"line" reading functions as per your system. Then once you are through your simulation on the code, these two codes plot the outputs "phi" and "line" for better
understanding.








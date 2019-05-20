# Juliorus - An isometric immersion of the flat torus in Euclidean space
This program is used to output either .wrl or .xml (yafaray-xml format) of an isometric immersion of the flat torus in Euclidean 3-space. These formats can be used by a raytracer to generate the following image:

![alt text](https://i.imgur.com/Jfdf24f.jpg)
**NB:** The following torus was calculated on a grid size 8500 and with a pixel resolution of 9000x9000. The default settings are much lower, but can easily be changed.

## Dependencies
Before running this program, you need to install the following packages in Julia:
```
LinearAlgebra
ProgressMeter
StaticArrays
DifferentialEquations
Plots
SpecialFunctions
```
One installs these programs either manually, or by running
InstallPackages.jl in the src folder. To do this, one opens a terminal in the src folder and runs the command
```
julia InstallPackages.jl
```
The program is run by running the command in the terminal
```
julia -O3 Main.jl
```
in the src folder.

By default, the program does three convex integrations and outputs
3 .wrl files. These can be opened with the free program paraview, or another renderer that can view .wrl. Some windows PC's actually have an application pre-installed to open .wrl format files!

At default, the grid size is very low, but the first and second corrugations are in okay quality. To see what happens in the third quality, a higher grid resolution is necessary.

Alternatively, one can set `save_xml = true` in `Main.jl`, to output 3 .xml files. One generates a raytraced a png by installing Yafaray-v3 running
```
yafaray-xml -f png {your .xml file here}
```

## Reading the code
In the comments, we refer to the following source for the algorithm

[B] Borrelli, Vincent; Jabrane, Saïd; Lazarus, Francis; Thibert, Boris Isometric embeddings of the square flat torus in ambient space. Ensaios Matemáticos [Mathematical Surveys], 24. Sociedade Brasileira de Matemática, Rio de Janeiro, 2013. ii+91 pp. ISBN: 978-85-8337-011-6 (Reviewer: Ivan Izmestiev) 53C42 (53A05) 

This source has a chapter on the implementation. We follow this to generate the points on the torus.

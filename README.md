# curve-my-triangle
CGDI project of Corto Cristofoli and Max Royer based on the article *Curved PN 
Triangles* by Alex VLACHOS, Jörg PETERS, Chas BOYD and Jason L. MITCHELL.

Voir le [git](https://github.com/corto-cristofoli/curve-my-triangle)

## Goal
The idea of the project is to take a 3d *lowpoly* mesh and to apply a technique 
called **curved PN triangles tesselation** to smooth the mesh and make it more 
organic.

## How to compile
- The project uses **Geometry Central** and **Polyscope** libraries. Running 
    the commands below should install all the dependencies.
- `cmake -B build` in the root folder
- `make` in the build folder

## How to run the code
- `$alpha` default value is -0.5
- for `$lod` you should start by trying 2
- `./build/src/main -s $source_file_path -o $output_file_path -a $alpha -l $lod`

## Structure
The report is in the `rapport` folder, and we also provided some `obj` file in 
the `data` folder. The article we based our work on is inside of `docs`, under 
the name `VlachosPeters.pdf`.


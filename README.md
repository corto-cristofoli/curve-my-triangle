# curve-my-triangle
CGDI project of Corto Cristofoli and Max Royer based on the article *Curved PN 
Triangles* by Alex VLACHOS, Jörg PETERS, Chas BOYD and Jason L. MITCHELL.

## How to compile
- `cmake -B build` in the root folder
- `make` in the build folder

## How to run the code
- `./build/src/main -s $source_file_path -o $output_file_path -a $alpha -l $lod`

The report is in the `rapport` folder, and we also provided some `obj` file in the `data` folder.

## Goal
- The idea of the project is to take a 3d *lowpoly* mesh and to apply a technique called **curved PN triangles
tesselation** to smooth the mesh and make it more organic.


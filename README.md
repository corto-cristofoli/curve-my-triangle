# curve-my-triangle
CGDI project with Max Royer based on the article *Curved PN Triangles* by Alex VLACHOS, Jörg PETERS, Chas BOYD
and Jason L. MITCHELL.

## Goal
The idea of the project is to take a 3d *lowpoly* mesh and to apply a technique called **curved PN triangles
tesselation** to smooth the mesh and make it more organic.

## Some notes
- "Normal component of **PN triangle** is independently specified from geometric component".
- $b_{ijk}$: geometric point in the triangle. $i+j+k=3$.
    - $b_{300}, b_{030}, b_{003}$: vertex coefficients
    - $b_{210}, b_{120}, b_{021}, b_{012}, b_{201}, b_{102}$: tangent coefficients
    - $b_{111}$: center coefficient
- `lod`: Level of detail, i. e. number of control point inside an edge.

## Todos


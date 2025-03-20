# curve-my-triangle
CGDI project with Max Royer based on the article *Curved PN Triangles* by Alex VLACHOS, Jörg PETERS, Chas BOYD
and Jason L. MITCHELL.

## How to compile
- `cmake -B build` in the root folder
- `make` in the build folder

## Goal
- The idea of the project is to take a 3d *lowpoly* mesh and to apply a technique called **curved PN triangles
tesselation** to smooth the mesh and make it more organic.

## Some notes
- "Normal component of **PN triangle** is independently specified from geometric component".
- $b_{ijk}$: geometric point in the triangle. $i+j+k=3$.
    - $b_{300}, b_{030}, b_{003}$: vertex coefficients
    - $b_{210}, b_{120}, b_{021}, b_{012}, b_{201}, b_{102}$: tangent coefficients
    - $b_{111}$: center coefficient
- `lod`: Level of detail, i. e. number of control point inside an edge.

## Algorithm
0. Place the $b_{ijk}$ at pos $(iP_1 + jP_2 + kP_3)/3$.
1. Vertex coefficient stay at the same place so they match corner pos.
2. For the 2 closest tangent coef of a corner $l\in\{i,j,k\}$, project them into the plan defined by $N_l$.
3. Move central point from start pos $V$ to average of all 6 tangent point position, then continue for $1/2$
   of the movement.

## Todos
- Maybe we start with only geometric component and we do normal later ?
- What do we do if we want some sharp edges ?


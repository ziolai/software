# My Simple Geometry
Poisson Equation on unit interval / square / cube domain.

## Section 1: Unit Square Domain 
Here we assume $0 < x < 1$. 

### Domain with three (3) cells 
Generate mesh with 3 cells using blockMesh 

Check mesh using checkMesh: shows 3 hex elements using 

blocks
(
  hex (0 1 2 3 4 5 6 7) (3 1 1) edgeGrading (1 1 1 1 2 2 2 2 3 3 3 3)
);

Cell centered finite volume representation with resulting matrix 
$A = \begin{pmatrix} 1.5 & -1.5 & 0 \\ -1.5 & 3 & -1.5 \\ 0 & -1.5 & 1.5 \end{pmatrix}$

Resulting Compressed Sparse Row (CSR) representation 

$\text{rows} = \begin{pmatrix} 0 & 2 & 5 & 7 \end{pmatrix}$

$\text{cols} = \begin{pmatrix} 0 & 1 & 0 & 1 & 2 & 1 & 2 \end{pmatrix}$

$\text{vals} = \begin{pmatrix} 1.5 & -1.5 & -1.5 & 3 & -1.5 & -1.5 & 1.5 \end{pmatrix}$

### Domain with four (4) cells 
Generate mesh using blockMesh with 4 cells ($ncells = 4$) with 3 internal faces ($nfaces = 3$)

Check mesh using checkMesh: shows 4 hex elements using 

blocks
(
  hex (0 1 2 3 4 5 6 7) (4 1 1) edgeGrading (1 1 1 1 2 2 2 2 3 3 3 3)
);

Cell centered finite volume representation with resulting matrix of size 
$ncells \times ncells = 4 \times 4$ and $nfaces$ non-zero elements in the upper and lower triangle:  
$A = \begin{pmatrix} 2 & -2 & 0  & 0 \\ 
                    -2 &  4 & -2 & 0 \\ 
                     0 & -2 & 4  & -2 \\
                     0 &  0 & -2 & 2 \end{pmatrix}$

Resulting Compressed Sparse Row (CSR) representation 

$\text{rows} = \begin{pmatrix} 0 & 2 & 5 & 8 & 10 \end{pmatrix}$

$\text{cols} = \begin{pmatrix} 0 & 1 & 0 & 1 & 2 & 1 & 2 & 3 & 2 & 3 \end{pmatrix}$

$\text{vals} = \begin{pmatrix} 2 & -2 & -2 & 4 & -2 & -2 & 4 & -2 & -2 & 2 \end{pmatrix}$

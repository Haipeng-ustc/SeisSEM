# SeisFEM

This project is for our class project for Computational Seismology at USTC.

The team members include: Haipeng Li, Hongbo Han, and Wancong Zhao.

I will continue to update this repository for more features later.

Contact: haipengl@mail.ustc.edu.cn

## How to use

#### For regular mesh,
```bash
run: Case1_SEM_inner_mesh.m
```
#### For irregular mesh:

```bash
run: Case2_SEM_use_Gmsh_mesh.m	

Before that, you need to:

 1. Use the topography data (x and z) to generate mymesh.geo file by running Write_Gmsh.m;
 2. Load the mymesh.geo file into Gmsh and click: Mesh->2D. Then save the mesh by clicking: File->Save Mesh
 3. Run Case2_SEM_use_Gmsh_mesh.m 

Please notice that the order of the mesh generated by Gmsh should be the same with that in the main code. Here, P=4 is used.
```

## References
```bash
1. Pozrikidis, C. (2005). Finite and spectral element methods using Matlab. University of California San Diego, USA.
2. Caltech Lecture notes: GE 263 – Computational Geophysics The spectral element method by Jean-Paul Ampuero
3. Specfem2D from: https://github.com/geodynamics/specfem2d
```

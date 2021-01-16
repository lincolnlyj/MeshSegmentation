# MeshSegmentation
- Divide the mesh objects automatically into parts

- There's two version: 01-means segmentation and k-means segmentation

- The input file must be an obj file, and the amount of Faces is better to be under 4k. Since the algorism isn't that perfect, if there's too many faces, the program can be very very slow, and it can consume quite an amount of memory. (The release version can be way faster than the debug version)

- Reference: *Hierarchical Mesh Decomposition using Fuzzy Clustering and Cuts*

&nbsp;

**caution:** 
- *It's just for fun, there's plenty bugs in it, and the code style can be very ugly*
- *Please use visual studio to compile it, since the file is encoded as GB2312*
# MeshSegmentation
- Divide the mesh objects automatically into parts

- There's two version: 01-means segmentation and k-means segmentation, k-means segmentation can be found on the branch "KMeans"

- The input file must be an obj file, and the amount of Faces is better to be under 4k. Since the algorism isn't that perfect, if there's too many faces, the program can be very very slow.

- Reference: *Hierarchical Mesh Decomposition using Fuzzy Clustering and Cuts*

**caution:** *It's just for fun, there's plenty bugs in it, and the code style can be very ugly*
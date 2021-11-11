#ifndef MESH_H
#define MESH_H

//================================================================================
// TYPES for MESH
//================================================================================

typedef struct MeshGroupOfVertices   MeshGroupOfVertices;
typedef struct MeshGroupOfElements   MeshGroupOfElements;
typedef struct Mesh                  Mesh;

struct MeshGroupOfVertices
{
  int     num;
  int*    tag;
  double* x;
  double* y;
  double* z;
};

struct MeshGroupOfElements
{
  int   num;
  int*  tag;
  int** vertices;  // Nelem x NvertPerElem
};

struct Mesh
{
  MeshGroupOfVertices* vertices;
  MeshGroupOfElements* elements;
};

//================================================================================
// FUNCTIONS for MESH
//================================================================================

// Read the mesh from a gmsh-file (.msh) and store in a mesh-structure 'm'
Mesh* readMsh(char nameFile[128]);

// Save a solution 'u' in in mesh file (.msh)
void saveToMsh(double *val, Mesh *mesh, char name[], char fileName[]);

// Free space
void destroyMesh(Mesh *M);

#endif // MESH_H

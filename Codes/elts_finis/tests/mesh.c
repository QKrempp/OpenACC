#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"

//================================================================================
// Read the mesh from a gmsh-file (.msh) -- Mesh format 2.2
//================================================================================

Mesh* readMsh(char nameFile[128]){
  
  // Open file
  FILE *file = NULL;
  file = fopen(nameFile, "r");
  
  // Read mesh format
  float version;
  char line[128] = "";
  fgets(line, 128, file);
  if(strstr(line, "$MeshFormat") != NULL){
    fscanf(file, "%f\n", &version);
   // printf("   -> mesh format: %.1f\n", version);
  }
  else{
   // printf("ERROR: Impossible to read the mesh format.\n");
    exit(1);
  }
  if(version != 2.2f){
   // printf("ERROR: Mesh reader not available for mesh format %.1f.\n", version);
    exit(1);
  }
  
  // Read mesh
  Mesh *mesh = malloc(sizeof(Mesh));
  while(fgets(line, 128, file) != NULL){
  
    // scan VERTICES in gmsh file
    if(strstr(line, "$Nodes") != NULL){
      mesh->vertices = malloc(sizeof(MeshGroupOfVertices));
      
      // read number of vertices
      int allVertices;
      fscanf(file, "%i\n", &allVertices);
      mesh->vertices->num = allVertices;
      
      // read vertices
      mesh->vertices->tag = malloc(allVertices*sizeof(int));
      mesh->vertices->x = malloc(allVertices*sizeof(double));
      mesh->vertices->y = malloc(allVertices*sizeof(double));
      mesh->vertices->z = malloc(allVertices*sizeof(double));
      for(int v=0; v<allVertices; ++v){
        fscanf(file, "%i %lf %lf %lf",
          &mesh->vertices->tag[v],
          &mesh->vertices->x[v],
          &mesh->vertices->y[v],
          &mesh->vertices->z[v]);
      }
      
     // printf("   -> %i nodes\n", mesh->vertices->num);
    }
    
    // scan ELEMENTS in gmsh file
    if(strstr(line, "$Elements") != NULL){
      mesh->elements = malloc(sizeof(MeshGroupOfElements));
      
      // read number of elements
      int allTriangles;
      fscanf(file, "%i\n", &allTriangles);
      mesh->elements->num = allTriangles;
      
      // read elements
      mesh->elements->tag      = malloc(allTriangles * sizeof(int));
      mesh->elements->vertices = malloc(allTriangles * sizeof(int*));
      for(int e=0; e<allTriangles; e++) {
        int eGmshNum, eGmshType, eInfos, eTagPhi, eTagGeo, dummy;
        fscanf(file, "%i %i %i %i %i", &eGmshNum, &eGmshType, &eInfos, &eTagPhi, &eTagGeo);
        if(eGmshType != 2){  // 2 is for Triangles
         // printf("ERROR when reading .msh: mesh reader only for TRIANGULAR elements!\n");
        exit(1);
        }
        switch(eInfos){
          case 2: break;
          case 4: fscanf(file, "%i %i", &dummy, &dummy); break;
          default:// printf("ERROR when reading .msh: eInfos too large!\n");
				  exit(1); 
					break;
        }
        mesh->elements->tag[e] = eGmshNum;
        mesh->elements->vertices[e] = malloc(3 * sizeof(int*));
        for(int v=0; v<3; v++){
          fscanf(file, "%i", &dummy);
          mesh->elements->vertices[e][v] = dummy-1; // (gmsh is 1-index, here is 0-index)
        }
      }
      
     // printf("   -> %i elements\n", mesh->elements->num);
    }
  }
  
  fclose(file);
 
  return mesh;
}

//================================================================================
// Save a solution 'u' in a mesh file (.msh)
//================================================================================

void saveToMsh(double *val, Mesh *mesh, char name[], char fileName[])
{
  FILE* file = NULL;
  file = fopen(fileName, "w");
  fprintf(file, "$MeshFormat\n");
  fprintf(file, "2.1 0 8\n");
  fprintf(file, "$EndMeshFormat\n");
  
  fprintf(file, "$ElementNodeData\n");
  fprintf(file, "%i\n", 2);
  fprintf(file, "\"%s\"\n", name); // name of the view
  fprintf(file, "%i\n", 0);
  fprintf(file, "%i\n", 1);
  fprintf(file, "%i\n", 0); // ("Time")
  fprintf(file, "%i\n", 3);
  fprintf(file, "%i\n", 0); // ("timeStep")
  fprintf(file, "%i\n", 1); // ("numComp")
  fprintf(file, "%i\n", mesh->elements->num);  // total number of elementNodeData in this file
  
  for(int iTri=0; iTri<mesh->elements->num; iTri++){
    fprintf(file, "%i %i", mesh->elements->tag[iTri], 3);
    for(int v=0; v<3; v++){
      int i = mesh->elements->vertices[iTri][v];
      fprintf(file, " %g", val[i]);
    }
    fprintf(file, "\n");
  }
  fprintf(file, "$EndElementNodeData\n");
  fclose(file);
  
  return;
}

//================================================================================
// Free space
//================================================================================

void destroyMesh(Mesh *M){
  free(M->vertices->tag);
  free(M->vertices->x);
  free(M->vertices->y);
  free(M->vertices->z);
  free(M->vertices);
  
  free(M->elements->tag);
  for(int v=0; v<M->elements->num; v++){
    free(M->elements->vertices[v]);
  }
  free(M->elements->vertices);
  free(M->elements);
  
  free(M);
};

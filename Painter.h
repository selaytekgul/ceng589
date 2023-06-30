#define HAVE_INT8_T
#include "Mesh.h"
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoCone.h>



class Painter
{
public:
	SoSeparator* getShapeSep(Mesh* mesh);
	void triangleNormalVectors(Mesh* mesh);
	float calculateLength(int v[]);
	void normalize(int v[], float nor[]);
	void Painter::crossProductFunction(int v_A[], int v_B[], int c_P[]);
	float Painter::rayIntersectsTriangle(float* p, float* d, float* v0, float* v1, float* v2);
};

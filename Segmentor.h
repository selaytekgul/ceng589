#pragma once

#define HAVE_INT8_T
#include "Mesh.h"
//#include <Inventor/nodes/SoMaterial.h>
//#include <Inventor/nodes/SoCoordinate3.h>
//#include <Inventor/nodes/SoIndexedFaceSet.h>
//#include <Inventor/nodes/SoSeparator.h>
//#include <Inventor/nodes/SoShapeHints.h>
//#include <Inventor/nodes/SoCone.h>
#include <vector>
#include <array>
#include <cmath>



class Segmentor
{
public:
	//SoSeparator* getShapeSep(Mesh* mesh);
	void assignLengthValuesOfVertices(Mesh* mesh);
private:
	inline void vector(float v_A[], const float v_B[], const float v_C[]);
	inline void crossProduct(const float v_A[], const float v_B[], float CP[]);
	inline float calculateLengthOfVector(const float v[]);
	inline float innerProduct(const float v[], const float q[]);
	float rayIntersectsTriangle(float* p, float* d, float* v0, float* v1, float* v2);
	//void normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr);
};

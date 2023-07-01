#define HAVE_INT8_T
#include "Mesh.h"
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoCone.h>
#include <vector>
#include <array>



class Painter
{
public:
	SoSeparator* getShapeSep(Mesh* mesh);
	void triangleNormalVectors(Mesh* mesh);
private:
	float calculateLength(const float v[]);
	void normalize(const float v[], float nor[]);
	void crossProductFunction(const int v_A[], const int v_B[], int c_P[]);
	float rayIntersectsTriangle(float* p, float* d, float* v0, float* v1, float* v2);
	float normalizedValueInRange(double value, double min, double max);
	void normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr);
};

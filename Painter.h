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
	void crossProductFunction(const float v_A[], const float v_B[], float CP[]);
	float rayIntersectsTriangle(float* p, float* d, float* v0, float* v1, float* v2);
	void normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr);
};

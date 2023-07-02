#pragma once

#define HAVE_INT8_T
#include "Mesh.h"

class Segmentor
{
public:
	Segmentor(Mesh* mesh);
	~Segmentor();
	void assignLengthValuesOfVertices();
	void setColorValuesToVertices();
private:
	Mesh* mesh;
	void calculateCrossProduct(float crossProductVector[3], const triVertsCoords& coordinatesOfVerticesOfTriangle, const size_t selectedVertexNumber);
	void calculateShortestDiameter(const int triangleIndex, const int selectedVertexNumber, const triVertsIds vertexIdsOfTriangle, const float p[3], const float d[3]);
};

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
	void calculateShortestDiameter(int triangleIndex, int selectedVertexNumber, triVertsIds vertexIdsOfTriangle, float p[3], float d[3]);
};

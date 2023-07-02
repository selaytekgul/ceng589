#pragma once

#define HAVE_INT8_T
#include "Mesh.h"

class Segmentor
{
public:
	void assignLengthValuesOfVertices(Mesh* mesh);
	void setColorValuesToVertices(Mesh* mesh);
private:
	void Segmentor::selectATargetTriangle(Mesh* mesh, int triangleIndex, int selectedVertexNumber, triVertsIds vertexIdsOfTriangle, float p[3], float d[3]);
};

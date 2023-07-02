#pragma once

#define HAVE_INT8_T
#include "Mesh.h"

class Segmentor
{
public:
	void assignLengthValuesOfVertices(Mesh* mesh);
	void setColorValuesToVertices(Mesh* mesh);
private:
	void Segmentor::selectATargetTriangle(Mesh* mesh, int triangleIndex, int selectedVertexNumber, std::array<int, 3> vertexIdsOfTriangle, float p[3], float d[3]);
};

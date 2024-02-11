#pragma once

#include "Mesh.h"

class KMeans
{
public:
	KMeans(Mesh* mesh);
	~KMeans();
	void assignClusterIdsOfVertices();
	void setColorValuesToVertices();
private:
	Mesh* mesh;
	void calculateNormalVector(float crossProductVector[3], const triVertsCoords& coordinatesOfVerticesOfTriangle, const size_t selectedVertexNumber);
	void calculateShortestDiameter(const int triangleIndex, const int selectedVertexNumber, const triVertsIds& vertexIdsOfTriangle, const float p[3], const float d[3]);
};


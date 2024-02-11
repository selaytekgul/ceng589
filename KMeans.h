#pragma once

#include "Mesh.h"

struct Point {
	float x, y, z; // Assuming 3D points
};
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
	float distanceBetweenPoints3D(const float* point1, const float* point2);
	Vertex* findNearestPoint(const std::vector<Vertex*>& vertices, const Point& meanPoint);
	Point computeMeanPoint(const std::vector<Vertex*>& vertices);


};


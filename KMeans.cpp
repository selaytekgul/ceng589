#include "KMeans.h"
#include <random>

KMeans::KMeans(Mesh* mesh) 
	: mesh(mesh)
{}

KMeans::~KMeans() = default;

void KMeans::assignClusterIdsOfVertices()
{
	auto numberOfVerts = mesh->verts.size();


	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distrib(0, numberOfVerts);

	std::array<int, 6> centroids = {};
	for (int i = 0; i < 6; ++i) {
		centroids[i] = distrib(gen);
	}


}

void KMeans::setColorValuesToVertices()
{
	constexpr int numberOfColors = 6;
	constexpr float colors[numberOfColors][3] = {
		{0, 1, 0}, //green
		{0, 1, 1}, //cyan
		{0, 0, 1}, //blue
		{1, 0, 0}, //red
		{1, 0, 1}, //magenta
		{1, 1, 0} //yellow
	};

	int usedNumberOfColors = 6; //number of segmented groups

	//assign color values to the vertices according to their diameter values
	for (Vertex* vertex : mesh->verts)
	{
		// setting the color for a vertex
		int colorIndex = static_cast<int>(vertex->clusterId);

		//index out of bounds error check for colors array
		if (colorIndex > usedNumberOfColors - 1)
			colorIndex = usedNumberOfColors - 1;

		TD::fillWith(vertex->color, colors[colorIndex], 3);
	}
}



void KMeans::calculateNormalVector(float crossProductVector[3], const triVertsCoords& coordinatesOfVerticesOfTriangle, const size_t selectedVertexNumber)
{}

void KMeans::calculateShortestDiameter(const int triangleIndex, const int selectedVertexNumber, const triVertsIds& vertexIdsOfTriangle, const float p[3], const float d[3])
{}
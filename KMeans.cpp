#include "KMeans.h"
#include <cmath>
#include <random>
#include <limits>


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
	std::array<std::vector<Vertex*>, 6> vertexClusters = {};

	for (Vertex* vertex : mesh->verts)
	{
		auto smallest_distance = std::numeric_limits<int>::max();
		int nearest_cluster_id = 0;
		for (size_t i = 0; i < centroids.size(); i++)
		{
			auto distance = distanceBetweenPoints3D(mesh->verts.at(centroids[i])->coords, vertex->coords);

			//smallest_distance, nearest_cluster_id = smallest_distance > distance ? distance, i: smallest_distance , nearest_cluster_id;
			if (smallest_distance > distance) {
				smallest_distance = distance;
				nearest_cluster_id = i;
			}
		}
		vertexClusters[nearest_cluster_id].push_back(vertex);
		vertex->clusterId = nearest_cluster_id;
	}

	auto centroidsChanged = true;
	for (size_t epoch = 0; epoch < 1 && centroidsChanged; epoch++)
	{

		for (size_t i = 0; i < centroids.size(); i++)
		{
			Point meanPointCalculated = computeMeanPoint(vertexClusters[i]);
			Vertex* newCentroid = findNearestPoint(vertexClusters[i], meanPointCalculated);
			if (centroids[i] != newCentroid->idx)
				centroidsChanged = true;
			centroids[i] = newCentroid->idx;
		}

		for (Vertex* vertex : mesh->verts)
		{
			auto smallest_distance = std::numeric_limits<int>::max();
			int nearest_cluster_id = 0;
			for (size_t i = 0; i < centroids.size(); i++)
			{
				auto distance = distanceBetweenPoints3D(mesh->verts.at(centroids[i])->coords, vertex->coords);

				//smallest_distance, nearest_cluster_id = smallest_distance > distance ? distance, i: smallest_distance , nearest_cluster_id;
				if (smallest_distance > distance) {
					smallest_distance = distance;
					nearest_cluster_id = i;
				}
			}
			vertexClusters[nearest_cluster_id].push_back(vertex);
			vertex->clusterId = nearest_cluster_id;
		}
	}

}


// Function to compute the mean point from a list of vertices
Point KMeans::computeMeanPoint(const std::vector<Vertex*>& vertices) {
	Point mean = { 0.0f, 0.0f, 0.0f };
	int count = vertices.size();

	for (const auto& vertex : vertices) {
		mean.x += vertex->coords[0];
		mean.y += vertex->coords[1];
		mean.z += vertex->coords[2];
	}

	if (count > 0) {
		mean.x /= count;
		mean.y /= count;
		mean.z /= count;
	}

	return mean;
}

// Function to find the nearest point to the mean point from a vertex array
Vertex* KMeans::findNearestPoint(const std::vector<Vertex*>& vertices, const Point& meanPoint) {
	Vertex* nearestPoint = vertices.at(0);
	float minDistance = std::numeric_limits<float>::max();

	for (const auto& vertex : vertices) {
		float* coords = new float[3];
		coords[0] = meanPoint.x;
		coords[1] = meanPoint.y;
		coords[2] = meanPoint.z;
		float dist = distanceBetweenPoints3D(vertex->coords, coords);
		if (dist < minDistance) {
			minDistance = dist;
			nearestPoint = vertex;
		}
		delete[] coords;
	}

	return nearestPoint;
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


float KMeans::distanceBetweenPoints3D(const float* point1, const float* point2) {
	float deltaX = point2[0] - point1[0];
	float deltaY = point2[1] - point1[1];
	float deltaZ = point2[2] - point1[2];

	return fabs(sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ));
}

void KMeans::calculateNormalVector(float crossProductVector[3], const triVertsCoords& coordinatesOfVerticesOfTriangle, const size_t selectedVertexNumber)
{}

void KMeans::calculateShortestDiameter(const int triangleIndex, const int selectedVertexNumber, const triVertsIds& vertexIdsOfTriangle, const float p[3], const float d[3])
{}
#pragma once
#include <array>
#include "Mesh.h"

namespace TriangleMeshMath{
	inline std::array<int, 3> getVertexIdsOfTriangleAsStdArray(const Mesh* mesh, const int triangleIndex);
	inline std::array<std::array<float, 3>, 3> getCoordsOfOfTriangleAsStdArrayOfStdArray(const Mesh* mesh, std::array<int, 3> vertexIdsOfTriangle);
	inline std::array<std::array<float, 3>, 2> getOtherCoordsOfOfTriangleAsStdArrayOfStdArray(const std::array<std::array<float, 3>, 3> coordinatesOfVerticesOfTriangle, const int selectedVertexNumber);
	inline std::array<std::array<float, 3>, 2> getVectorsToTheOtherVertices(const std::array<std::array<float, 3>, 2> coordinatesOfOtherVertices, const int selectedVertexNumber, const std::array<std::array<float, 3>, 3> coordinatesOfVerticesOfTriangle);

	std::array<int, 3> getVertexIdsOfTriangleAsStdArray(const Mesh* mesh, const int triangleIndex)
	{
		std::array<int, 3> vertexIdsOfTriangle;
		vertexIdsOfTriangle[0] = mesh->tris[triangleIndex]->v1i;
		vertexIdsOfTriangle[1] = mesh->tris[triangleIndex]->v2i;
		vertexIdsOfTriangle[2] = mesh->tris[triangleIndex]->v3i;
		return vertexIdsOfTriangle;
	}
	
	std::array<std::array<float, 3>, 3> getCoordsOfOfTriangleAsStdArrayOfStdArray(const Mesh* mesh, const std::array<int, 3> vertexIdsOfTriangle)
	{
		std::array<std::array<float, 3>, 3> coordinatesOfVerticesOfTriangle;
		for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
		{
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				coordinatesOfVerticesOfTriangle[vertexNumber][coordinate] = mesh->verts[vertexIdsOfTriangle[vertexNumber]]->coords[coordinate];
			}
		}
		return coordinatesOfVerticesOfTriangle;
	}

	std::array<std::array<float, 3>, 2> getOtherCoordsOfOfTriangleAsStdArrayOfStdArray(const std::array<std::array<float, 3>, 3> coordinatesOfVerticesOfTriangle, const int selectedVertexNumber)
	{
		std::array<std::array<float, 3>, 2> coordinatesOfOtherVertices;
		int number = 0;
		for (size_t otherVertexNumber = 0; otherVertexNumber < 3; otherVertexNumber++)
		{
			if (selectedVertexNumber != otherVertexNumber)
			{
				//fill the coordinate values of other 2 vertices' array
				for (size_t coordinate = 0; coordinate < 3; coordinate++) {
					coordinatesOfOtherVertices[number][coordinate] = coordinatesOfVerticesOfTriangle[otherVertexNumber][coordinate];
				}
				number++;
			}
		}//find the other 2 vertices of the triangle
		return coordinatesOfOtherVertices;
	}

	std::array<std::array<float, 3>, 2> getVectorsToTheOtherVertices(const std::array<std::array<float, 3>, 2> coordinatesOfOtherVertices, const int selectedVertexNumber, const std::array<std::array<float, 3>, 3> coordinatesOfVerticesOfTriangle) {
		std::array<std::array<float, 3>, 2> vectorsToTheOtherVertices;
		for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
		{
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				vectorsToTheOtherVertices[otherVertexNumber][coordinate] =
					coordinatesOfOtherVertices[otherVertexNumber][coordinate]
					- coordinatesOfVerticesOfTriangle[selectedVertexNumber][coordinate];
			}
		}
		return vectorsToTheOtherVertices;
	}
}
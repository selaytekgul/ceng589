#pragma once
#include <array>
#include "Mesh.h"

namespace TriangleMeshMath{
	inline std::array<int, 3> getVertexIdsOfTriangleAsStdArray(const Mesh* mesh, const int triangleIndex);

	inline std::array<int, 3> getVertexIdsOfTriangleAsStdArray(const Mesh* mesh, const int triangleIndex)
	{
		std::array<int, 3> vertexIdsOfTriangle;
		vertexIdsOfTriangle[0] = mesh->tris[triangleIndex]->v1i;
		vertexIdsOfTriangle[1] = mesh->tris[triangleIndex]->v2i;
		vertexIdsOfTriangle[2] = mesh->tris[triangleIndex]->v3i;
		return vertexIdsOfTriangle;
	}
}
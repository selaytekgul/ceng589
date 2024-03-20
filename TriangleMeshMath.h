#pragma once
#include "Mesh.h"
#include "TypeDefinitions.h"

namespace TriangleMeshMath
{
	inline triVertsIds getVertexIdsOfTriangle(const Mesh* mesh, const int triangleIndex);
	inline edgeVertsIds getVertexIdsOfEdge(const Mesh* mesh, const int edgeIndex);
	inline triVertsCoords getCoordsOfTriangle(const Mesh* mesh, triVertsIds vertexIdsOfTriangle);
	inline triOtherVertsCoords getOtherCoordsOfTriangle(const triVertsCoords& coordsOfVerticesOfTriangle, const int selectedVertexNumber);
	inline triOtherVertsCoords getVectorsToTheOtherVertices(const triOtherVertsCoords& coordsOfOtherVertices, const int selectedVertexNumber, const triVertsCoords& coordsOfVerticesOfTriangle);
	inline int getOtherVertexIdOfTriangle(const Mesh* mesh, const int triangleIndex, const int edgeId);

	triVertsIds getVertexIdsOfTriangle(const Mesh* mesh, const int triangleIndex)
	{
		triVertsIds vertexIdsOfTriangle;
		vertexIdsOfTriangle[0] = mesh->tris[triangleIndex]->v1i;
		vertexIdsOfTriangle[1] = mesh->tris[triangleIndex]->v2i;
		vertexIdsOfTriangle[2] = mesh->tris[triangleIndex]->v3i;
		return vertexIdsOfTriangle;
	}
	
	edgeVertsIds getVertexIdsOfEdge(const Mesh* mesh, const int edgeIndex)
	{
		edgeVertsIds vertexIdsOfEdge;
		vertexIdsOfEdge[0] = mesh->edges[edgeIndex]->v1i;
		vertexIdsOfEdge[1] = mesh->edges[edgeIndex]->v2i;
		return vertexIdsOfEdge;
	}
	
	triVertsCoords getCoordsOfTriangle(const Mesh* mesh, const triVertsIds vertexIdsOfTriangle)
	{
		triVertsCoords coordinatesOfVerticesOfTriangle;
		for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
		{
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				coordinatesOfVerticesOfTriangle[vertexNumber][coordinate] = mesh->verts[vertexIdsOfTriangle[vertexNumber]]->coords[coordinate];
			}
		}
		return coordinatesOfVerticesOfTriangle;
	}

	triOtherVertsCoords getOtherCoordsOfTriangle(const triVertsCoords& coordsOfVerticesOfTriangle, const int selectedVertexNumber)
	{
		triOtherVertsCoords coordsOfOtherVertices;
		int number = 0;
		for (size_t otherVertexNumber = 0; otherVertexNumber < 3; otherVertexNumber++)
		{
			if (selectedVertexNumber == otherVertexNumber)
				continue;
			//fill the coordinate values of other 2 vertices' array
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				coordsOfOtherVertices[number][coordinate] = coordsOfVerticesOfTriangle[otherVertexNumber][coordinate];
			}
			number++;
		}
		return coordsOfOtherVertices;
	}

	triOtherVertsCoords getVectorsToTheOtherVertices(const triOtherVertsCoords& coordsOfOtherVertices, const int selectedVertexNumber, const triVertsCoords& coordsOfVerticesOfTriangle) {
		triOtherVertsCoords vectorsToTheOtherVertices;
		for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
		{
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				vectorsToTheOtherVertices[otherVertexNumber][coordinate] =
					coordsOfOtherVertices[otherVertexNumber][coordinate]
					- coordsOfVerticesOfTriangle[selectedVertexNumber][coordinate];
			}
		}
		return vectorsToTheOtherVertices;
	}

	int getOtherVertexIdOfTriangle(const Mesh* mesh, const int triangleIndex, const int edgeId)
	{
		triVertsIds vertIdxsInTri = getVertexIdsOfTriangle(mesh, triangleIndex);
		edgeVertsIds edgeIdxs = getVertexIdsOfEdge(mesh, edgeId);
		int retVertId = 0;
		for (size_t i = 0; i < 3; i++)
		{
			const int v = vertIdxsInTri[i];
			//auto it = std::find(vertIdxsInTri.begin(), vertIdxsInTri.end(), v);
			//if (it == vertIdxsInTri.end()) {
			//	retVertId = v;
			//}
			auto it = std::find(edgeIdxs.begin(), edgeIdxs.end(), v);
			if (it == edgeIdxs.end()) {
				retVertId = v;
			}
		}
		return retVertId;
	}
}
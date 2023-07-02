#pragma once
#include <array>
#include "Mesh.h"
#include "TypeDefinitions.h"

namespace TriangleMeshMath{
	inline triVertsIds getVertexIdsOfTriangleAsStdArray(const Mesh* mesh, const int triangleIndex);
	inline triVertsCoords getCoordsOfOfTriangleAsStdArrayOfStdArray(const Mesh* mesh, triVertsIds vertexIdsOfTriangle);
	inline triOtherVertsCoords getOtherCoordsOfOfTriangleAsStdArrayOfStdArray(const triVertsCoords coordinatesOfVerticesOfTriangle, const int selectedVertexNumber);
	inline triOtherVertsCoords getVectorsToTheOtherVertices(const triOtherVertsCoords coordinatesOfOtherVertices, const int selectedVertexNumber, const triVertsCoords coordinatesOfVerticesOfTriangle);

	triVertsIds getVertexIdsOfTriangleAsStdArray(const Mesh* mesh, const int triangleIndex)
	{
		triVertsIds vertexIdsOfTriangle;
		vertexIdsOfTriangle[0] = mesh->tris[triangleIndex]->v1i;
		vertexIdsOfTriangle[1] = mesh->tris[triangleIndex]->v2i;
		vertexIdsOfTriangle[2] = mesh->tris[triangleIndex]->v3i;
		return vertexIdsOfTriangle;
	}
	
	triVertsCoords getCoordsOfOfTriangleAsStdArrayOfStdArray(const Mesh* mesh, const triVertsIds vertexIdsOfTriangle)
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

	triOtherVertsCoords getOtherCoordsOfOfTriangleAsStdArrayOfStdArray(const triVertsCoords coordinatesOfVerticesOfTriangle, const int selectedVertexNumber)
	{
		triOtherVertsCoords coordinatesOfOtherVertices;
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

	triOtherVertsCoords getVectorsToTheOtherVertices(const triOtherVertsCoords coordinatesOfOtherVertices, const int selectedVertexNumber, const triVertsCoords coordinatesOfVerticesOfTriangle) {
		triOtherVertsCoords vectorsToTheOtherVertices;
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
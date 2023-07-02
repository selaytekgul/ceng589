#include "Segmentor.h"
#include "VectorMath.h"
#include "TriangleMeshMath.h"

void Segmentor::assignLengthValuesOfVertices(Mesh* mesh)
{
	//loop through the triangles to trace each vertex, find normals, draw rays, calculate and add lengths to vertex's length attribute
	for (size_t triangleIndex = 0; triangleIndex < mesh->tris.size(); triangleIndex++)
	{
		//get the vertex index number of the vertices of the triangle at hand
		triVertsIds vertexIdsOfTriangle = TriangleMeshMath::getVertexIdsOfTriangle(mesh, triangleIndex);

		//get the coordinates of the vertices of the triangle at hand (by using the vertex index numbers)
		triVertsCoords coordinatesOfVerticesOfTriangle = TriangleMeshMath::getCoordsOfTriangle(mesh, vertexIdsOfTriangle);
		
		//select a vertex from 3 vertices of triangle
		for (size_t selectedVertexNumber = 0; selectedVertexNumber < 3; selectedVertexNumber++)
		{
			//find the other 2 vertices of the triangle
			triOtherVertsCoords coordinatesOfOtherVertices = TriangleMeshMath::getOtherCoordsOfTriangle(coordinatesOfVerticesOfTriangle, selectedVertexNumber);

			//create two vectors from the selected vertex
			triOtherVertsCoords vectorsToTheOtherVertices = TriangleMeshMath::getVectorsToTheOtherVertices(coordinatesOfOtherVertices, selectedVertexNumber, coordinatesOfVerticesOfTriangle);

			//store two vectors from the selected vertices in a 2D array
			float vectorsToTheOtherVerticesArray[2][3];
			for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
			{
				TD::fillWith(vectorsToTheOtherVerticesArray[otherVertexNumber], vectorsToTheOtherVertices[otherVertexNumber], 3);
			}

			//calculate cross product of the two vectors
			float crossProductVector[3];
			VectorMath::crossProduct(vectorsToTheOtherVerticesArray[0], vectorsToTheOtherVerticesArray[1], crossProductVector);

			//p is the selected vertex of the base triangle
			float p[3];
			TD::fillWith(p, coordinatesOfVerticesOfTriangle[selectedVertexNumber], coordinatesOfVerticesOfTriangle[selectedVertexNumber].size());

			//d is the normal vector from the selected vertex of the base triangle drawn according to the other two vertices
			float d[3];
			TD::fillWith(d, crossProductVector, 3);

			calculateShortestDiameter(mesh, triangleIndex, selectedVertexNumber, vertexIdsOfTriangle, p, d);
		}//select a vertex from 3 vertices of triangle
	}//loop through the triangles to trace each vertex, find normals, draw rays, calculate and add lengths to vertex's length attribute
	mesh->setMinMaxLenghts();
	mesh->discardInfAndNegativeLenghts();
}

void Segmentor::calculateShortestDiameter(Mesh* mesh, int triangleIndex, int selectedVertexNumber, triVertsIds vertexIdsOfTriangle, float p[3], float d[3]) {
	Vertex* selectedVertex = mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]];
	//for each of the target triangles:
	for (size_t targetTriangleIndex = 0; targetTriangleIndex < mesh->tris.size(); targetTriangleIndex++)
	{
		//make sure that the target triangle is not the same as the base triangle
		if (triangleIndex == targetTriangleIndex)
			continue;

		triVertsIds vertexIdsOfTargetTriangle = TriangleMeshMath::getVertexIdsOfTriangle(mesh, targetTriangleIndex);

		//get the coordinates of the vertices of the target triangle at hand (by using the vertex index numbers)
		triVertsCoords coordinatesOfVerticesOfTargetTriangle = TriangleMeshMath::getCoordsOfTriangle(mesh, vertexIdsOfTargetTriangle);

		//fill the 2D raw array of targetTriangleVertices with the coordinate values
		float targetTriangleVertices[3][3];
		for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
			TD::fillWith(targetTriangleVertices[vertexNumber], coordinatesOfVerticesOfTargetTriangle[vertexNumber], 3);

		//calculate the diameters
		const float intersectionLength = VectorMath::rayTriangleIntersectLength(p, d, targetTriangleVertices[0], targetTriangleVertices[1], targetTriangleVertices[2]);
		const float previousLength = selectedVertex->diameter;
		if (intersectionLength > 0 && (previousLength <= 0 || intersectionLength < previousLength))
			selectedVertex->diameter = intersectionLength;
	}
}
void Segmentor::setColorValuesToVertices(Mesh* mesh) {
	//assign color values to the vertices according to their diameter values
	for (int i = 0; i < (int)mesh->verts.size(); i++)
	{
		float k = Vertex::maxDiameter / 8.0f;
		if (mesh->verts[i]->diameter < k)
		{
			mesh->verts[i]->color[0] = 0;
			mesh->verts[i]->color[1] = 1;
			mesh->verts[i]->color[2] = 0;
		}
		else if (mesh->verts[i]->diameter >= k && mesh->verts[i]->diameter < 2 * k)
		{
			mesh->verts[i]->color[0] = 0;
			mesh->verts[i]->color[1] = 1;
			mesh->verts[i]->color[2] = 1;
		}
		else if (mesh->verts[i]->diameter >= 2 * k && mesh->verts[i]->diameter < 3 * k)
		{
			mesh->verts[i]->color[0] = 0;
			mesh->verts[i]->color[1] = 0;
			mesh->verts[i]->color[2] = 1;
		}
		else if (mesh->verts[i]->diameter >= 3 * k && mesh->verts[i]->diameter < 4 * k)
		{
			mesh->verts[i]->color[0] = 1;
			mesh->verts[i]->color[1] = 0;
			mesh->verts[i]->color[2] = 0;
		}
		else
		{
			mesh->verts[i]->color[0] = 1;
			mesh->verts[i]->color[1] = 0;
			mesh->verts[i]->color[2] = 1;
		}
	}
}

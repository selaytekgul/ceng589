#include "Segmentor.h"
#include "VectorMath.h"
#include "TriangleMeshMath.h"

Segmentor::Segmentor(Mesh* mesh)
	: mesh(mesh)
{}

Segmentor::~Segmentor() = default;

void Segmentor::assignDiameterValuesOfVertices()
{
	//loop through the triangles to trace each vertex, find normals, draw rays, calculate and add diameters to vertex's diameter attribute
	for (size_t triangleIndex = 0; triangleIndex < mesh->tris.size(); triangleIndex++)
	{
		//get the vertex index number of the vertices of the triangle at hand
		triVertsIds vertexIdsOfTriangle = TriangleMeshMath::getVertexIdsOfTriangle(mesh, triangleIndex);

		//get the coordinates of the vertices of the triangle at hand (by using the vertex index numbers)
		triVertsCoords coordinatesOfVerticesOfTriangle = TriangleMeshMath::getCoordsOfTriangle(mesh, vertexIdsOfTriangle);
		
		//for each vertex of base triangle
		for (size_t selectedVertexNumber = 0; selectedVertexNumber < 3; selectedVertexNumber++)
		{
			//d is the normal vector from the selected vertex of the base triangle drawn according to the other two vertices
			float d[3];
			calculateCrossProduct(d, coordinatesOfVerticesOfTriangle, selectedVertexNumber);

			//p is the selected vertex of the base triangle
			float p[3];
			TD::fillWith(p, coordinatesOfVerticesOfTriangle[selectedVertexNumber], coordinatesOfVerticesOfTriangle[selectedVertexNumber].size());

			calculateShortestDiameter(triangleIndex, selectedVertexNumber, vertexIdsOfTriangle, p, d);
		}
	}//loop through the triangles to trace each vertex, find normals, draw rays, calculate and add diameters to vertex's diameter attribute
	mesh->setMinMaxDiameters();
	mesh->discardInfAndNegativeDiameters();
}
void Segmentor::calculateCrossProduct(float crossProductVector[3], const triVertsCoords& coordinatesOfVerticesOfTriangle, const size_t selectedVertexNumber)
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
	VectorMath::crossProduct(crossProductVector, vectorsToTheOtherVerticesArray[0], vectorsToTheOtherVerticesArray[1]);
}

void Segmentor::calculateShortestDiameter(const int triangleIndex, const int selectedVertexNumber, const triVertsIds& vertexIdsOfTriangle, const float p[3], const float d[3])
{
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
		
		//set the diameter attribute
		const float previousLength = selectedVertex->diameter;
		if (intersectionLength > 0 && (previousLength <= 0 || intersectionLength < previousLength))
			selectedVertex->diameter = intersectionLength;
	}
}

void Segmentor::setColorValuesToVertices()
{
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

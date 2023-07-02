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
				for (size_t coordinate = 0; coordinate < 3; coordinate++) {
					vectorsToTheOtherVerticesArray[otherVertexNumber][coordinate] = vectorsToTheOtherVertices[otherVertexNumber][coordinate];
				}
			}//store two vectors from the selected vertices in a 2D array

			//calculate cross product of the two vectors
			float crossProductVector[3];
			VectorMath::crossProduct(vectorsToTheOtherVerticesArray[0], vectorsToTheOtherVerticesArray[1], crossProductVector);

			//p is the selected vertex of the base triangle
			float p[3];
			p[0] = coordinatesOfVerticesOfTriangle[selectedVertexNumber][0];
			p[1] = coordinatesOfVerticesOfTriangle[selectedVertexNumber][1];
			p[2] = coordinatesOfVerticesOfTriangle[selectedVertexNumber][2];

			//d is the normal vector from the selected vertex of the base triangle drawn according to the other two vertices
			float d[3];
			d[0] = crossProductVector[0];
			d[1] = crossProductVector[1];
			d[2] = crossProductVector[2];

			//add the normal vector coordinates to the mesh->verts[].normals attribute
			mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->normalList.push_back(d); //NOT USED

			selectATargetTriangle(mesh, triangleIndex, selectedVertexNumber, vertexIdsOfTriangle, p, d);
		}//select a vertex from 3 vertices of triangle
	}//loop through the triangles to trace each vertex, find normals, draw rays, calculate and add lengths to vertex's length attribute
	mesh->setMinMaxLenghts();
	mesh->discardInfAndNegativeLenghts();
}

void Segmentor::selectATargetTriangle(Mesh* mesh, int triangleIndex, int selectedVertexNumber, triVertsIds vertexIdsOfTriangle, float p[3], float d[3]) {
	//select a target triangle
	for (size_t targetTriangleIndex = 0; targetTriangleIndex < mesh->tris.size(); targetTriangleIndex++)
	{
		//make sure that the target triangle is not the same as the base triangle
		if (triangleIndex != targetTriangleIndex) {
			triVertsIds vertexIdsOfTargetTriangle = TriangleMeshMath::getVertexIdsOfTriangle(mesh, targetTriangleIndex);

			//get the coordinates of the vertices of the target triangle at hand (by using the vertex index numbers)
			triVertsCoords coordinatesOfVerticesOfTargetTriangle;
			for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
			{
				for (size_t coordinate = 0; coordinate < 3; coordinate++) {
					coordinatesOfVerticesOfTargetTriangle[vertexNumber][coordinate] = mesh->verts[vertexIdsOfTargetTriangle[vertexNumber]]->coords[coordinate];
				}
			}//get the coordinates of the vertices of the target triangle at hand (by using the vertex index numbers)

			//fill the 2D array of targetTriangleVertices with the coordinate values
			float targetTriangleVertices[3][3];
			for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
			{
				for (size_t coordinate = 0; coordinate < 3; coordinate++) {
					targetTriangleVertices[vertexNumber][coordinate] = coordinatesOfVerticesOfTargetTriangle[vertexNumber][coordinate];
				}
			}//fill the 2D array of targetTriangleVertices with the coordinate values

			//calculate the lengths
			const float intersectionLength = VectorMath::rayIntersectsTriangle(p, d, targetTriangleVertices[0], targetTriangleVertices[1], targetTriangleVertices[2]);
			const float previousLength = mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->length;
			if (intersectionLength > 0
				&& (previousLength <= 0 || intersectionLength < previousLength))
			{
				mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->length = intersectionLength;
				mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->numberOfLenghtsContributed++;
				//long length vertices are marked
				if (intersectionLength > 80) {
					mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->hasLongLength = true;
					mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->intersectionTriangleIdsList.push_back(targetTriangleIndex);
					mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->intersectionNormalsList.push_back(d);
					mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->intersectionTrianglesVertexIdsList.push_back(vertexIdsOfTargetTriangle);
				}
			}
		}//make sure that the target triangle is not the same as the base triangle
	}//select a target triangle
}
void Segmentor::setColorValuesToVertices(Mesh* mesh) {
	//assign color values to the vertices according to their length values
	for (int i = 0; i < (int)mesh->verts.size(); i++)
	{
		float k = Vertex::maxLength / 8.0f;
		if (mesh->verts[i]->length < k)
		{
			mesh->verts[i]->color[0] = 0;
			mesh->verts[i]->color[1] = 1;
			mesh->verts[i]->color[2] = 0;
		}
		else if (mesh->verts[i]->length >= k && mesh->verts[i]->length < 2 * k)
		{
			mesh->verts[i]->color[0] = 0;
			mesh->verts[i]->color[1] = 1;
			mesh->verts[i]->color[2] = 1;
		}
		else if (mesh->verts[i]->length >= 2 * k && mesh->verts[i]->length < 3 * k)
		{
			mesh->verts[i]->color[0] = 0;
			mesh->verts[i]->color[1] = 0;
			mesh->verts[i]->color[2] = 1;
		}
		else if (mesh->verts[i]->length >= 3 * k && mesh->verts[i]->length < 4 * k)
		{
			mesh->verts[i]->color[0] = 1;
			mesh->verts[i]->color[1] = 0;
			mesh->verts[i]->color[2] = 0;
		}
		else //if (mesh->verts[i]->length >= 4 * k && mesh->verts[i]->length < 5 * k)
		{
			mesh->verts[i]->color[0] = 1;
			mesh->verts[i]->color[1] = 0;
			mesh->verts[i]->color[2] = 1;
		}
		//else // if (mesh->verts[i]->length >= 4 * k && mesh->verts[i]->length < 5 * k)
		//{
		//	mesh->verts[i]->color[0] = 1;
		//	mesh->verts[i]->color[1] = 1;
		//	mesh->verts[i]->color[2] = 0;
		//}
	}
}

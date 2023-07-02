#include "Segmentor.h"
#include "VectorMath.h"

void Segmentor::assignLengthValuesOfVertices(Mesh* mesh)
{
	//fill the length attribute of each of the mesh->verts with -5
	for (size_t triangleIndex = 0; triangleIndex < mesh->tris.size(); triangleIndex++)
	{
		//get the vertex index number of the vertices of the triangle at hand
		std::array<int, 3> vertexIdsOfTriangle;
		vertexIdsOfTriangle[0] = mesh->tris[triangleIndex]->v1i;
		vertexIdsOfTriangle[1] = mesh->tris[triangleIndex]->v2i;
		vertexIdsOfTriangle[2] = mesh->tris[triangleIndex]->v3i;

		//fill the lenghts with -5
		for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
		{
			mesh->verts[vertexIdsOfTriangle[vertexNumber]]->length = -5.0;
		}
	}

	//loop through triangles to trace each vertex, find normals, draw rays, calculate and add lengths to vertex's length attribute
	for (size_t triangleIndex = 0; triangleIndex < mesh->tris.size(); triangleIndex++)
	{
		//get the vertex index number of the vertices of the triangle at hand
		std::array<int, 3> vertexIdsOfTriangle;
		vertexIdsOfTriangle[0] = mesh->tris[triangleIndex]->v1i;
		vertexIdsOfTriangle[1] = mesh->tris[triangleIndex]->v2i;
		vertexIdsOfTriangle[2] = mesh->tris[triangleIndex]->v3i;

		//get the coordinates of the vertices of the triangle at hand (by using the vertex index numbers)
		std::array<std::array<float,3>, 3> coordinatesOfVerticesOfTriangle;
		for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
		{
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				coordinatesOfVerticesOfTriangle[vertexNumber][coordinate] = mesh->verts[vertexIdsOfTriangle[vertexNumber]]->coords[coordinate];
			}
		}//get the coordinates of the vertices of the triangle at hand (by using the vertex index numbers)
		
		//select a vertex from 3 vertices of triangle
		for (size_t selectedVertexNumber = 0; selectedVertexNumber < 3; selectedVertexNumber++)
		{
			//find the other 2 vertices of the triangle
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

			//create two vectors from the selected vertex
			std::array<std::array<float, 3>, 2> vectorsToTheOtherVertices;
			for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
			{
				for (size_t coordinate = 0; coordinate < 3; coordinate++) {
					vectorsToTheOtherVertices[otherVertexNumber][coordinate] =
						coordinatesOfOtherVertices[otherVertexNumber][coordinate]
						- coordinatesOfVerticesOfTriangle[selectedVertexNumber][coordinate];
				}
			}//create two vectors from the selected vertex

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
			//printf("NEW(B-A): x=%f y=%f z=%f\n", vectorsToTheOtherVerticesArray[0][0], vectorsToTheOtherVerticesArray[0][1], vectorsToTheOtherVerticesArray[0][2]);
			//printf("NEW(C-A): x=%f y=%f z=%f\n", vectorsToTheOtherVerticesArray[1][0], vectorsToTheOtherVerticesArray[1][1], vectorsToTheOtherVerticesArray[1][2]);
			//printf("NEWCrossProduct x=%f y=%f z=%f\n", crossProductVector[0], crossProductVector[1], crossProductVector[2]);

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

			//select a target triangle
			for (size_t targetTriangleIndex = 0; targetTriangleIndex < mesh->tris.size(); targetTriangleIndex++)
			{
				//make sure that the target triangle is not the same as the base triangle
				if (triangleIndex != targetTriangleIndex) {
					std::array<int, 3> vertexIdsOfTargetTriangle;
					vertexIdsOfTargetTriangle[0] = mesh->tris[targetTriangleIndex]->v1i;
					vertexIdsOfTargetTriangle[1] = mesh->tris[targetTriangleIndex]->v2i;
					vertexIdsOfTargetTriangle[2] = mesh->tris[targetTriangleIndex]->v3i;

					//get the coordinates of the vertices of the target triangle at hand (by using the vertex index numbers)
					std::array<std::array<float, 3>, 3> coordinatesOfVerticesOfTargetTriangle;
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
		}//select a vertex from 3 vertices of triangle
	}//loop through triangles to trace each vertex, find normals, draw rays, calculate and add lengths to vertex's length attribute
	
	//loop through vertices (mesh->verts):
	//find min & max values
	Vertex::minLength = std::numeric_limits<float>::max();
	Vertex::maxLength = std::numeric_limits<float>::min();
	for (size_t selectedVertexIndex = 0; selectedVertexIndex < mesh->verts.size(); selectedVertexIndex++)
	{
		//find min & max values
		if (std::isfinite(mesh->verts[selectedVertexIndex]->length) && mesh->verts[selectedVertexIndex]->length > 0.0001f)
		{
			Vertex::maxLength = Vertex::maxLength > mesh->verts[selectedVertexIndex]->length ? Vertex::maxLength : mesh->verts[selectedVertexIndex]->length;
			Vertex::minLength = Vertex::minLength < mesh->verts[selectedVertexIndex]->length ? Vertex::minLength : mesh->verts[selectedVertexIndex]->length;
		}
	}//select a vertex from mesh->verts to calculate the lenghts as average of recorded lengths

	//loop through vertices (mesh->verts):
	//select a vertex from mesh->verts to discard the inf values
	for (size_t selectedVertexIndex = 0; selectedVertexIndex < mesh->verts.size(); selectedVertexIndex++)
	{
		if (mesh->verts[selectedVertexIndex]->length < 0.0f)
		{
			mesh->verts[selectedVertexIndex]->length = Vertex::minLength;
		}
		if (std::isinf(mesh->verts[selectedVertexIndex]->length))
		{
			mesh->verts[selectedVertexIndex]->length = Vertex::maxLength;
		}
	}//select a vertex from mesh->verts to discard the inf values
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

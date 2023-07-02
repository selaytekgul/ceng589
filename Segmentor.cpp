#include "Segmentor.h"

/* a = b - c */
void Segmentor::vector(float v_A[], const float v_B[],const float v_C[]) {
	v_A[0] = v_B[0] - v_C[0];
	v_A[1] = v_B[1] - v_C[1];
	v_A[2] = v_B[2] - v_C[2];
}
void Segmentor::crossProduct(const float v_A[], const float v_B[], float CP[]) {
	CP[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
	CP[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
	CP[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

float Segmentor::innerProduct(const float v[], const float q[]) {
	return v[0] * q[0] + v[1] * q[1] + v[2] * q[2];
}

float Segmentor::calculateLengthOfVector(const float v[]) {
	float square = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	float root = sqrt(square);
	return root;
}

float Segmentor::rayIntersectsTriangle(float* p, float* d, float* v0, float* v1, float* v2) {
	float e1[3], e2[3], h[3], s[3], q[3];
	float a, f, u, v;
	vector(e1, v1, v0);
	vector(e2, v2, v0);
	crossProduct(d, e2, h);
	a = innerProduct(e1, h);
	if (a > -0.00001 && a < 0.00001)
		//return(false);
		return(-1);

	f = 1 / a;
	vector(s, p, v0);
	u = f * (innerProduct(s, h));

	if (u < 0.0 || u > 1.0)
		//return(false);
		return(-1);

	crossProduct(s, e1, q);
	v = f * innerProduct(d, q);

	if (v < 0.0 || u + v > 1.0)
		//return(false);
		return(-1);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	float t = f * innerProduct(e2, q);

	if (t > 0.00001) // ray intersection
	{
		const float directionVector[3] = {d[0], d[1], d[2]};
		const float lenghtOfDirectionVector = calculateLengthOfVector(directionVector);
		//return(true);
		return(t * lenghtOfDirectionVector);
	}
	else // this means that there is a line intersection
		 // but not a ray intersection
		//return (false);
		return (-2);
}

//void Segmentor::normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr) {
//	float minValue = std::numeric_limits<float>::max();
//	float maxValue = std::numeric_limits<float>::min();
//
//	// Find the minimum and maximum values in the input vector
//	for (float value : inputArr) {
//		if (value < minValue)
//			minValue = value;
//		if (value > maxValue)
//			maxValue = value;
//	}
//	int i = 0;
//	// Normalize the input vector values and store them in the output vector
//	for (float value : inputArr) {
//		float normalizedValue = (value - minValue) / (maxValue - minValue);
//		outputArr[i] = normalizedValue;
//		i++;
//	}
//}

void Segmentor::assignLengthValuesOfVertices(Mesh* mesh)
{
	int numberOfIntersections = 0;
	int numberOfNOTIntersections = 0;

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
			crossProduct(vectorsToTheOtherVerticesArray[0], vectorsToTheOtherVerticesArray[1], crossProductVector);
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
					const float intersectionLength = rayIntersectsTriangle(p, d, targetTriangleVertices[0], targetTriangleVertices[1], targetTriangleVertices[2]);
					const float previousLength = mesh->verts[vertexIdsOfTriangle[selectedVertexNumber]]->length;
					if (intersectionLength > 0
						&& (previousLength <= 0 || intersectionLength < previousLength))
					{
						numberOfIntersections++;
						//printf("p(%d) =%f, %f, %f\n", vertexIdsOfTriangle[selectedVertexNumber], p[0], p[1], p[2]);
						//printf("d =%f, %f, %f\n", d[0], d[1], d[2]);
						//printf("v0 =%f, %f, %f\n", targetTriangleVertices[0][0], targetTriangleVertices[0][1], targetTriangleVertices[0][2]);
						//printf("v1 =%f, %f, %f\n", targetTriangleVertices[1][0], targetTriangleVertices[1][1], targetTriangleVertices[1][2]);
						//printf("v2 =%f, %f, %f\n", targetTriangleVertices[2][0], targetTriangleVertices[2][1], targetTriangleVertices[2][2]);
						//printf("Intersects = %f \n\n\n\n", intersects);

						//add the lengths and keep the number of the intersections occured for each of the vertices
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
					else
					{
						//printf("NOT INT p(%d) =%f, %f, %f\n", vertexIdsOfTriangle[selectedVertexNumber], p[0], p[1], p[2]);
						//printf("NOT INT d =%f, %f, %f\n", d[0], d[1], d[2]);
						//printf("NOT INT v0 =%f, %f, %f\n", v0[0], v0[1], v0[2]);
						//printf("NOT INT v1 =%f, %f, %f\n", v1[0], v1[1], v1[2]);
						//printf("NOT INT v2 =%f, %f, %f\n", v2[0], v2[1], v2[2]);
						//printf("NOT INT Intersects = %f \n\n\n\n", intersects);
						numberOfNOTIntersections++;
					}
				}//make sure that the target triangle is not the same as the base triangle
			}//select a target triangle
		}//select a vertex from 3 vertices of triangle
	}//loop through triangles to trace each vertex, find normals, draw rays, calculate and add lengths to vertex's length attribute
	
	//loop through vertices (mesh->verts):
	//select a vertex from mesh->verts to calculate the lenghts as average of recorded lengths
	//find min & max values
	//float minLength = std::numeric_limits<float>::max();
	//float maxLength = std::numeric_limits<float>::min();
	Vertex::minLength = std::numeric_limits<float>::max();
	Vertex::maxLength = std::numeric_limits<float>::min();
	for (size_t selectedVertexIndex = 0; selectedVertexIndex < mesh->verts.size(); selectedVertexIndex++)
	{
		////calculate the lenghts as average of recorded lengths
		//mesh->verts[selectedVertexIndex]->length /= mesh->verts[selectedVertexIndex]->numberOfLenghtsContributed;
		//find min & max values
		if (std::isfinite(mesh->verts[selectedVertexIndex]->length) && mesh->verts[selectedVertexIndex]->length > 0.0001f)
		{
			Vertex::maxLength = Vertex::maxLength > mesh->verts[selectedVertexIndex]->length ? Vertex::maxLength : mesh->verts[selectedVertexIndex]->length;
			Vertex::minLength = Vertex::minLength < mesh->verts[selectedVertexIndex]->length ? Vertex::minLength : mesh->verts[selectedVertexIndex]->length;
		}
		//printf("p(%d) =%f\n", selectedVertexIndex, mesh->verts[selectedVertexIndex]->length);
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
		//printf("DISCARDED p(%d) =%f\n", selectedVertexIndex, mesh->verts[selectedVertexIndex]->length);
	}//select a vertex from mesh->verts to discard the inf values


	//printf("numberofintersections = %d\n", numberOfIntersections);
	//printf("numberofNOTintersections = %d\n", numberOfNOTIntersections);
}

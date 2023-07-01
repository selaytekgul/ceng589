#include "Painter.h"

/* a = b - c */
#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];

#define innerProduct(v,q) \
((v)[0] * (q)[0] + \
	(v)[1] * (q)[1] + \
	(v)[2] * (q)[2])

#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

void Painter::crossProductFunction(const float v_A[], const float v_B[], float CP[]) {
	CP[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
	CP[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
	CP[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

float Painter::rayIntersectsTriangle(float* p, float* d, float* v0, float* v1, float* v2) {
	float e1[3], e2[3], h[3], s[3], q[3];
	float a, f, u, v;
	vector(e1, v1, v0);
	vector(e2, v2, v0);
	crossProduct(h, d, e2);
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

	crossProduct(q, s, e1);
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

float Painter::calculateLengthOfVector(const float v[]) {
	float square = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	float root = sqrt(square);
	return root;
}

//void Painter::normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr) {
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

void Painter::assignLengthValuesOfVertices(Mesh* mesh)
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
			crossProductFunction(vectorsToTheOtherVerticesArray[0], vectorsToTheOtherVerticesArray[1], crossProductVector);
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

SoSeparator* Painter::getShapeSep(Mesh* mesh)
{
	SoSeparator* res = new SoSeparator();

	//transformation
	//not needed

	//color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(0, 1, 0); //paint all vertices with this color
	//mat->transparency = 0.5f : 0.0f; //0 makes it completely opaque, the default

	//std::vector<float> inputArray((int)mesh->verts.size());
	//std::vector<float> outputArray((int)mesh->verts.size());
	//for (int i = 0; i < (int)mesh->verts.size(); i++)
	//{
	//	float value = mesh->verts[i]->length;
	//	inputArray[i] = value;
	//	outputArray[i] = -9.0f;
	//}
	//normalizeArray(inputArray, outputArray);

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

	bool youWantToPaintEachVertexDifferently = false;
	youWantToPaintEachVertexDifferently = true;
	if (youWantToPaintEachVertexDifferently)
		for (int i = 0; i < (int)mesh->verts.size(); i++) //i = 0 obj->color above overwritten here
			mat->diffuseColor.set1Value(i, mesh->verts[i]->color); //vert color according to its x-y-z coord (for mesh1) and to the transferred color (for mesh2)
	res->addChild(mat);

	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14f;
	res->addChild(hints); //Gouraud shading

	if (youWantToPaintEachVertexDifferently)
	{
		SoMaterialBinding* materialBinding = new SoMaterialBinding; //for 2+ diffuse color usage on the same mesh
		materialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
		res->addChild(materialBinding);
	}

	//shape
	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->verts.size(); c++)
		coords->point.set1Value(c, mesh->verts[c]->coords[0], mesh->verts[c]->coords[1], mesh->verts[c]->coords[2]);
	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	for (int c = 0; c < mesh->tris.size(); c++)
	{
		faceSet->coordIndex.set1Value(c*4, mesh->tris[c]->v1i);
		faceSet->coordIndex.set1Value(c*4 + 1, mesh->tris[c]->v2i);
		faceSet->coordIndex.set1Value(c*4 + 2, mesh->tris[c]->v3i);
		faceSet->coordIndex.set1Value(c*4 + 3, -1);

		if (youWantToPaintEachVertexDifferently)
		{
			faceSet->materialIndex.set1Value(0 + 4*c, mesh->tris[c]->v1i);
			faceSet->materialIndex.set1Value(1 + 4*c, mesh->tris[c]->v2i);
			faceSet->materialIndex.set1Value(2 + 4*c, mesh->tris[c]->v3i);
			faceSet->materialIndex.set1Value(3 + 4*c, -1);
		}
	}
	res->addChild(coords);
	res->addChild(faceSet);

	return res;
}



/* stuff below are from my old projects; should run fine and be useful in your development

if (drawThickEdges) //draw thick edges (may be useful in geodesic path drawing)
	{
		SoSeparator* thickEdgeSep = new SoSeparator;
		//material
		SoMaterial* ma = new SoMaterial;
		ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 1.0f);
		thickEdgeSep->addChild(ma);
		SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 5.0f;	thickEdgeSep->addChild(sty);

		//shape
		SoIndexedLineSet* ils = new SoIndexedLineSet;
		SoCoordinate3* co = new SoCoordinate3;

		//assumes no edge in sedges is removed
		for (unsigned int se = 0; se < mesh->sedges.size(); se++)
		{
			SbVec3f end1 = mesh->verts[ mesh->sedges[se]->v1i ]->coords + SbVec3f(deltaX, 0.0f, 0.0f),
					end2 = mesh->verts[ mesh->sedges[se]->v2i ]->coords + SbVec3f(deltaX, 0.0f, 0.0f);
			co->point.set1Value(2*se, end1);
			co->point.set1Value(2*se + 1, end2);
		}

		for (unsigned int ci = 0; ci < mesh->sedges.size(); ci++)
		{
			ils->coordIndex.set1Value(3*ci, 2*ci);	ils->coordIndex.set1Value(3*ci + 1, 2*ci + 1);
			ils->coordIndex.set1Value(3*ci + 2, -1); //end this edge with -1
		}
		thickEdgeSep->addChild(co);	thickEdgeSep->addChild(ils);
		obj->sep->addChild(thickEdgeSep);
	}
	
	
SoSeparator* Painter::get1PointSep(ScreenObject* obj, int pnt, int drawWhat, float deltaX, float deltaY, float scale)
{
	//renders only 1 pnt in blue, w/ drawWhat = 1 for spectral coords, = 2 for spatial coords, = 5 for coord written here

	Mesh* mesh = obj->getMesh();

	SoSeparator* pntSep = new SoSeparator;
	//material
	SoMaterial* mat = new SoMaterial;
	if (mesh->targetMesh)
		mat->diffuseColor.setValue(SbColor(1.0f, 0.0f, 0.0f)); //red
	else
		mat->diffuseColor.setValue(SbColor(0.0f, 1.0f, 0.0f)); //green
if (pnt == 594) mat->diffuseColor.setValue(SbColor(1.0f, 0.0f, 1.0f)); //magenta
//if (pnt == 6916) mat->diffuseColor.setValue(SbColor(0.0f, 1.0f, 1.0f));

	pntSep->addChild(mat);
	SoDrawStyle* style = new SoDrawStyle;
	style->pointSize = 17.0f;
	pntSep->addChild(style);
	
	//shape
	SoVertexProperty* vp = new SoVertexProperty;
	if (drawWhat == 2)
		vp->vertex.set1Value(0, scale*mesh->verts[pnt]->coords[0]+deltaX, scale*mesh->verts[pnt]->coords[1]+deltaY, scale*mesh->verts[pnt]->coords[2]);
	else if (drawWhat == 5)
		vp->vertex.set1Value(0, 0.016721f, -0.000984876f, 0.0f);
	else
		vp->vertex.set1Value(0, scale*mesh->verts[pnt]->spectralK[0]+deltaX, scale*mesh->verts[pnt]->spectralK[1]+deltaX, scale*mesh->verts[pnt]->spectralK[2]+deltaX);
	SoPointSet* pSet = new SoPointSet;
	pSet->numPoints = 1;
	pSet->vertexProperty = vp;
	pntSep->addChild(pSet);

//cout << pnt << " ------> " << mesh->verts[pnt]->matchIdx << endl;
	return pntSep;
}

SoSeparator* Painter::getPointsSep(Mesh* mesh, SbColor c)
{
	//renders grid points, i.e. voxel centers

	SoSeparator* pntsSep = new SoSeparator;

	//material
	SoMaterial* mat = new SoMaterial;
	mat->diffuseColor.set1Value(0, c); //SbColor(199.0f/255.0f, 166.0f/255.0f, 1.0f));
	pntsSep->addChild(mat);
	SoDrawStyle* style = new SoDrawStyle;
	style->pointSize = 7.0f;
	pntsSep->addChild(style);

	//shape
	SoVertexProperty* vp = new SoVertexProperty;	
	int nPnts = (int) mesh->verts.size(), nAdds = 0;
	for (int p = 0; p < nPnts; p++)
	{
//if (selection[p]->center[1] > 120) continue; //just draw allowed voxels see its volume/thickness better
		vp->vertex.set1Value(nAdds++, mesh->verts[p]->coords);
	}
	SoPointSet* pSet = new SoPointSet;
	pSet->numPoints = nAdds;
	pSet->vertexProperty = vp;
	pntsSep->addChild(pSet);

	return pntsSep;
}

SoSeparator* Painter::getSpheresSep(Mesh* mesh, float deltaX, float deltaY, float scale)
{
	//returns a set of spheres to highlight each mesh.samples[i]

	SoSeparator* spheresSep = new SoSeparator();

	float radius = 50.0f;

	for (int i = 0; i < (int) mesh->samples.size(); i++)
	{
		//1 sphere for this sample
		SoSeparator* sphere1Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		tra->translation.setValue(scale*mesh->verts[ mesh->samples[i] ]->coords[0]+deltaX, scale*mesh->verts[ mesh->samples[i] ]->coords[1]+deltaY, scale*mesh->verts[ mesh->samples[i] ]->coords[2]);
		sphere1Sep->addChild(tra);

		//material
		SoMaterial* ma = new SoMaterial;
		if (i == 0)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.0f, 0.7f));
		else if (i == 1)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.0f, 0.0f));
		else if (i == 2)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.7f, 0.0f));
		else if (i == 3)
			ma->diffuseColor.setValue(SbColor(0.7f, 0.0f, 0.7f));
		else if (i == 4)
			ma->diffuseColor.setValue(SbColor(0.7f, 0.7f, 0.0f));
		else
			ma->diffuseColor.setValue(SbColor(0.7f, 0.0f, 0.0f));

		sphere1Sep->addChild(ma);

		//shape
		SoSphere* sph1 = new SoSphere();
		sph1->radius = radius;
		sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above

		spheresSep->addChild(sphere1Sep);
	}
	
	return spheresSep;
}
*/

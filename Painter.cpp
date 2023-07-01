#include "Painter.h"
struct Point3Di
{
	int x, y, z;
};

struct Point3Df
{
	float x, y, z;
};

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

void Painter::crossProductFunction(const int v_A[], const int v_B[], int c_P[]) {
	c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
	c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
	c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

float Painter::rayIntersectsTriangle(float* p, float* d,
	float* v0, float* v1, float* v2) {
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
		const float lenghtOfDirectionVector = calculateLength(directionVector);
		//return(true);
		return(t * lenghtOfDirectionVector);
	}
	else // this means that there is a line intersection
		 // but not a ray intersection
		//return (false);
		return (-2);
}


float Painter::calculateLength(const float v[]) {
	int square = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	float root = sqrt(square);
	return root;
}

void Painter::normalizeArray(const std::vector<float>& inputArr, std::vector<float>& outputArr) {
	float minValue = std::numeric_limits<float>::max();
	float maxValue = std::numeric_limits<float>::min();

	// Find the minimum and maximum values in the input vector
	for (float value : inputArr) {
		if (value < minValue)
			minValue = value;
		if (value > maxValue)
			maxValue = value;
	}
	int i = 0;
	// Normalize the input vector values and store them in the output vector
	for (float value : inputArr) {
		float normalizedValue = (value - minValue) / (maxValue - minValue);
		outputArr[i] = normalizedValue;
		i++;
	}
}

void Painter::triangleNormalVectors(Mesh* mesh)
{
	int numberOfIntersections = 0;
	int numberOfNOTIntersections = 0;
	for (size_t i = 0; i < mesh->tris.size(); i++)
	{
		int vertex1ID = mesh->tris[i]->v1i;
		int vertex2ID = mesh->tris[i]->v2i;
		int vertex3ID = mesh->tris[i]->v3i;
		mesh->verts[vertex1ID]->length = -5.0;
		mesh->verts[vertex2ID]->length = -5.0;
		mesh->verts[vertex3ID]->length = -5.0;
	}
	for (size_t i = 0; i < mesh->tris.size(); i++)
	{
		int vertex1ID = mesh->tris[i]->v1i;
		int vertex2ID = mesh->tris[i]->v2i;
		int vertex3ID = mesh->tris[i]->v3i;

		std::array<int, 3> vertexIdsOfTriangle;
		vertexIdsOfTriangle[0] = mesh->tris[i]->v1i;
		vertexIdsOfTriangle[1] = mesh->tris[i]->v2i;
		vertexIdsOfTriangle[2] = mesh->tris[i]->v3i;

		std::array<std::array<int,3>, 3> coordinatesOfVerticesOfTriangle;
		for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
		{
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				coordinatesOfVerticesOfTriangle[vertexNumber][coordinate] = mesh->verts[vertexIdsOfTriangle[vertexNumber]]->coords[coordinate];
			}
		}
		
		//select a vertex
		for (size_t selectedVertexNumber = 0; selectedVertexNumber < 3; selectedVertexNumber++)
		{
			//find the other 2 vertices of the triangle
			std::array<std::array<int, 3>, 2> coordinatesOfOtherVertices;
			int number = 0;
			for (size_t otherVertexNumber = 0; otherVertexNumber < 3; otherVertexNumber++)
			{
				if (selectedVertexNumber != otherVertexNumber)
				{
					for (size_t coordinate = 0; coordinate < 3; coordinate++) {
						coordinatesOfOtherVertices[number][coordinate] = coordinatesOfVerticesOfTriangle[otherVertexNumber][coordinate];
					}
					number++;
				}
			}

			//create two vectors from the selected vertex
			std::array<std::array<int, 3>, 2> vectorsToTheOtherVertices;
			for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
			{
				for (size_t coordinate = 0; coordinate < 3; coordinate++) {
					vectorsToTheOtherVertices[otherVertexNumber][coordinate] =
						coordinatesOfOtherVertices[otherVertexNumber][coordinate]
						- coordinatesOfVerticesOfTriangle[selectedVertexNumber][coordinate];
				}
			}
			//store two vectors from the selected vertices in a 2D array
			int vectorsToTheOtherVerticesArray[2][3];
			for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
			{
				for (size_t coordinate = 0; coordinate < 3; coordinate++) {
					vectorsToTheOtherVerticesArray[otherVertexNumber][coordinate] = vectorsToTheOtherVertices[otherVertexNumber][coordinate];
				}
			}

			//calculate cross product of the two vectors
			int crossProductVector[3];
			crossProductFunction(vectorsToTheOtherVerticesArray[0], vectorsToTheOtherVerticesArray[1], crossProductVector);
			printf("NEW(B-A): x=%d y=%d z=%d\n", vectorsToTheOtherVerticesArray[0][0], vectorsToTheOtherVerticesArray[0][1], vectorsToTheOtherVerticesArray[0][2]);
			printf("NEW(C-A): x=%d y=%d z=%d\n", vectorsToTheOtherVerticesArray[1][0], vectorsToTheOtherVerticesArray[1][1], vectorsToTheOtherVerticesArray[1][2]);
			printf("NEWCrossProduct x=%d y=%d z=%d\n", crossProductVector[0], crossProductVector[1], crossProductVector[2]);
		}

		int x1_A = mesh->verts[vertex1ID]->coords[0];
		int y1_A = mesh->verts[vertex1ID]->coords[1];
		int z1_A = mesh->verts[vertex1ID]->coords[2];

		int x2_B = mesh->verts[vertex2ID]->coords[0];
		int y2_B = mesh->verts[vertex2ID]->coords[1];
		int z2_B = mesh->verts[vertex2ID]->coords[2];

		int x3_C = mesh->verts[vertex3ID]->coords[0];
		int y3_C = mesh->verts[vertex3ID]->coords[1];
		int z3_C = mesh->verts[vertex3ID]->coords[2];

		int x_B_A = x2_B - x1_A;
		int y_B_A = y2_B - y1_A;
		int z_B_A = z2_B - z1_A;

		int x_C_A = x3_C - x1_A;
		int y_C_A = y3_C - y1_A;
		int z_C_A = z3_C - z1_A;

		int v_B_A[3] = { x_B_A, y_B_A, z_B_A };
		int v_C_A[3] = { x_C_A, y_C_A, z_C_A };
		int c_P[3];
		crossProductFunction(v_B_A, v_C_A, c_P);
		printf("(B-A): x=%d y=%d z=%d\n", x_B_A, y_B_A, z_B_A);
		printf("(C-A): x=%d y=%d z=%d\n", x_C_A, y_C_A, z_C_A);
		printf("CrossProduct x=%d y=%d z=%d\n", c_P[0], c_P[1], c_P[2]);

		float p[3];
		p[0] = x1_A;
		p[1] = y1_A;
		p[2] = z1_A;

		float d[3];
		d[0] = c_P[0];
		d[1] = c_P[1];
		d[2] = c_P[2];

		for (size_t j = 0; j < mesh->tris.size(); j++)
		{
			if (i != j) {
				int idxofVertex0inTriangle = mesh->tris[j]->v1i;
				int idxofVertex1inTriangle = mesh->tris[j]->v2i;
				int idxofVertex2inTriangle = mesh->tris[j]->v3i;
				float v0[3];
				float v1[3];
				float v2[3];

				float x0 = mesh->verts[idxofVertex0inTriangle]->coords[0];
				float x1 = mesh->verts[idxofVertex1inTriangle]->coords[0];
				float x2 = mesh->verts[idxofVertex2inTriangle]->coords[0];
				v0[0] = x0;
				v1[0] = x1;
				v2[0] = x2;


				float y0 = mesh->verts[idxofVertex0inTriangle]->coords[1];
				float y1 = mesh->verts[idxofVertex1inTriangle]->coords[1];
				float y2 = mesh->verts[idxofVertex2inTriangle]->coords[1];
				v0[1] = y0;
				v1[1] = y1;
				v2[1] = y2;

				float z0 = mesh->verts[idxofVertex0inTriangle]->coords[2];
				float z1 = mesh->verts[idxofVertex1inTriangle]->coords[2];
				float z2 = mesh->verts[idxofVertex2inTriangle]->coords[2];
				v0[2] = z0;
				v1[2] = z1;
				v2[2] = z2;


				float intersects = rayIntersectsTriangle(p, d, v0, v1, v2);
				if (intersects > 0) {
					numberOfIntersections++;
					printf("p =%f, %f, %f\n", p[0], p[1], p[2]);
					printf("d =%f, %f, %f\n", d[0], d[1], d[2]);
					printf("v0 =%f, %f, %f\n", v0[0], v0[1], v0[2]);
					printf("v1 =%f, %f, %f\n", v1[0], v1[1], v1[2]);
					printf("v2 =%f, %f, %f\n", v2[0], v2[1], v2[2]);
					printf("Intersects = %f \n\n\n\n", intersects);
					if (mesh->verts[vertex1ID]->length <= intersects) {
						mesh->verts[vertex1ID]->length = intersects;
						mesh->verts[vertex2ID]->length = intersects;
						mesh->verts[vertex3ID]->length = intersects;
					}
				}
				else
				{
					numberOfNOTIntersections++;
				}
			}
		}
	}
	printf("numberofintersections = %d\n", numberOfIntersections);
	printf("numberofNOTintersections = %d\n", numberOfNOTIntersections);
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

	std::vector<float> inputArray((int)mesh->verts.size());
	std::vector<float> outputArray((int)mesh->verts.size());
	for (int i = 0; i < (int)mesh->verts.size(); i++)
	{
		float value = mesh->verts[i]->length;
		inputArray[i] = value;
		outputArray[i] = -9.0f;
	}
	normalizeArray(inputArray, outputArray);

	//for (int i = 0; i < (int)mesh->verts.size(); i++) //i = 0 obj->color above overwritten here
	for (int i = 0; i < (int)outputArray.size(); i++)
	{
		//const float r = i % 2 == 0 ? 1 : 0;
		//const float g = i % 2 != 0 ? 1 : 0;
		//const float b = i % 3 == 0 ? 1 : 0;
		//float r = 0;
		//float g = 1;
		//float b = 0;
		////if (i < (int)mesh->verts.size()/2)
		////if (mesh->verts[i]->length < 0)
		//if (outputArray[i] < 0)
		//{
		//	r = 1;
		//	g = 1;
		//	b = 1;
		//}
		////else if (mesh->verts[i]->length >= 0 && mesh->verts[i]->length < 0.5)
		//else if (outputArray[i] >= 0 && outputArray[i] < 0.02)
		//{
		//	r = 0;
		//	g = 1;
		//	b = 0;
		//}
		//else if (outputArray[i] >= 0.02 && outputArray[i] < 0.05)
		//{
		//	r = 0;
		//	g = 1;
		//	b = 0;
		//}

		//else if (outputArray[i] >= 0.05 && outputArray[i] < 0.15)
		//{
		//	r = 0;
		//	g = 0;
		//	b = 1;
		//}
		//else if (outputArray[i] >= 0.15 && outputArray[i] < 0.25)
		//{
		//	r = 0;
		//	g = 0;
		//	b = 1;
		//}
		//else if (outputArray[i] >= 0.25 && outputArray[i] < 0.5)
		//{
		//	r = 0;
		//	g = 0;
		//	b = 1;
		//}
		//else if (outputArray[i] >= 0.5 && outputArray[i] < 0.75)
		//{
		//	r = 1;
		//	g = 0;
		//	b = 0;
		//}
		//else
		//{
		//	r = 1;
		//	g = 0;
		//	b = 0;
		//}
		//mesh->verts[i]->color[0] = r;
		//mesh->verts[i]->color[1] = g;
		//mesh->verts[i]->color[2] = b;

		mesh->verts[i]->color[0] = 0;
		mesh->verts[i]->color[1] = outputArray[i];
		mesh->verts[i]->color[2] = 0;
	}
	bool youWantToPaintEachVertexDifferently = false;
	youWantToPaintEachVertexDifferently = true;
	if (youWantToPaintEachVertexDifferently)
		for (int i = 0; i < (int)mesh->verts.size(); i++) //i = 0 obj->color above overwritten here
			mat->diffuseColor.set1Value(i, mesh->verts[i]->color); //vert color according to its x-y-z coord (for mesh1) and to the transferred color (for mesh2)

	//float green[3] = {0, 1, 0};
	//mat->diffuseColor.set1Value(0, green); //vert color according to its x-y-z coord (for mesh1) and to the transferred color (for mesh2)

	res->addChild(mat);

	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
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

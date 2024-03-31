#include "Mesh.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
float Vertex::minDiameter = std::numeric_limits<float>::max();
float Vertex::maxDiameter = std::numeric_limits<float>::min();

void Mesh::loadOff(std::string name)
{
	FILE* fPtr = fopen(name.c_str(), "r");
	char str[334];

	fscanf(fPtr, "%s", str);

	int nVerts, nTris, n, i = 0;
	float x, y, z;

	fscanf(fPtr, "%d %d %d\n", &nVerts, &nTris, &n);
	while (i++ < nVerts)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addVertex(x, y, z);
	}

	while (fscanf(fPtr, "%d", &i) != EOF)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addTriangle((int) x, (int) y, (int) z);
	}

	fclose(fPtr);
}

void Mesh::createCube(float sideLen)
{
	//coordinates
	float flbc[3] = {0, 0, 0}, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
			case 1:
				deltaX = sideLen;
				break;
			case 2:
				deltaZ = -sideLen;
				break;
			case 3:
				deltaX = 0;
				break;
			case 4:
				deltaZ = 0;
				deltaY = sideLen;
				break;
			case 5:
				deltaX = sideLen;
				break;
			case 6:
				deltaZ = -sideLen;
				break;
			default:
				deltaX = 0;;
				break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addTriangle(0, 2, 1);
	addTriangle(0, 3, 2);

	addTriangle(1, 2, 5);
	addTriangle(2, 6, 5);

	addTriangle(2, 3, 6);
	addTriangle(3, 7, 6);

	addTriangle(3, 4, 7);
	addTriangle(3, 0, 4);

	addTriangle(4, 5, 6);
	addTriangle(4, 6, 7);

	addTriangle(0, 1, 5);
	addTriangle(0, 5, 4);
}

void Mesh::createDoubleOpenCube(float sideLen)
{
	//coordinates
	float flbc[3] = {0, 0, 0}, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
			case 1:
				deltaX = sideLen;
				break;
			case 2:
				deltaZ = -sideLen;
				break;
			case 3:
				deltaX = 0;
				break;
			case 4:
				deltaZ = 0;
				deltaY = sideLen;
				break;
			case 5:
				deltaX = sideLen;
				break;
			case 6:
				deltaZ = -sideLen;
				break;
			default:
				deltaX = 0;;
				break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addTriangle(0, 2, 1);
	addTriangle(0, 3, 2);

	addTriangle(1, 2, 5);
	addTriangle(2, 6, 5);

	addTriangle(2, 3, 6);
	addTriangle(3, 7, 6);

	addTriangle(3, 4, 7);
	addTriangle(3, 0, 4);

	addTriangle(4, 5, 6);
	addTriangle(4, 6, 7);

	//addTriangle(0, 1, 5);
	//addTriangle(0, 5, 4);
}

void Mesh::createOpenCube(float sideLen)
{
	//coordinates
	float flbc[3] = {0, 0, 0}, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
			case 1:
				deltaX = sideLen;
				break;
			case 2:
				deltaZ = -sideLen;
				break;
			case 3:
				deltaX = 0;
				break;
			case 4:
				deltaZ = 0;
				deltaY = sideLen;
				break;
			case 5:
				deltaX = sideLen;
				break;
			case 6:
				deltaZ = -sideLen;
				break;
			default:
				deltaX = 0;;
				break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addTriangle(0, 2, 1);
	addTriangle(0, 3, 2);

	addTriangle(1, 2, 5);
	addTriangle(2, 6, 5);

	addTriangle(2, 3, 6);
	addTriangle(3, 7, 6);

	addTriangle(3, 4, 7);
	addTriangle(3, 0, 4);

	addTriangle(4, 5, 6);
	addTriangle(4, 6, 7);

	//addTriangle(0, 1, 5);
	//addTriangle(0, 5, 4);
}

void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	tris.push_back( new Triangle(idx, v1, v2, v3) );

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	if (! makeVertsNeighbor(v1, v2, idx))
		addEdge(v1, v2, idx);

	if (!makeVertsNeighbor(v1, v3, idx))
		addEdge(v1, v3, idx);

	if (!makeVertsNeighbor(v2, v3, idx))
		addEdge(v2, v3, idx);

}

bool Mesh::makeVertsNeighbor(int v1i, int v2i, int triIdx)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < verts[v1i]->vertList.size(); i++)
	{
		if (verts[v1i]->vertList[i] == v2i)
		{
			modifyEdge(v1i, v2i, triIdx);
			return true;
		}
	}


	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return false;
}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = verts.size();
	float* c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;

	verts.push_back( new Vertex(idx, c) );
}

void Mesh::addEdge(int v1, int v2, int triIdx)
{
	int idx = edges.size();
	Edge* edge = new Edge(idx, v1, v2);
	edge->existedTriangeNumber++;
	edge->triList.push_back(triIdx);
	edges.push_back(edge);

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

void Mesh::modifyEdge(int v1, int v2, int triIdx)
{
	//int edges_size = edges.size();
	//for (size_t i = 0; i < edges_size; i++)
	//{
	//	if ((edges[i]->v1i == v1 && edges[i]->v2i == v2) || (edges[i]->v2i == v1 && edges[i]->v1i == v2))
	//	{
	//		edges[i]->existedTriangeNumber++;
	//	}
	//}
	for (const auto neighEdge : verts[v1]->edgeList)
	{
		if (edges[neighEdge]->v1i == v1 && edges[neighEdge]->v2i == v2
			|| edges[neighEdge]->v2i == v1 && edges[neighEdge]->v1i == v2)
		{
			edges[neighEdge]->existedTriangeNumber++;
			edges[neighEdge]->triList.push_back(triIdx);
		}
	}
}

void Mesh::setMinMaxDiameters()
{
	//loop through vertices (this->verts): find min & max values
	for (Vertex* selectedVertex : this->verts)
	{
		//find min & max values
		if (std::isfinite(selectedVertex->diameter) && selectedVertex->diameter > 0.0001f)
		{
			Vertex::maxDiameter = std::max(Vertex::maxDiameter, selectedVertex->diameter);
			Vertex::minDiameter = std::min(Vertex::minDiameter, selectedVertex->diameter);
		}
	}
}

void Mesh::discardInfAndNegativeDiameters()
{
	//loop through vertices (this->verts): select a vertex from this->verts to discard the negative & inf values
	for(Vertex* selectedVertex : this->verts)
	{
		if (selectedVertex->diameter < 0.0f)
			selectedVertex->diameter = Vertex::minDiameter;


		if (std::isinf(selectedVertex->diameter))
			selectedVertex->diameter = Vertex::maxDiameter;
	}
}


void Mesh::computeLength(int edgeIdx)
{
	int endP1 = edges[edgeIdx]->v1i;
	int endP2 = edges[edgeIdx]->v2i;
	const float* v1coords = verts[endP1]->coords;
	const float* v2coords = verts[endP2]->coords;
	edges[edgeIdx]->length = VectorMath::distanceBetweenVectors(v1coords, v2coords);
}

void Mesh::windingNumberByYusufSahillioglu(Vertex* pnt)
{
	//computes generalized winding number (eq. 5 in Jacobson'13) of pnt to see whether it is inside the mesh (winding=1) or not (winding=0); holds for watertight meshes (still well-behaved otherwise)

	double a[3], b[3], c[3], aLen, bLen, cLen, twoPI = 2.0 * M_PI;
	pnt->winding = 0.0;
	for (int t = 0; t < (int)tris.size(); t++)
	{
		for (int ci = 0; ci < 3; ci++)
		{
			a[ci] = verts[tris[t]->v1i]->coords[ci] - pnt->coords[ci];
			b[ci] = verts[tris[t]->v2i]->coords[ci] - pnt->coords[ci];
			c[ci] = verts[tris[t]->v3i]->coords[ci] - pnt->coords[ci];
		}
		aLen = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		bLen = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
		cLen = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
		double atan2val= atan2(//determinant(a, b, c), writing the formula below instead of calling a function
			a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]),
			aLen * bLen * cLen + cLen * (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) + aLen * (c[0] * b[0] + c[1] * b[1] + c[2] * b[2]) + bLen * (a[0] * c[0] + a[1] * c[1] + a[2] * c[2]));
		//pnt->winding += atan2(5.0, 2.0); //get the value of tan-1(5.0 / 2.0) which is 1.19029 (return in radians: [-pi,pi], 0 if both 2 params=0.0)
		pnt->winding += atan2val;
	}
	float compare1 = pnt->winding;
	float compare2 = twoPI;
	//if (pnt->winding >= twoPI)
	if (abs(pnt->winding - twoPI) < 0.0001)
		pnt->winding = 1.0; //inside
	else
		pnt->winding = 0.0; //outside
}
#include "Mesh.h"

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
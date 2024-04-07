#include "Mesh.h"
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline triVertsIds getVertexIdsOfTriangleMesh(const Mesh* mesh, const int triangleIndex);
inline triVertsCoords getCoordsOfTriangleMesh(const Mesh* mesh, triVertsIds vertexIdsOfTriangle);
inline triOtherVertsCoords getOtherCoordsOfTriangleMesh(const triVertsCoords& coordsOfVerticesOfTriangle, const int selectedVertexNumber);
inline triOtherVertsCoords getVectorsToTheOtherVerticesMesh(const triOtherVertsCoords& coordsOfOtherVertices, const int selectedVertexNumber, const triVertsCoords& coordsOfVerticesOfTriangle);

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

void Mesh::computeLength(int edgeIdx)
{
	int endP1 = edges[edgeIdx]->v1i;
	int endP2 = edges[edgeIdx]->v2i;
	const float* v1coords = verts[endP1]->coords;
	const float* v2coords = verts[endP2]->coords;
	edges[edgeIdx]->length = VectorMath::distanceBetweenVectors(v1coords, v2coords);
}

void Mesh::computeDistFromEdgeMidToEndPntsTangPla(int edgeIdx)
{
	float mid[3] = { 0.0f, 0.0f, 0.0f };
	VectorMath::midpoint(mid, verts[edges[edgeIdx]->v1i]->coords, verts[edges[edgeIdx]->v2i]->coords);

	const float dist1 = VectorMath::calculateDistanceToTangentPlane(mid, verts[edges[edgeIdx]->v1i]->tangent_plane_normal, verts[edges[edgeIdx]->v1i]->tangent_plane_origin);
	const float dist2 = VectorMath::calculateDistanceToTangentPlane(mid, verts[edges[edgeIdx]->v2i]->tangent_plane_normal, verts[edges[edgeIdx]->v2i]->tangent_plane_origin);

	//std::cout << "edge id: " << edges[edgeIdx]->edge_idx << " distance: " << dist1 + dist2 << std::endl;
	edges[edgeIdx]->midToEndPointTangentPlanesDist = dist1 + dist2;
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

void Mesh::collapseEdgeTo(Edge* edge, int tovi)
{
	//delete edge
	if (edges[edge->edge_idx]->deleted == false)
	{
		numDeletedEdge++;
		edges[edge->edge_idx]->deleted = true;
	}

	//delete triangles
	for (size_t t = 0; t < edge->triList.size(); t++)
	{
		int trid = edge->triList[t];
		if (tris[trid]->deleted == false)
		{
			numDeletedTri++;
			tris[trid]->deleted = true;
		}
	}

	const int endP1i = edge->v1i;
	const int endP2i = edge->v2i;

	int fromvid;
	int tovid;

	if (endP1i == tovi)
	{
		fromvid = endP2i;
		tovid = endP1i;
	}
	else //if (endP2i == tovi)
	{
		fromvid = endP1i;
		tovid = endP2i;
	}

	//delete vertex
	if (verts[fromvid]->deleted == false)
	{
		verts[fromvid]->deleted = true;
		numDeletedVert++;
	}
	
	//modify connected triangles
	for (size_t t = 0; t < verts[fromvid]->triList.size(); t++)
	{
		int trid = verts[fromvid]->triList[t];
		if (tris[trid]->deleted == true)
			continue;

		if (fromvid == tris[trid]->v1i)
		{
			tris[trid]->v1i = tovid;
		}
		else if (fromvid == tris[trid]->v2i)
		{
			tris[trid]->v2i = tovid;
		}
		else // if (fromvid == tris[trid]->v3i)
		{
			tris[trid]->v3i = tovid;
		}
	}

	//modify connected edges
	for (size_t e = 0; e < verts[fromvid]->edgeList.size(); e++)
	{
		int edgeid = verts[fromvid]->edgeList[e];
		if (edges[edgeid]->deleted == true)
			continue;
		if (fromvid == edges[edgeid]->v1i)
		{
			edges[edgeid]->v1i = tovid;
		}
		else //if (fromvid == edges[edgeid]->v2i)
		{
			edges[edgeid]->v2i = tovid;
		}
	}
};

bool Mesh::collapseEdge(Edge* edge, std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>>* minHeap)
{
	int fromvid;
	int tovid;
	float mid_coords[3];
	
	{
		const int endP1i = edge->v1i;
		const int endP2i = edge->v2i;
		VectorMath::midpoint(mid_coords, verts[endP1i]->coords, verts[endP2i]->coords);
		fromvid = endP1i;
		tovid = endP2i;

		if (tovid == fromvid)
		{
			int b = 8;
			edge->deleted = true;
			return false;
		}
	}

	{
		int numSameNeighVerts = 0;
		for (size_t i = 0; i < verts[fromvid]->vertList.size(); i++)
		{
			int neighIdx = verts[fromvid]->vertList[i];
			for (size_t j = 0; j < verts[tovid]->vertList.size(); j++)
			{
				int neighIdx2 = verts[tovid]->vertList[j];
				if (neighIdx == neighIdx2 && verts[neighIdx2]->deleted == false && verts[neighIdx]->deleted == false)
				{
					numSameNeighVerts++;
					break;
				}
			}
		}

		if (numSameNeighVerts >= 3)
			return false;
	}

	//delete vertex
	if (verts[fromvid]->deleted == false)
	{
		verts[fromvid]->deleted = true;
		numDeletedVert++;
	}

	//delete edge
	if (edges[edge->edge_idx]->deleted == false)
	{
		numDeletedEdge++;
		edges[edge->edge_idx]->deleted = true;
	}
	else {
		int alreadydeletededge = 9;
	}

	//delete triangles
	for (size_t t = 0; t < edge->triList.size(); t++)
	{
		int trid = edge->triList[t];
		if (tris[trid]->deleted == false)
		{
			numDeletedTri++;
			tris[trid]->deleted = true;
		}
	}

	//modify remaining vertex coords
	verts[tovid]->coords[0] = mid_coords[0];
	verts[tovid]->coords[1] = mid_coords[1];
	verts[tovid]->coords[2] = mid_coords[2];


	//modify connected triangles
	for (size_t t = 0; t < verts[fromvid]->triList.size(); t++)
	{
		int trid = verts[fromvid]->triList[t];
		if (tris[trid]->deleted == true)
			continue;

		if (fromvid == tris[trid]->v1i)
		{
			tris[trid]->v1i = tovid;
		}
		else if (fromvid == tris[trid]->v2i)
		{
			tris[trid]->v2i = tovid;
		}
		else if (fromvid == tris[trid]->v3i)
		{
			tris[trid]->v3i = tovid;
		}
		else {
			int g = 0;
			continue;
		}
		auto it = std::find(verts[tovid]->triList.begin(), verts[tovid]->triList.end(), trid);
		if (it == verts[tovid]->triList.end()) {
			verts[tovid]->triList.push_back(trid);
		}
	}
	
	//modify connected edges
	for (size_t e = 0; e < verts[fromvid]->edgeList.size(); e++)
	{
		int edgeid = verts[fromvid]->edgeList[e];
		if (edges[edgeid]->deleted == true)
			continue;
		if (fromvid == edges[edgeid]->v1i)
		{
			edges[edgeid]->v1i = tovid;
		}
		else if (fromvid == edges[edgeid]->v2i)
		{
			edges[edgeid]->v2i = tovid;
		}
		else {
			int a = 9;
			continue;
		}
		auto it = std::find(verts[tovid]->edgeList.begin(), verts[tovid]->edgeList.end(), edgeid);
		if (it == verts[tovid]->edgeList.end()) {
			verts[tovid]->edgeList.push_back(edgeid);
		}
		//computeLength(edgeid);
		computeDistFromEdgeMidToEndPntsTangPla(edgeid);
		//minHeap->push({ edges[edgeid]->length, edgeid});
		minHeap->push({ edges[edgeid]->midToEndPointTangentPlanesDist, edges[edgeid]->edge_idx });
	}

	//modify connected edges
	for (size_t e = 0; e < verts[tovid]->edgeList.size(); e++)
	{
		int edgeid = verts[tovid]->edgeList[e];
		if (edges[edgeid]->deleted == true)
			continue;		
		//computeLength(edgeid);
		computeDistFromEdgeMidToEndPntsTangPla(edgeid);
		//minHeap->push({ edges[edgeid]->length, edgeid });
		minHeap->push({ edges[edgeid]->midToEndPointTangentPlanesDist, edges[edgeid]->edge_idx });
	}

	for (size_t i = 0; i < verts[fromvid]->vertList.size(); i++)
	{
		int vertid = verts[fromvid]->vertList[i];
		if (verts[vertid]->deleted)
			continue;
		{
			auto it = std::find(verts[tovid]->vertList.begin(), verts[tovid]->vertList.end(), vertid);
			if (it == verts[tovid]->vertList.end() && tovid != vertid) {
				verts[tovid]->vertList.push_back(vertid);
			}
		}
		
		{
			auto it2 = std::find(verts[vertid]->vertList.begin(), verts[vertid]->vertList.end(), tovid);
			if (it2 == verts[vertid]->vertList.end() && tovid != vertid) {
				verts[vertid]->vertList.push_back(tovid);
			}
		}
	}

	return true;
}

void Mesh::toOFF(const std::string& filename)
{
	const int numVerts = verts.size() - numDeletedVert;
	const int numEdges = edges.size() - numDeletedEdge;
	const int numTris = tris.size() - numDeletedTri;

	std::ofstream ofs(filename);
	if (!ofs.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	// Write the header for OFF file
	ofs << "OFF" << std::endl;
	//ofs << numVerts << " " << numTris << " 0" << std::endl;
	ofs << verts.size() << " " << numTris << " 0" << std::endl;

	for (size_t i = 0; i < verts.size(); i++)
	{
		ofs << verts[i]->coords[0] << " " << verts[i]->coords[1] << " " << verts[i]->coords[2] << std::endl;
	}
	
	for (size_t i = 0; i < tris.size(); i++)
	{
		if (tris[i]->deleted)
			continue;
		ofs << "3 " << tris[i]->v1i << " " << tris[i]->v2i << " " << tris[i]->v3i << std::endl;
	}

	ofs.close();
	std::cout << "OFF file saved: " << filename << std::endl;
}

void Mesh::inflatePoint(Mesh* mesh, Vertex* vert)
{
	if (vert->deleted)
		return;
	mesh->windingNumberByYusufSahillioglu(vert);
	if (vert->winding == 0.0f)
		return;
	do
	{
		//float* normal = returnPointNormal(vert);
		calculateVertexNormal(vert);
		float* normal = vert->point_normal;
		float alpha = 0.025;
		//float normal_to_be_added[3] = {0.0f, 0.0f, 0.0f};
		//normal_to_be_added[0] = normal[0] - vert->coords[0];
		//normal_to_be_added[1] = normal[1] - vert->coords[1];
		//normal_to_be_added[2] = normal[2] - vert->coords[2];
		//vert->coords[0] -= normal_to_be_added[0] * alpha;
		//vert->coords[1] -= normal_to_be_added[1] * alpha;
		//vert->coords[2] -= normal_to_be_added[2] * alpha;

		vert->coords[0] += normal[0] * alpha;
		vert->coords[1] += normal[1] * alpha;
		vert->coords[2] += normal[2] * alpha;

		mesh->windingNumberByYusufSahillioglu(vert);
		int a = 5;
	} while (vert->winding == 1.0f);
}

void Mesh::calculateNormalVectorMesh(float crossProductVector[3], const triVertsCoords& coordinatesOfVerticesOfTriangle, const size_t selectedVertexNumber)
{
	//find the other 2 vertices of the triangle
	triOtherVertsCoords coordinatesOfOtherVertices = getOtherCoordsOfTriangleMesh(coordinatesOfVerticesOfTriangle, selectedVertexNumber);

	//create two vectors from the selected vertex
	triOtherVertsCoords vectorsToTheOtherVertices = getVectorsToTheOtherVerticesMesh(coordinatesOfOtherVertices, selectedVertexNumber, coordinatesOfVerticesOfTriangle);

	//store two vectors from the selected vertices in a 2D array
	float vectorsToTheOtherVerticesArray[2][3];
	for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
		TD::fillWith(vectorsToTheOtherVerticesArray[otherVertexNumber], vectorsToTheOtherVertices[otherVertexNumber], 3);

	//calculate cross product of the two vectors
	VectorMath::crossProduct(crossProductVector, vectorsToTheOtherVerticesArray[0], vectorsToTheOtherVerticesArray[1]);
}


float* Mesh::returnPointNormal(Vertex* point)
{
	float normal[3] = { 0.0f, 0.0f, 0.0f };
	//loop through the triangles to trace each vertex, find normals, draw rays, calculate and add diameters to vertex's diameter attribute
	for (size_t t= 0; t< point->triList.size(); t++)
	{
		if (tris[t]->deleted)
			continue;
		int triangleIndex = tris[t]->tri_idx;
		//get the vertex index number of the vertices of the triangle at hand
		triVertsIds vertexIdsOfTriangle = getVertexIdsOfTriangleMesh(this, triangleIndex);

		int selectedVertexNumber = -1;
		//from the selected vertex of base triangle
		for (size_t v = 0; v < 3; v++)
		{
			if (point->idx == vertexIdsOfTriangle[v])
			{
				selectedVertexNumber = v;
				break;
			}
		}
		//get the coordinates of the vertices of the triangle at hand (by using the vertex index numbers)
		triVertsCoords coordinatesOfVerticesOfTriangle = getCoordsOfTriangleMesh(this, vertexIdsOfTriangle);
		//d is the normal vector from the selected vertex of the base triangle
		if (selectedVertexNumber == -1)
		{
			int a = 5;
			continue;
		}
		float d[3] = { 0.0f, 0.0f, 0.0f };
		calculateNormalVectorMesh(d, coordinatesOfVerticesOfTriangle, selectedVertexNumber);
		normal[0] += d[0];
		normal[1] += d[1];
		normal[2] += d[2];
	}
	float zero[3] = { 0.0f, 0.0f, 0.0f };
	float distance = VectorMath::distanceBetweenVectors(zero, normal);
	if (distance - 0.0f < 0.001)
		distance = 0.001;
	normal[0] /= distance;
	normal[1] /= distance;
	normal[2] /= distance;

	return normal;
}


triVertsIds getVertexIdsOfTriangleMesh(const Mesh* mesh, const int triangleIndex)
{
	triVertsIds vertexIdsOfTriangle;
	vertexIdsOfTriangle[0] = mesh->tris[triangleIndex]->v1i;
	vertexIdsOfTriangle[1] = mesh->tris[triangleIndex]->v2i;
	vertexIdsOfTriangle[2] = mesh->tris[triangleIndex]->v3i;
	return vertexIdsOfTriangle;
}

triVertsCoords getCoordsOfTriangleMesh(const Mesh* mesh, const triVertsIds vertexIdsOfTriangle)
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

triOtherVertsCoords getOtherCoordsOfTriangleMesh(const triVertsCoords& coordsOfVerticesOfTriangle, const int selectedVertexNumber)
{
	triOtherVertsCoords coordsOfOtherVertices;
		int number = 0;
		for (size_t otherVertexNumber = 0; otherVertexNumber < 3; otherVertexNumber++)
		{
			if (selectedVertexNumber == otherVertexNumber)
				continue;
			//fill the coordinate values of other 2 vertices' array
			for (size_t coordinate = 0; coordinate < 3; coordinate++) {
				coordsOfOtherVertices[number][coordinate] = coordsOfVerticesOfTriangle[otherVertexNumber][coordinate];
			}
			number++;
		}
	return coordsOfOtherVertices;
}

triOtherVertsCoords getVectorsToTheOtherVerticesMesh(const triOtherVertsCoords& coordsOfOtherVertices, const int selectedVertexNumber, const triVertsCoords& coordsOfVerticesOfTriangle) {
	triOtherVertsCoords vectorsToTheOtherVertices;
	for (size_t otherVertexNumber = 0; otherVertexNumber < 2; otherVertexNumber++)
	{
		for (size_t coordinate = 0; coordinate < 3; coordinate++) {
			vectorsToTheOtherVertices[otherVertexNumber][coordinate] =
				coordsOfOtherVertices[otherVertexNumber][coordinate]
				- coordsOfVerticesOfTriangle[selectedVertexNumber][coordinate];
		}
	}
	return vectorsToTheOtherVertices;
}

void Mesh::calculateTriangleNormal(Triangle* triangle)
{
	// Get the coordinates of the vertices of the triangle
	float* v1 = verts[triangle->v1i]->coords;
	float* v2 = verts[triangle->v2i]->coords;
	float* v3 = verts[triangle->v3i]->coords;

	// Calculate the vectors representing two edges of the triangle
	float edge1[3] = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
	float edge2[3] = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };

	// Calculate the cross product of the two edges to get the normal vector
	float normal[3] = {
		edge1[1] * edge2[2] - edge1[2] * edge2[1],
		edge1[2] * edge2[0] - edge1[0] * edge2[2],
		edge1[0] * edge2[1] - edge1[1] * edge2[0]
	};

	// Normalize the normal vector
	float length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	triangle->triangle_normal[0] = normal[0] / length;
	triangle->triangle_normal[1] = normal[1] / length;
	triangle->triangle_normal[2] = normal[2] / length;
}


void Mesh::calculateVertexNormals()
{
	// Iterate over each vertex
	for (Vertex* vertex : verts) {
		// Initialize the sum of normals to zero
		float sumNormals[3] = { 0.0f, 0.0f, 0.0f };

		// Iterate over each triangle in the vertex's triangle list
		for (int triIdx : vertex->triList) {
			// Get the triangle object
			Triangle* triangle = tris[triIdx];

			// Calculate the normal vector for the triangle
			calculateTriangleNormal(triangle);

			// Add the triangle's normal to the sum
			sumNormals[0] += triangle->triangle_normal[0];
			sumNormals[1] += triangle->triangle_normal[1];
			sumNormals[2] += triangle->triangle_normal[2];
		}

		// Normalize the sum of normals to get the vertex normal
		float length = sqrt(sumNormals[0] * sumNormals[0] + sumNormals[1] * sumNormals[1] + sumNormals[2] * sumNormals[2]);
		vertex->point_normal[0] = sumNormals[0] / length;
		vertex->point_normal[1] = sumNormals[1] / length;
		vertex->point_normal[2] = sumNormals[2] / length;
	}
}

void Mesh::calculateVertexNormal(Vertex* vertex)
{
	float sumNormals[3] = { 0.0f, 0.0f, 0.0f };

	// Iterate over each triangle in the vertex's triangle list
	for (int triIdx : vertex->triList) {
		// Get the triangle object
		Triangle* triangle = tris[triIdx];

		// Calculate the normal vector for the triangle
		calculateTriangleNormal(triangle);

		// Add the triangle's normal to the sum
		sumNormals[0] += triangle->triangle_normal[0];
		sumNormals[1] += triangle->triangle_normal[1];
		sumNormals[2] += triangle->triangle_normal[2];
	}

	// Normalize the sum of normals to get the vertex normal
	float length = sqrt(sumNormals[0] * sumNormals[0] + sumNormals[1] * sumNormals[1] + sumNormals[2] * sumNormals[2]);
	vertex->point_normal[0] = sumNormals[0] / length;
	vertex->point_normal[1] = sumNormals[1] / length;
	vertex->point_normal[2] = sumNormals[2] / length;
}

// Calculate the tangent plane for a triangle
void Mesh::calculateTriangleTangentPlane(Triangle* triangle) {
	// Choose one vertex of the triangle as the origin point of the tangent plane
	float* origin = verts[triangle->v1i]->coords;

	// Calculate the normal vector of the triangle
	calculateTriangleNormal(triangle);
	float* normal = triangle->triangle_normal;

	// Store the tangent plane data in the triangle
	std::copy(normal, normal + 3, triangle->tangent_plane_normal);
	std::copy(origin, origin + 3, triangle->tangent_plane_origin);
}

// Calculate the tangent plane for a vertex
void Mesh::calculateVertexTangentPlane(Vertex* vertex) {
	// Initialize sum of normals to calculate average
	float sumNormals[3] = { 0.0f, 0.0f, 0.0f };
	int numAdjacentTriangles = 0;

	// Iterate over each triangle adjacent to the vertex
	for (int triIdx : vertex->triList) {
		// Calculate the normal vector for the triangle
		Triangle* triangle = tris[triIdx];
		calculateTriangleNormal(triangle);
		float* normal = triangle->triangle_normal;

		// Add the normal vector to the sum
		sumNormals[0] += normal[0];
		sumNormals[1] += normal[1];
		sumNormals[2] += normal[2];

		numAdjacentTriangles++;
	}

	// Compute the average normal vector
	float avgNormal[3] = {
		sumNormals[0] / numAdjacentTriangles,
		sumNormals[1] / numAdjacentTriangles,
		sumNormals[2] / numAdjacentTriangles
	};

	// Use the vertex coordinates as the origin point of the tangent plane
	float* origin = vertex->coords;

	// Store the tangent plane data in the vertex
	std::copy(avgNormal, avgNormal + 3, vertex->tangent_plane_normal);
	std::copy(origin, origin + 3, vertex->tangent_plane_origin);
}

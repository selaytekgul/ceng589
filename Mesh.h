#pragma once

#include <iostream>
#include <vector>
#include <map>
#include "TypeDefinitions.h"
#include "VectorMath.h"

struct Vertex
{
	//float* coords, * normals; //3d coordinates etc
	float* coords, color[3] = {0.0f, 0.0f, 0.0f}, diameter = -5.0f; //3d coordinates etc
	int clusterId = 3;
	static float minDiameter;
	static float maxDiameter;
	bool isItInLongestBoundary = false;
	int idx; //who am i; verts[idx]
	std::vector< int > vertList; //adj vertices;
	std::vector< int > triList;
	std::vector< int > edgeList;
	Vertex(int i, float* c) : idx(i), coords(c) {};
};

struct Edge
{
	int edge_idx; //edges[idx]
	int v1i, v2i; //endpnts
	float length;
	int existedTriangeNumber = 0;
	bool isItBoundary = false;
	bool isItPathPart = false;
	bool isItInLongestBoundary = false;
	bool isItTraversed = false;
	bool isInShortestPath = false;
	std::vector< int > triList;
	Edge(int id, int v1, int v2) : edge_idx(id), v1i(v1), v2i(v2) { computeLength(); };

	void computeLength()
	{
		length = 1;
	}
};

struct Triangle
{
	int tri_idx; //tris[idx]
	int v1i, v2i, v3i;
	Triangle(int id, int v1, int v2, int v3) : tri_idx(id), v1i(v1), v2i(v2), v3i(v3) {};
};

class Mesh
{
private:
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2, int triIdx);
	void modifyEdge(int v1, int v2, int triIdx);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i, int triIdx);
public:
	void setMinMaxDiameters();
	void discardInfAndNegativeDiameters();
	std::vector< Vertex* > verts;
	std::vector< Triangle* > tris;
	std::vector< Edge* > edges;
	std::vector< int > samples;
	std::map<int, int> boundIndexToMeshId;
	std::map<int, int> meshIndexToBoundId;

	Mesh() {} ;
	void createCube(float side);
	void createOpenCube(float side);
	void createDoubleOpenCube(float side);
	void loadOff(char* name);
};

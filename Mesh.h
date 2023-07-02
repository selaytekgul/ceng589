#pragma once

#include <iostream>
#include <vector>
#include <array>
#include "TypeDefinitions.h"

struct Vertex
{
	//float* coords, * normals; //3d coordinates etc
	float* coords, color[3] = {0.0f, 0.0f, 0.0f}, length = -5.0f; //3d coordinates etc
	static float minLength;
	static float maxLength;
	int idx; //who am i; verts[idx]
	std::vector< int > vertList; //adj vertices;
	std::vector< int > triList;
	std::vector< int > edgeList;
	std::vector< float* > normalList;
	bool hasLongLength = false;

	Vertex(int i, float* c) : idx(i), coords(c) {};
};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	float length;
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2) { computeLength(); };

	void computeLength()
	{
		length = 7;
	}
};

struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {};
};

class Mesh
{
private:
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i);
public:
	void setMinMaxLenghts();
	void discardInfAndNegativeLenghts();
	std::vector< Vertex* > verts;
	std::vector< Triangle* > tris;
	std::vector< Edge* > edges;


	Mesh() {} ;
	void createCube(float side);
	void loadOff(char* name);
};

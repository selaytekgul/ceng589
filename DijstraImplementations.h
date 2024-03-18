#pragma once
#include "AdjListGraph.h"
#include "Mesh.h"
#include "VectorMath.h"
//***C++11 Style:***
//https://stackoverflow.com/a/27739925
#include <chrono>

namespace Dijstra
{
	Graph meshToGraph(const Mesh* mesh)
	{
		const int V = mesh->verts.size();
		Graph gmesh(V);

		for (const Edge* edge : mesh->edges)
		{
			const int v1Idx = edge->v1i;
			const int v2Idx = edge->v2i;
			Vertex* vert1 = mesh->verts[v1Idx];
			Vertex* vert2 = mesh->verts[v2Idx];
			const float* coords1 = vert1->coords;
			const float* coords2 = vert2->coords;
			const float weight = VectorMath::distanceBetweenVectors(coords1, coords2);
			gmesh.addEdge(v1Idx, v2Idx, weight);
		}
		return gmesh;
	}
	
	void fprinting(Graph g)
	{
		const int numNodes = g.V;
		for (size_t i = 0; i < numNodes; i++)
		{
			g.shortestPath(i, true, -1);
		}
	}

	void timing(Graph g, int source, int dest)
	{
		std::cout << "Time is started." << std::endl;
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		g.shortestPath(source, false, dest);
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
	}

	void pathDrawing(Mesh* mesh, Graph g, int source, int dest)
	{
		std::vector<int> path = g.shortestPath(source, false, dest);
		for (Edge* edge : mesh->edges)
		{
			const int v1Idx = edge->v1i;
			const int v2Idx = edge->v2i;
			auto iter1 = std::find(path.begin(), path.end(), v1Idx);
			auto iter2 = std::find(path.begin(), path.end(), v2Idx);
			if (iter1 != path.end() || iter2 != path.end())
				edge->isItTraversed = true;
		}
		bool sourceIsFound = false;
		int index = dest;
		while (!sourceIsFound)
		{
			int wasAt = path[index];
			for (const auto neighEdge : mesh->verts[index]->edgeList)
			{
				if (mesh->edges[neighEdge]->v1i == index && mesh->edges[neighEdge]->v2i == wasAt
					|| mesh->edges[neighEdge]->v2i == index && mesh->edges[neighEdge]->v1i == wasAt)
				{
					mesh->edges[neighEdge]->isInShortestPath = true;
				}
			}
			index = wasAt;
			if (index == source)
				sourceIsFound = true;
		}
	}
}


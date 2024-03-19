#pragma once
#include "AdjListGraph.h"
#include "Mesh.h"
#include "VectorMath.h"
//***C++11 Style:***
//https://stackoverflow.com/a/27739925
#include <chrono>

namespace Dijkstra
{
	enum MinFindMethod
	{
		MIN_HEAP,
		ARRAY,
		FIB_HEAP
	};
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
	
	void fprinting(Graph g, MinFindMethod method)
	{
		const int numNodes = g.V;
		for (size_t i = 0; i < numNodes; i++)
		{
			switch (method)
			{
			case Dijkstra::MIN_HEAP:
				g.shortestPath(i, true, -1);
				break;
			case Dijkstra::ARRAY:
				g.shortestPathArray(i, true, -1);
				break;
			case Dijkstra::FIB_HEAP:
				break;
			default:
				//g.shortestPathFibHeap(i, true, -1);
				break;
			}
		}
	}
		
	void fprintingOnce(Graph g, int source, int dest, MinFindMethod method)
	{
		switch (method)
		{
		case Dijkstra::MIN_HEAP:
			g.shortestPath(source, true, dest);
			break;
		case Dijkstra::ARRAY:
			//g.shortestPathArray(source, true, dest);
			break;
		case Dijkstra::FIB_HEAP:
			break;
		default:
			//g.shortestPathFibHeap(source, true, dest);
			break;
		}
	}

	void timing(Graph g, int source, int dest, MinFindMethod method)
	{
		std::chrono::steady_clock::time_point begin;
		std::chrono::steady_clock::time_point end;
		switch (method)
		{
		case Dijkstra::MIN_HEAP:
			std::cout << "Time is started." << std::endl;
			begin = std::chrono::steady_clock::now();
			g.shortestPath(source, false, dest);
			end = std::chrono::steady_clock::now();
			break;
		case Dijkstra::ARRAY:
			std::cout << "Time is started." << std::endl;
			begin = std::chrono::steady_clock::now();
			g.shortestPathArray(source, false, dest);
			end = std::chrono::steady_clock::now();			
			break;
		case Dijkstra::FIB_HEAP:
			std::cout << "Time is started." << std::endl;
			begin = std::chrono::steady_clock::now();
			//g.shortestPathFibHeap(source, false, dest);
			end = std::chrono::steady_clock::now();			
			break;
		default:
			break;
		}
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
	}

	void pathDrawing(Mesh* mesh, Graph g, int source, int dest, MinFindMethod method)
	{
		std::vector<int> path = {};
		switch (method)
		{
		case Dijkstra::MIN_HEAP:
			path = g.shortestPath(source, false, dest);
			break;
		case Dijkstra::ARRAY:
			path = g.shortestPathArray(source, false, dest);
			break;
		case Dijkstra::FIB_HEAP:
			break;
		default:
			//path = g.shortestPathFibHeap(source, false, dest);
			break;
		}
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


#pragma once

#include <utility>
#include <list>
#include <vector>
#include <queue>
#include "FibHeap.hpp"
// C++ Program to find Dijkstra's shortest path using
// https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
// priority_queue in STL
using namespace std;
#define INF 0x3f3f3f3f

// iPair ==> Integer Pair
typedef pair<int, float> iPair;

// This class represents a directed graph using
// adjacency list representation
class Graph {

	// In a weighted graph, we need to store vertex
	// and weight pair for every edge
	list<pair<int, float> >* adj;

public:
	int V; // No. of vertices
	Graph(int V); // Constructor

	// function to add an edge to graph
	void addEdge(int u, int v, float w);

	// prints shortest path from s
	vector<int> Graph::shortestPath(int src, bool printOpen, int dest);
	vector<int> Graph::shortestPathFibHeap(int src, bool printOpen, int dest);
	vector<int> Graph::shortestPathArray(int src, bool printOpen, int dest);

};

// Allocates memory for adjacency list
Graph::Graph(int V)
{
	this->V = V;
	adj = new list<iPair>[V];
}

void Graph::addEdge(int u, int v, float w)
{
	adj[u].push_back(make_pair(v, w));
	adj[v].push_back(make_pair(u, w));
}

// Prints shortest paths from src to all other vertices
vector<int> Graph::shortestPath(int src, bool printOpen, int dest)
{
	// Create a priority queue to store vertices that
	// are being preprocessed. This is weird syntax in C++.
	// Refer below link for details of this syntax
	// https://www.geeksforgeeks.org/implement-min-heap-using-stl/
	priority_queue<iPair, vector<iPair>, greater<iPair> >
		pq;

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<float> dist(V, INF);
	vector<int> arrivedFrom(V, -1);

	// Insert source itself in priority queue and initialize
	// its distance as 0.
	pq.push(make_pair(0, src));
	dist[src] = 0;
	arrivedFrom[src] = -5;
	/* Looping till priority queue becomes empty (or all
	distances are not finalized) */
	while (!pq.empty()) {
		// The first vertex in pair is the minimum distance
		// vertex, extract it from priority queue.
		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		int u = pq.top().second;
		// mark the node u as visited
		// early termination possible here in the query timing
		if (u == dest && !printOpen)
		{
			//printf("dest is arrived.");
			return arrivedFrom;
		}
		pq.pop();

		// 'i' is used to get all adjacent vertices of a
		// vertex
		list<pair<int, float> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i) {
			// Get vertex label and weight of current
			// adjacent of u.
			int v = (*i).first;
			int weight = (*i).second;

			// If there is shorter path to v through u.
			if (dist[v] > dist[u] + weight) {
				// Updating distance of v
				dist[v] = dist[u] + weight;
				arrivedFrom[v] = u;
				pq.push(make_pair(dist[v], v));
			}
		}
	}
}

// Prints shortest paths from src to all other vertices
vector<int> Graph::shortestPathFibHeap(int src, bool printOpen, int dest)
{
	// Create a priority queue to store vertices that
	// are being preprocessed. This is weird syntax in C++.
	// Refer below link for details of this syntax
	// https://www.geeksforgeeks.org/implement-min-heap-using-stl/
	//priority_queue<iPair, vector<iPair>, greater<iPair> >
	//	pq;
	FH::FibHeap<int, float> h{};

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<float> dist(V, INF);
	vector<int> arrivedFrom(V, -1);

	// Insert source itself in priority queue and initialize
	// its distance as 0.
	//pq.push(make_pair(0, src));
	h.insert(src, 0);
	dist[src] = 0;
	arrivedFrom[src] = -5;
	/* Looping till priority queue becomes empty (or all
	distances are not finalized) */
	while (!h.getSize() == 0) {
		// The first vertex in pair is the minimum distance
		// vertex, extract it from priority queue.
		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		//int u = pq.top().second;
		int u = h.peekMax().value;
		// mark the node u as visited
		// early termination possible here in the query timing
		if (u == dest && !printOpen)
		{
			//printf("dest is arrived.");
			return arrivedFrom;
		}
		//pq.pop();
		h.popMax();
		// 'i' is used to get all adjacent vertices of a
		// vertex
		list<pair<int, float> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i) {
			// Get vertex label and weight of current
			// adjacent of u.
			int v = (*i).first;
			int weight = (*i).second;

			// If there is shorter path to v through u.
			if (dist[v] > dist[u] + weight) {
				// Updating distance of v
				dist[v] = dist[u] + weight;
				arrivedFrom[v] = u;
				//pq.push(make_pair(dist[v], v));
				h.insert(v, dist[v]);
			}
		}
	}

	// Print shortest distances stored in dist[]
	if (printOpen)
	{
		printf("Vertex Distance from Source %d\n", src);
		for (int i = 0; i < V; ++i)
			printf("%d \t\t %f \tfrom %d\n", i, dist[i], arrivedFrom[i]);
	}

	return arrivedFrom;
}

int returnMinWeightIndexInAnArray(const std::vector<float> weights, const std::vector<bool> visited)
{
	float minWeight = weights[0];
	int minIndex = 0;
	for (size_t i = 1; i < weights.size(); i++)
	{
		if (minWeight > weights[i] && visited[i] == false)
		{
			minWeight = weights[i];
			minIndex = i;
		}
	}
	return minIndex;
}

// Prints shortest paths from src to all other vertices
vector<int> Graph::shortestPathArray(int src, bool printOpen, int dest)
{
	// Create a priority queue to store vertices that
	// are being preprocessed. This is weird syntax in C++.
	// Refer below link for details of this syntax
	// https://www.geeksforgeeks.org/implement-min-heap-using-stl/
	//priority_queue<iPair, vector<iPair>, greater<iPair> >
	//	pq;
	//FH::FibHeap<int, float> h{};
	std::vector<float> weightsArr(V, std::numeric_limits<float>::max());
	std::vector<bool> visited(V, false);
	int numVisited = 0;

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<float> dist(V, INF);
	vector<int> arrivedFrom(V, -1);

	// Insert source itself in priority queue and initialize
	// its distance as 0.
	//pq.push(make_pair(0, src));
	//h.insert(src, 0);
	weightsArr[src] = 0;
	dist[src] = 0;
	arrivedFrom[src] = -5;
	/* Looping till priority queue becomes empty (or all
	distances are not finalized) */
	while (numVisited == V) {
		// The first vertex in pair is the minimum distance
		// vertex, extract it from priority queue.
		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		//int u = pq.top().second;
		//int u = h.peekMax().value;
		int u = returnMinWeightIndexInAnArray(weightsArr, visited);
		// mark the node u as visited
		// early termination possible here in the query timing
		if (u == dest && !printOpen)
		{
			//printf("dest is arrived.");
			return arrivedFrom;
		}
		//pq.pop();
		//h.popMax();
		visited[u] = false;
		// 'i' is used to get all adjacent vertices of a
		// vertex
		list<pair<int, float> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i) {
			// Get vertex label and weight of current
			// adjacent of u.
			int v = (*i).first;
			int weight = (*i).second;

			// If there is shorter path to v through u.
			if (dist[v] > dist[u] + weight) {
				// Updating distance of v
				dist[v] = dist[u] + weight;
				arrivedFrom[v] = u;
				//pq.push(make_pair(dist[v], v));
				//h.insert(v, dist[v]);
				weightsArr[v] = dist[v];
			}
		}
	}

	// Print shortest distances stored in dist[]
	if (printOpen)
	{
		printf("Vertex Distance from Source %d\n", src);
		for (int i = 0; i < V; ++i)
			printf("%d \t\t %f \tfrom %d\n", i, dist[i], arrivedFrom[i]);
	}

	return arrivedFrom;
}


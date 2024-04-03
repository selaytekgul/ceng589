#define HAVE_INT8_T
#include "GraphOperations.h"
#include <random>
#include <string>
#include <queue>

#include <Eigen/Dense>

int main(int, char ** argv)
{
	Mesh* mesh = new Mesh();
	Mesh* original_mesh = new Mesh();
	std::string fileName = { "equal_cube3" };
	mesh->loadOff(fileName + ".off");
	original_mesh->loadOff(fileName + ".off");

	// Create a priority queue of pairs of integers (max heap)
	std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> minHeap;

	for (size_t i = 0; i < mesh->edges.size(); i++)
	{
		// Insert key-value pairs into the max heap
		const int edgeidx = mesh->edges[i]->edge_idx;
		mesh->computeLength(edgeidx);
		minHeap.push({ mesh->edges[i]->length, edgeidx});
	}


	//CHANGE
	mesh->toOFF(fileName + "_.off");
	return 0;
}


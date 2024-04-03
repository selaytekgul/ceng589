#define HAVE_INT8_T
#include "GraphOperations.h"
#include <random>
#include <string>
#include <queue>

#include <Eigen/Dense>

void windingNumberTest(Mesh* mesh)
{
	float ptr3[3] = { -5.8959f, -16.9698f, 128.8068f };
	Vertex v3 = { (int)mesh->verts.size(), ptr3 };
	mesh->windingNumberByYusufSahillioglu(&v3);
	mesh->windingNumberByYusufSahillioglu(mesh->verts[0]);
	mesh->windingNumberByYusufSahillioglu(mesh->verts[1]);
	for (size_t i = 0; i < mesh->verts.size(); i++)
	{
		mesh->windingNumberByYusufSahillioglu(mesh->verts[i]);
		if (mesh->verts[i]->winding == 1)
		{
			std::cout << i << ": " << mesh->verts[i]->winding << std::endl;
		}
	}
}

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

	// //Print the elements of the max heap
	//std::cout << "Max heap elements (key-value pairs):" << std::endl;
	//while (!minHeap.empty() && mesh->numDeletedTri < mesh->verts.size()/2.0) {
	//	auto kvp = minHeap.top(); // Get the top element
	//	mesh->computeLength(kvp.second);
	//	if (kvp.first - mesh->edges[kvp.second]->length > 0.0001)
	//	{
	//		std::cout << "Skipped Key: " << kvp.first << ", Value: " << kvp.second << std::endl;
	//		minHeap.pop(); // Remove the top element
	//		continue;
	//	}
	//	std::cout << "Key: " << kvp.first << ", Value: " << kvp.second << std::endl;
	//	//mesh->collapseEdgeTo(mesh->edges[kvp.first], mesh->edges[kvp.first]->v1i);
	//	mesh->collapseEdge(mesh->edges[kvp.second], &minHeap);
	//	minHeap.pop(); // Remove the top element
	//}

	for (size_t i = 0; i < mesh->verts.size(); i++)
	{
		if (mesh->verts[i]->deleted)
			continue;
		original_mesh->inflatePoint(mesh->verts[i]);
	}

	mesh->toOFF(fileName + "inflate.off");
	return 0;
}


#define HAVE_INT8_T
#include "Mesh.h"
#include <string>
#include <queue>

int main(int, char ** argv)
{
	Mesh* mesh = new Mesh();
	Mesh* original_mesh = new Mesh();
	std::string fileName = { "cube_24" };
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



	////Print the elements of the max heap
	//std::cout << "Max heap elements (key-value pairs):" << std::endl;
	////while (!minHeap.empty() && mesh->numDeletedTri < mesh->verts.size()/4.0) {
	//while (!minHeap.empty() && mesh->numDeletedTri < 1) {
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

	//for (size_t i = 0; i < mesh->verts.size(); i++)
	//{
	//	if (mesh->verts[i]->deleted)
	//		continue;
	//	original_mesh->inflatePoint(mesh->verts[i]);
	//}
	//int i = 3;
	//original_mesh->inflatePoint(mesh->verts[i]);

	//CHANGE
	//mesh->toOFF(fileName + "_" + std::to_string(i) + "_inflate.off");
	mesh->toOFF(fileName + "_all_inflate_05.off");
	//mesh->toOFF(fileName + "_all" + "_collapse_inflate_original_mesh_winding.off");
	//mesh->toOFF(fileName + "_all_collapse_original_mesh_winding.off");
	return 0;
}


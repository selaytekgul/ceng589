#define HAVE_INT8_T
#include "Segmentor.h"
#include "KMeans.h"
#include "GraphOperations.h"
#include "DijstraImplementations.h"
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

	//std::string fileName = { "tr_reg_000" };
	//
	//std::string fileName = { "for timing/centaur" };
	//std::string fileName = { "for timing/man" };
	//
	//std::string fileName = { "for timing/weirdSphere" };
	//std::string fileName = { "for fprinting/horse0" };
	//std::string fileName = { "for fprinting/man0" };

	//std::string fileName = { "0" };
	//std::string fileName = { "1" };
	//std::string fileName = { "car" };
	//std::string fileName = { "coffeecup" };
	//std::string fileName = { "bunny" };
	std::string fileName = { "cube3" };
	//std::string fileName = { "equal_cube3" };
	
	//std::string fileName = { "faces/face" };
	//std::string fileName = { "faces/face-low" };
	//std::string fileName = { "faces/facem" };
	//std::string fileName = { "faces/facem-low" };
	
	//std::string fileName = { "doubleOpenCube3" };
	//std::string fileName = { "man0" };

	//mesh->createCube(20.0f);
	//mesh->createOpenCube(20.0f);
	//mesh->createDoubleOpenCube(20.0f);

	mesh->loadOff(fileName + ".off");
	original_mesh->loadOff(fileName + ".off");

	//int ind = 0;
	//while (mesh->numDeletedTri < 4)
	//{
	//	mesh->collapseEdgeTo(mesh->edges[ind], mesh->edges[ind++]->v1i);
	//}

	//for (size_t i = 0; i < mesh->verts.size(); i++)
	//{
	//	if (mesh->verts[i]->deleted)
	//		continue;
	//	original_mesh->inflatePoint(mesh->verts[i]);
	//}

	// Create a priority queue of pairs of integers (max heap)
	//std::priority_queue<std::pair<int, int>> maxHeap;

	//for (size_t i = 0; i < mesh->edges.size(); i++)
	//{
	//	// Insert key-value pairs into the max heap
	//	const int edgeidx = mesh->edges[i]->edge_idx;
	//	mesh->computeLength(edgeidx);
	//	maxHeap.push({ edgeidx , mesh->edges[i]->length});
	//}

	// Print the elements of the max heap
	//std::cout << "Max heap elements (key-value pairs):" << std::endl;
	//while (!maxHeap.empty() && mesh->numDeletedTri < 6) {
	//	auto kvp = maxHeap.top(); // Get the top element
	//	std::cout << "Key: " << kvp.first << ", Value: " << kvp.second << std::endl;
	//	mesh->collapseEdgeTo(mesh->edges[kvp.first], mesh->edges[kvp.first]->v1i);
	//	maxHeap.pop(); // Remove the top element
	//}
	mesh->collapseEdge(mesh->edges[0]);

	mesh->toOFF(fileName + "_collapse.off");
	return 0;
}


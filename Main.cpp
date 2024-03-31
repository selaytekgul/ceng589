#define HAVE_INT8_T
#include "Segmentor.h"
#include "KMeans.h"
#include "GraphOperations.h"
#include "DijstraImplementations.h"
#include <random>
#include <string>

#include <Eigen/Dense>

int main(int, char ** argv)
{
	Mesh* mesh = new Mesh();

	//std::string fileName = { "tr_reg_000" };
	//
	//std::string fileName = { "for timing/centaur" };
	//std::string fileName = { "for timing/man" };
	//
	//std::string fileName = { "for timing/weirdSphere" };
	//std::string fileName = { "for fprinting/horse0" };
	//std::string fileName = { "for fprinting/man0" };

	//std::string fileName = { "0" };
	//std::string fileName = { "car" };
	//std::string fileName = { "coffeecup" };
	//std::string fileName = { "bunny" };
	//std::string fileName = { "cube3" };
	
	//std::string fileName = { "faces/face" };
	//std::string fileName = { "faces/face-low" };
	//std::string fileName = { "faces/facem" };
	//std::string fileName = { "faces/facem-low" };
	
	//std::string fileName = { "doubleOpenCube3" };
	std::string fileName = { "man0" };

	//mesh->createCube(20.0f);
	//mesh->createOpenCube(20.0f);
	//mesh->createDoubleOpenCube(20.0f);

	mesh->loadOff(fileName + ".off");

	mesh->windingNumberByYusufSahillioglu(mesh->verts[0]);
	mesh->windingNumberByYusufSahillioglu(mesh->verts[1]);
	for (size_t i = 0; i < mesh->verts.size(); i++)
	{
		mesh->windingNumberByYusufSahillioglu(mesh->verts[i]);
		if (mesh->verts[i]->winding == 1)
		{
			int c = 7;
			std::cout << i << ": " << mesh->verts[i]->winding << std::endl;
		}
	}
	int a = 5;
	//float ptr[3] = {-15.f, 10.f, -5.f };
	//float ptr1[3] = {-45.f, 10.f, -5.f };
	//float ptr1[3] = {45.f, 10.f, -5.f };
	float ptr3[3] = {-5.8959f, -16.9698f, 128.8068f };
	//Vertex v = { (int)mesh->verts.size(), ptr };
	//mesh->windingNumberByYusufSahillioglu(&v);
	//Vertex v1 = { (int)mesh->verts.size(), ptr1 };
	//mesh->windingNumberByYusufSahillioglu(&v1);
	Vertex v3 = { (int)mesh->verts.size(), ptr3 };
	mesh->windingNumberByYusufSahillioglu(&v3);

	int b = 9;
	return 0;
}

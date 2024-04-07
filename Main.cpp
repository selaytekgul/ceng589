#include <cmath> // For mathematical functions like sqrt
#include <algorithm> // For algorithms like std::copy
#include <iostream> // For input/output operations
#include <string> // For string manipulation
#include <queue> // For priority_queue
#include <vector> // For std::vector
#include <fstream> // For file input/output operations

// Include your Mesh.h header file
#include "Mesh.h"


int main(int, char** argv) {
    Mesh* mesh = new Mesh();
    Mesh* original_mesh = new Mesh();
    std::string fileName = { "381" };
    mesh->loadOff(fileName + ".off");
    original_mesh->loadOff(fileName + ".off");

    // Calculate tangent planes for triangles
    for (size_t i = 0; i < mesh->tris.size(); i++) {
        mesh->calculateTriangleTangentPlane(mesh->tris[i]);
    }

    // Calculate tangent planes for vertices
    for (size_t i = 0; i < mesh->verts.size(); i++) {
        mesh->calculateVertexTangentPlane(mesh->verts[i]);
    }

    // Create a priority queue of pairs of floats (min heap)
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> minHeap;

    for (size_t i = 0; i < mesh->edges.size(); i++) {
        // Insert key-value pairs into the min heap
        const int edgeidx = mesh->edges[i]->edge_idx;
        mesh->computeLength(edgeidx);
        minHeap.push({ mesh->edges[i]->length, edgeidx });
    }

    // Print elements of the min heap
    std::cout << "Min heap elements (key-value pairs):" << std::endl;
    while (!minHeap.empty() && mesh->numDeletedTri < mesh->verts.size()/2) {
    //while (!minHeap.empty() && mesh->numDeletedTri < 50) {
        auto kvp = minHeap.top(); // Get the top element
        mesh->computeLength(kvp.second);
        if (kvp.first - mesh->edges[kvp.second]->length > 0.0001) {
            std::cout << "Skipped key: " << kvp.first << ", value: " << kvp.second << std::endl;
            minHeap.pop(); // Remove the top element
            continue;
        }
        std::cout << "Key: " << kvp.first << ", value: " << kvp.second << std::endl;
        mesh->collapseEdge(mesh->edges[kvp.second], &minHeap);
        minHeap.pop(); // Remove the top element
    }
    mesh->toOFF(fileName + "_+.off"); // Mesh after collapsing edges

    // Inflate points after collapsing edges
    for (size_t i = 0; i < mesh->verts.size(); i++) {
        if (mesh->verts[i]->deleted)
            continue;
        original_mesh->inflatePoint(mesh->verts[i]);
    }

    // Save the modified mesh to OFF files
    mesh->toOFF(fileName + "_all_inflate_05.off"); // Mesh after inflating points
    //mesh->toOFF(fileName + "_all" + "_collapse_inflate_original_mesh_winding.off");
    //mesh->toOFF(fileName + "_all_collapse_original_mesh_winding.off");

    return 0;
}

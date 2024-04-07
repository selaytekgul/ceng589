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

    //// Calculate tangent planes for triangles
    //for (size_t i = 0; i < mesh->tris.size(); i++) {
    //    mesh->calculateTriangleTangentPlane(mesh->tris[i]);
    //    original_mesh->calculateTriangleTangentPlane(mesh->tris[i]);
    //}

    //// Calculate tangent planes for vertices
    //for (size_t i = 0; i < mesh->verts.size(); i++) {
    //    mesh->calculateVertexTangentPlane(mesh->verts[i]);
    //    original_mesh->calculateVertexTangentPlane(mesh->verts[i]);
    //}

    // Create a priority queue of pairs of floats (min heap)
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> minHeap;

    //// Calculate distance_between tangent planes and edge midpoints
    //for (size_t i = 0; i < mesh->edges.size(); i++) {
    //    mesh->computeDistFromEdgeMidToEndPntsTangPla(i);
    //    original_mesh->computeDistFromEdgeMidToEndPntsTangPla(i);
    //    minHeap.push({ mesh->edges[i]->midToEndPointTangentPlanesDist, mesh->edges[i]->edge_idx });
    //}

    for (size_t i = 0; i < mesh->edges.size(); i++) {
        // Insert key-value pairs into the min heap
        const int edgeidx = mesh->edges[i]->edge_idx;
        mesh->computeLength(edgeidx);
        minHeap.push({ mesh->edges[i]->length, edgeidx });
    }

    int numSkippedKey = 0;
    //while (!minHeap.empty() && mesh->numDeletedTri < mesh->tris.size()/2 && numSkippedKey < mesh->tris.size()/2) {
    //while (!minHeap.empty() && mesh->numDeletedTri < 70 * mesh->tris.size()/100 && numSkippedKey < mesh->tris.size() / 2) {
    //while (!minHeap.empty() && mesh->numDeletedTri < 80 * mesh->tris.size()/100 && numSkippedKey < mesh->tris.size()/2) {
    //while (!minHeap.empty() && mesh->numDeletedTri < 90 * mesh->tris.size()/100 && numSkippedKey < mesh->tris.size()/2) {
    while (!minHeap.empty() && mesh->numDeletedTri < 95 * mesh->tris.size()/100 && numSkippedKey < mesh->tris.size()/3) {
        
        auto kvp = minHeap.top(); // Get the top element
        
        mesh->computeLength(kvp.second);
        //mesh->computeDistFromEdgeMidToEndPntsTangPla(kvp.second);

        if (kvp.first - mesh->edges[kvp.second]->length > 0.0001) {
        //if (kvp.first - mesh->edges[kvp.second]->midToEndPointTangentPlanesDist > 0.0001) {
            std::cout << "Skipped key: " << kvp.first << ", value: " << kvp.second << std::endl;
            minHeap.pop(); // Remove the top element
            numSkippedKey++;
            continue;
        }
        numSkippedKey = 0;
        std::cout << "Key: " << kvp.first << ", value: " << kvp.second << std::endl;
        if (mesh->edges[kvp.second]->deleted)
        {
            minHeap.pop(); // Remove the top element
            continue;
        }
        bool collapsed = mesh->collapseEdge(mesh->edges[kvp.second], &minHeap);
        
        if (!collapsed)
        {
            int a = 9;
        }
        minHeap.pop(); // Remove the top element
    }


    mesh->toOFF("5l.off"); // Mesh after collapsing edges

     //Inflate points after collapsing edges
    for (size_t i = 0; i < mesh->verts.size(); i++) {
        if (mesh->verts[i]->deleted)
            continue;
        mesh->inflatePoint(original_mesh, mesh->verts[i]);
    }
    //original_mesh->inflatePoint(mesh->verts[0]);

    // Save the modified mesh to OFF files
    mesh->toOFF("5li.off");

    return 0;
}

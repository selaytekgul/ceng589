#pragma once

#include <vector>
#include <array>
#include <set>
#include <utility>
#include <cmath>
#include "VectorMath.h"
#include "Mesh.h"
#include "TriangleMeshMath.h"

#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>
#include <map>
namespace GraphOperations
{
    enum ParameterizationMethod
    {
        UNIFORM = 0,
        HARMONIC = 1,
        MEAN = 2
    };

    inline std::pair<float*, float> findClosestVertex(const float* source, std::vector<float*> targetList);
    inline void generateDiskParameterization(Mesh* mesh, const ParameterizationMethod method);
    inline std::set<Edge*> findAnyBoundaries(const Mesh* mesh);
    inline std::vector<std::vector<int>> adjacencyMatrixFromEdges(Mesh* mesh, const std::set<Edge*>& boundaryEdges);
    inline std::vector<int> longestCycleConnectedComponent(const std::vector<std::vector<int>>& graph);
    inline void longestCycleConnectedComponentSetToEdgeInfo(Mesh* mesh, const std::vector<int>& longestBoundary);



    std::pair<float*, float> findClosestVertex(const float* source, std::vector<float*> targetList)
    {
        float* closest = NULL;
        float minDistance = FLT_MAX;

        //targets.forEach((target) = > {
        for (size_t i = 0; i < targetList.size(); i++)
        {
            const float currentDistance = VectorMath::distanceBetweenVectors(source, targetList[i]);
            if (currentDistance < minDistance) {
                closest = targetList[i];
                minDistance = currentDistance;
            }
        }
        return { closest, minDistance };
    }

    void generateDiskParameterization(Mesh* mesh, const ParameterizationMethod method)
    {
        std::set<Edge*> boundaryEdges = findAnyBoundaries(mesh);
        std::vector<std::vector<int>> graph = adjacencyMatrixFromEdges(mesh, boundaryEdges);
        std::vector<int> longestBoundary = longestCycleConnectedComponent(graph);
        longestCycleConnectedComponentSetToEdgeInfo(mesh, longestBoundary);
        int x = 0;
    }

    std::set<Edge*> findAnyBoundaries(const Mesh* mesh)
    {
        std::set<Edge*> boundaryEdges;
        for (size_t i = 0; i < mesh->edges.size(); i++)
        {
            if (mesh->edges[i]->existedTriangeNumber == 1)
            {
                boundaryEdges.insert(mesh->edges[i]);
                mesh->edges[i]->isItBoundary = true;
            }
        }
        return boundaryEdges;
    }

    std::vector<std::vector<int>> adjacencyMatrixFromEdges(Mesh* mesh, const std::set<Edge*>& boundaryEdges) {
        int index = 0;
        for (auto boundEdge : boundaryEdges)
        {
            if (mesh->meshIndexToBoundId.find(boundEdge->v1i) == mesh->meshIndexToBoundId.end())
            {
                mesh->boundIndexToMeshId[index] = boundEdge->v1i;
                mesh->meshIndexToBoundId[boundEdge->v1i] = index;
                index++;
            }
            
            if (mesh->meshIndexToBoundId.find(boundEdge->v2i) == mesh->meshIndexToBoundId.end())
            {
                mesh->boundIndexToMeshId[index] = boundEdge->v2i;
                mesh->meshIndexToBoundId[boundEdge->v2i] = index;
                index++;
            }
        }

        int numBoundEdges = boundaryEdges.size();
        std::vector<std::vector<int>> adjMatrix(numBoundEdges, std::vector<int>(numBoundEdges, false));
        for (auto boundEdge : boundaryEdges)
        {
            adjMatrix[mesh->meshIndexToBoundId[boundEdge->v1i]][mesh->meshIndexToBoundId[boundEdge->v2i]] =  true;
            adjMatrix[mesh->meshIndexToBoundId[boundEdge->v2i]][mesh->meshIndexToBoundId[boundEdge->v1i]] =  true;
        }
        return adjMatrix;
    }



    void dfs(const std::vector<std::vector<int>>& graph, int node, std::vector<bool>& visited, std::vector<int>& parent, std::vector<int>& cycle, std::vector<bool>& in_cycle) {
        visited[node] = true;
        for (int neighbor = 0; neighbor < graph.size(); ++neighbor) {
            if (graph[node][neighbor]) {
                if (!visited[neighbor]) {
                    parent[neighbor] = node;
                    dfs(graph, neighbor, visited, parent, cycle, in_cycle);
                }
                else if (!in_cycle[neighbor] && neighbor != parent[node]) {
                    int current = node;
                    while (current != neighbor) {
                        cycle.push_back(current);
                        current = parent[current];
                    }
                    cycle.push_back(neighbor);
                    in_cycle[neighbor] = true;
                    in_cycle[node] = true;
                    if (cycle.size() > 2 && (cycle.front() == cycle.back() || cycle[1] == cycle.back())) // To ensure that the cycle is closed
                        return;
                }
            }
        }
    }


    std::vector<int> longestCycleConnectedComponent(const std::vector<std::vector<int>>& graph) {
        int n = graph.size();
        std::vector<bool> visited(n, false);
        std::vector<int> parent(n, -1);
        std::vector<int> cycle;
        std::vector<bool> in_cycle(n, false);

        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                dfs(graph, i, visited, parent, cycle, in_cycle);
            }
        }

        return cycle;
    }
    
    void longestCycleConnectedComponentSetToEdgeInfo(Mesh* mesh, const std::vector<int>& longestBoundary) {
        for (size_t i = 0; i < longestBoundary.size(); i++)
        {
            int vertIdx1;
            int vertIdx2;
            if (i == longestBoundary.size() - 1)
            {
                vertIdx1 = longestBoundary[i];
                vertIdx2 = longestBoundary[0];
            }
            else
            {
                vertIdx1 = longestBoundary[i];
                vertIdx2 = longestBoundary[i + 1];
            }
            vertIdx1 = mesh->boundIndexToMeshId[vertIdx1];
            vertIdx2 = mesh->boundIndexToMeshId[vertIdx2];
            for (size_t j = 0; j < mesh->edges.size(); j++)
            {
                if (mesh->edges[j]->v1i == vertIdx1 && mesh->edges[j]->v2i == vertIdx2
                    ||
                    mesh->edges[j]->v1i == vertIdx2 && mesh->edges[j]->v2i == vertIdx1)
                {
                    mesh->edges[j]->isItInLongestBoundary = true;
                }
            }
        }
    }

}
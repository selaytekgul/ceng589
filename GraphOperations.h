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
    inline void generateDiskParameterization(const Mesh* mesh, const ParameterizationMethod method);
    inline std::set<Edge*> findAnyBoundaries(const Mesh* mesh);
    inline std::map<int, std::map<int, bool>> adjacencyMatrixFromEdges(const std::set<Edge*>& boundaryEdges);


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

    void generateDiskParameterization(const Mesh* mesh, const ParameterizationMethod method)
    {
        std::set<Edge*> boundaryEdges = findAnyBoundaries(mesh);
        std::map<int, std::map<int, bool>> graph = adjacencyMatrixFromEdges(boundaryEdges);
        //std::vector<std::vector<int>> longestBoundary = longestCycleConnectedComponent(graph);
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
            }
        }
        return boundaryEdges;
    }

    std::map<int, std::map<int, bool>> adjacencyMatrixFromEdges(const std::set<Edge*>& boundaryEdges) {
        // Determine the number of nodes
        //int numNodes = boundaryEdges.size();
       
        // Initialize adjacency matrix with zeros
        //std::vector<std::vector<int>> graph(numNodes, std::vector<int>(numNodes, 0));
        std::map<int, std::map<int, bool>> m_graph;

        // Fill adjacency matrix based on boundary edges
        int index = 0;
        for (const auto& edge : boundaryEdges)
        {
            m_graph[edge->v1i][edge->v2i] = true;
            m_graph[edge->v2i][edge->v1i] = true;
            //graph[edge->v1i][edge->node2] = 1;
            //graph[edge->node2][edge->node1] = 1;
            m_graph[0][1] = true;
            m_graph[0][4] = true;
        }

        //m[24][45] = "123";


        return m_graph;
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

}
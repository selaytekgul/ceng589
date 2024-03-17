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
    inline std::vector<Edge*> findAnyBoundaries(const Mesh* mesh);
    inline std::vector<Edge*> findOrderedEntireBoundaryList(std::vector<Edge*> boundaryEdges);
    inline std::vector<std::vector<Edge*>> findOrderedEntireBoundaryListImproved(std::vector<Edge*> boundaryEdges);
    inline std::vector<Edge*> returnLongestBoundary(std::vector<std::vector<Edge*>> boundaryEdges);


    inline std::vector<std::vector<int>> adjacencyMatrixFromEdges(Mesh* mesh, const std::vector<Edge*>& boundaryEdges);
    inline void printVectorOfVectors(const std::vector<std::vector<int>>& vec);
    inline void longestCycleConnectedComponentSetToEdgeInfo(Mesh* mesh, const std::vector<int>& longestBoundary);

    inline std::vector<int> longestCycleConnectedComponent(const std::vector<std::vector<int>>& graph);




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
        std::vector<Edge*> boundaryEdges = findAnyBoundaries(mesh);
        //std::vector<Edge*> boundaryList = findOrderedEntireBoundaryList(boundaryEdges);
        std::vector<std::vector<Edge*>> boundaryListImproved = findOrderedEntireBoundaryListImproved(boundaryEdges);

        //std::vector<Edge*> longestBoundaryList = returnLongestBoundary();

        std::vector<std::vector<int>> graph = adjacencyMatrixFromEdges(mesh, boundaryEdges);
        std::vector<int> longestBoundary = longestCycleConnectedComponent(graph);
        longestCycleConnectedComponentSetToEdgeInfo(mesh, longestBoundary);
        int x = 0;
    }

    std::vector<Edge*> findAnyBoundaries(const Mesh* mesh)
    {
        std::vector<Edge*> boundaryEdges;
        for (size_t i = 0; i < mesh->edges.size(); i++)
        {
            if (mesh->edges[i]->existedTriangeNumber == 1)
            {
                boundaryEdges.push_back(mesh->edges[i]);
                mesh->edges[i]->isItBoundary = true;
            }
        }
        return boundaryEdges;
    }

    void printVectorOfVectors(const std::vector<std::vector<int>>& vec)
    {
        for (const auto& innerVec : vec) {
            for (int i : innerVec) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
    }

    std::vector<std::vector<int>> adjacencyMatrixFromEdges(Mesh* mesh, const std::vector<Edge*>& boundaryEdges) {
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
        //printVectorOfVectors(adjMatrix);
        return adjMatrix;
    }


    std::vector<Edge*> findOrderedEntireBoundaryList(std::vector<Edge*> boundaryEdges)
    {
        std::vector<int> visitedBoundaryEdgeIndexes = {};
        std::vector<Edge*> listedEdges = {};
        if (boundaryEdges.size() > 0)
        {
            auto previousEdge = boundaryEdges[0];
            auto firstEdge = previousEdge;
            visitedBoundaryEdgeIndexes.push_back(previousEdge->edge_idx);
            listedEdges.push_back(previousEdge);
            while (listedEdges.size() < boundaryEdges.size())
            {
                for (size_t i = 1; i < boundaryEdges.size(); i++)
                {
                    Edge* edge = boundaryEdges[i];
                    int edgeIdx = edge->edge_idx;
                    auto iter = std::find(visitedBoundaryEdgeIndexes.begin(), visitedBoundaryEdgeIndexes.end(), edgeIdx);
                    if (iter == visitedBoundaryEdgeIndexes.end()
                        && (edge->v1i == previousEdge->v1i
                            || edge->v2i == previousEdge->v2i
                            || edge->v1i == previousEdge->v2i 
                            || edge->v2i == previousEdge->v1i))
                    {
                        listedEdges.push_back(edge);
                        visitedBoundaryEdgeIndexes.push_back(edge->edge_idx);
                        previousEdge = edge;
                    }
                }
            }
        }
        return listedEdges;
    }

    std::vector<std::vector<Edge*>> findOrderedEntireBoundaryListImproved(std::vector<Edge*> boundaryEdges)
    {
        std::vector<std::vector<Edge*>> listedBounds;
        std::vector<int> visitedBoundaryEdgeIndexes = {};
        std::vector<Edge*> listedEdges = {};
        if (boundaryEdges.size() == 0)
            return {};

        auto previousEdge = boundaryEdges[0];
        auto firstEdge = previousEdge;
        visitedBoundaryEdgeIndexes.push_back(previousEdge->edge_idx);
        listedEdges.push_back(previousEdge);

        bool entireEdgesTried = false;
        bool cycleClosed = false;
        int currentCycleSize = 1;
        int totalCyclesSize = 1;
        while (totalCyclesSize < boundaryEdges.size())
        {
            for (size_t i = 0; i < boundaryEdges.size(); i++)
            {
                Edge* edge = boundaryEdges[i];
                const int edgeIdx = edge->edge_idx;
                const auto iter = std::find(visitedBoundaryEdgeIndexes.begin(), visitedBoundaryEdgeIndexes.end(), edgeIdx);
                if (cycleClosed)
                {
                    if (iter == visitedBoundaryEdgeIndexes.end())
                    {
                        listedEdges.push_back(edge);
                        visitedBoundaryEdgeIndexes.push_back(edgeIdx);
                        previousEdge = edge;
                        firstEdge = edge;
                        cycleClosed = false;
                        currentCycleSize++;
                        totalCyclesSize++;
                    }
                }
                else
                {
                    if (iter == visitedBoundaryEdgeIndexes.end()
                        && (edge->v1i == previousEdge->v1i
                            || edge->v2i == previousEdge->v2i
                            || edge->v1i == previousEdge->v2i
                            || edge->v2i == previousEdge->v1i))
                    {
                        listedEdges.push_back(edge);
                        visitedBoundaryEdgeIndexes.push_back(edgeIdx);
                        previousEdge = edge;
                        currentCycleSize++;
                        totalCyclesSize++;
                        if (currentCycleSize > 2)
                        {
                            if ((edge->v1i == firstEdge->v1i
                                || edge->v2i == firstEdge->v2i
                                || edge->v1i == firstEdge->v2i
                                || edge->v2i == firstEdge->v1i))
                            {
                                cycleClosed = true;
                                currentCycleSize = 0;
                                break;
                            }
                        }
                    }
                }
            }
            if (cycleClosed)
            {
                listedBounds.push_back(listedEdges);
            }
        }

        return listedBounds;
    }

    std::vector<Edge*> returnLongestBoundary(std::vector<std::vector<Edge*>> boundaryEdges)
    {
        std::vector<Edge*> listedEdges = {};
        int maxLengtListIndex = 0;
        int maxLength = 0;
        int boundaryOptionIndex = 0;
        for (auto boundaryOption : boundaryEdges)
        {
            int length = 0;
            for (auto edge : boundaryOption)
            {
                length++;
            }
            if (length > maxLength)
            {
                maxLength = length;
                maxLengtListIndex = boundaryOptionIndex;
            }
            boundaryOptionIndex++;
        }
        listedEdges = boundaryEdges[maxLengtListIndex];
        return listedEdges;
    }

    void printVec(std::vector<int> vec)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            std::cout << vec.at(i) << " ";
        }
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
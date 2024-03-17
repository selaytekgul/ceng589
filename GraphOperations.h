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
    inline std::vector<Edge*> returnLongestBoundary(Mesh* mesh, std::vector<std::vector<Edge*>> boundaryEdges);

    inline void printVectorOfVectors(const std::vector<std::vector<int>>& vec);


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
        std::vector<Edge*> longestBoundaryList = returnLongestBoundary(mesh, boundaryListImproved);
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
                listedEdges = {};
            }
        }

        return listedBounds;
    }

    std::vector<Edge*> returnLongestBoundary(Mesh* mesh, std::vector<std::vector<Edge*>> boundaryEdges)
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

        std::vector<int> listedEdgeIdx = {};
        for (auto edge : listedEdges)
        {
            const int edgeIdx = edge->edge_idx;
            listedEdgeIdx.push_back(edgeIdx);
        }

        for (auto edge : mesh->edges)
        {
            const int edgeIdx = edge->edge_idx;
            const auto iter = std::find(listedEdgeIdx.begin(), listedEdgeIdx.end(), edgeIdx);
            if (iter != listedEdgeIdx.end())
            {
                edge->isItInLongestBoundary = true;
            }
        }

        return listedEdges;
    }

    void printVec(std::vector<int> vec)
    {
        for (size_t i = 0; i < vec.size(); i++)
        {
            std::cout << vec.at(i) << " ";
        }
    }
}
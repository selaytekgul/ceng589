#pragma once
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
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
#include <string>

#include <fstream>
#include <iomanip>

#include <Eigen/Dense>
#include<Eigen/SparseLU> 
namespace GraphOperations
{
    enum ParameterizationMethod
    {
        UNIFORM = 0,
        HARMONIC = 1,
        MEAN = 2
    };

    inline std::pair<float*, float> findClosestVertex(const float* source, std::vector<float*> targetList);
    inline void generateDiskParameterization(Mesh* mesh, const ParameterizationMethod method, std::string fileName);
    inline std::vector<Edge*> findAnyBoundaries(const Mesh* mesh);
    inline std::vector<Edge*> findOrderedEntireBoundaryList(std::vector<Edge*> boundaryEdges);
    inline std::vector<std::vector<Edge*>> findOrderedEntireBoundaryListImproved(std::vector<Edge*> boundaryEdges);
    inline std::vector<Edge*> returnLongestBoundary(Mesh* mesh, std::vector<std::vector<Edge*>> boundaryEdges);

    inline void parameterize(Mesh* mesh, const ParameterizationMethod method, std::vector<Edge*> boundaryEdges, std::string fileName);
    inline std::vector<std::vector<float>> createb(const Mesh* mesh, int coordinate, std::vector<Edge*> boundaryEdges);
    inline std::vector<std::vector<float>> createWUniform(const Mesh* mesh, const float weight);
    inline std::vector<std::vector<float>> createWHarmonic(Mesh* mesh);
    inline std::vector<std::vector<float>> createWMean(Mesh* mesh);

    inline Eigen::VectorXf calculateXDense(std::vector<std::vector<float>> W, std::vector<std::vector<float>> b);
    inline Eigen::VectorXf calculateXSparse(std::vector<std::vector<float>> W, std::vector<std::vector<float>> b);

    inline void printVec(std::vector<int> vec);
    inline void printMatrixToAFile(std::vector<std::vector<float>> W);
    inline void printVectorOfVectors(const std::vector<std::vector<float>>& vec);

    inline Eigen::MatrixXf createEigenMatrix(const std::vector<std::vector<float>>& data);

    inline float findAngle(Mesh* mesh, const int edgeIdx, const int vertIdx);
    inline float findAngleForTan(Mesh* mesh, const int triIdx, const int vertIdxOfi);

    inline float cot(float angle);
    inline float returnWeight(Mesh* mesh, const int vertIdxi, const int vertIdxj, const ParameterizationMethod method);

    inline void writeOFF(const Eigen::VectorXf& verticesx, const Eigen::VectorXf& verticesy, const std::string& filename);
    inline void manipulateFirstNLines(const Eigen::VectorXf& verticesx, const Eigen::VectorXf& verticesy, const std::string& readfilename, const std::string& writefilename);

    inline float calculateCotForEdge(Mesh* mesh, const int edgeIdx);
    inline float calculateTanHalfForEdge(Mesh* mesh, const int edgeIdx, int vertIdxi);


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

    void generateDiskParameterization(Mesh* mesh, const ParameterizationMethod method, std::string fileName)
    {
        std::vector<Edge*> boundaryEdges = findAnyBoundaries(mesh);
        //std::vector<Edge*> boundaryList = findOrderedEntireBoundaryList(boundaryEdges);
        std::vector<std::vector<Edge*>> boundaryListImproved = findOrderedEntireBoundaryListImproved(boundaryEdges);
        std::vector<Edge*> longestBoundaryList = returnLongestBoundary(mesh, boundaryListImproved);
        parameterize(mesh, method, longestBoundaryList, fileName);

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

    void printVectorOfVectors(const std::vector<std::vector<float>>& vec)
    {
        for (const auto& innerVec : vec) {
            for (int i : innerVec) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
            std::cout << std::endl;
            std::cout << "----------------------------------" << std::endl;
            std::cout << std::endl;
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
        if (boundaryEdges.size() == 0)
            return {};
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
                mesh->verts[edge->v1i]->isItInLongestBoundary = true;
                mesh->verts[edge->v2i]->isItInLongestBoundary = true;
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

    void parameterize(Mesh* mesh, const ParameterizationMethod method, std::vector<Edge*> boundaryEdges, std::string fileName)
    {
        std::vector<std::vector<float>> bx = createb(mesh, 0, boundaryEdges);
        std::vector<std::vector<float>> by = createb(mesh, 1, boundaryEdges);
        int length = mesh->verts.size();
        std::vector<std::vector<float>> W = {};
        switch (method)
        {
        case ParameterizationMethod::UNIFORM:
            {
                float weight = 1.0;
                W = createWUniform(mesh, weight);
                //printVectorOfVectors(W);
                //printMatrixToAFile(W);
                int a = 1;
            }
            break;
        case ParameterizationMethod::HARMONIC:
            {
                W = createWHarmonic(mesh);
                //printVectorOfVectors(W);
                //printMatrixToAFile(W);
                int a = 1;
            }
            break;
        case ParameterizationMethod::MEAN:
            {
                W = createWMean(mesh);
                //printVectorOfVectors(W);
                //printMatrixToAFile(W);
                int a = 1;
            }
            break;
        default:
            break;
        }

        const Eigen::VectorXf& verticesx = calculateXDense(W, bx);
        const Eigen::VectorXf& verticesy = calculateXDense(W, by);
        std::string readfn = fileName + ".off";
        std::string writefn = fileName + std::to_string(method) + ".off";
        manipulateFirstNLines(verticesx, verticesy, readfn, writefn);
    }

    Eigen::VectorXf calculateXSparse(std::vector<std::vector<float>> W, std::vector<std::vector<float>> b)
    {
        ////Eigen::MatrixXf matrixW = createEigenMatrix(W);
        //Eigen::MatrixXf matrixW = createEigenMatrix(W).transpose();

        //Eigen::VectorXf eigenVector(W.size());
        //Eigen::VectorXf x(W.size());

        //// Fill the Eigen vector with data from std::vector<std::vector<float>>
        //for (int i = 0; i < W.size(); ++i) {
        //    eigenVector(i) = b[i][0];
        //}

        ////std::cout << std::endl << eigenVector << std::endl;
        //Eigen::VectorXf x = matrixW.colPivHouseholderQr().solve(eigenVector);


        ////VectorXd x(n), b(n);
        //SparseMatrix<double> A;
        //SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
        //// fill A and b;
        //// Compute the ordering permutation vector from the structural pattern of A
        //solver.analyzePattern(A);
        //// Compute the numerical factorization 
        //solver.factorize(A);
        ////Use the factors to solve the linear system 
        //x = solver.solve(b);
        return {};
    }
    
    std::vector<std::vector<float>> createWUniform(const Mesh* mesh, const float weight)
    {
        int length = mesh->verts.size();
        std::vector<std::vector<float>> W(length, std::vector<float>(length, 0.0));
        for (size_t i = 0; i < length; i++)
        {
            Vertex* vert = mesh->verts[i];
            if (vert->isItInLongestBoundary)
            {
                W[i][i] = 1.0;
            }
            else
            {
                int numberOfNeighbors = vert->vertList.size();
                for (size_t j = 0; j < numberOfNeighbors; j++)
                {
                    int vertIdx = vert->vertList[j];
                    W[i][vertIdx] = weight;
                }
                W[i][i] = -numberOfNeighbors;
            }
        }
        return W;
    }
    
    std::vector<std::vector<float>> createWHarmonic(Mesh* mesh)
    {
        int length = mesh->verts.size();
        std::vector<std::vector<float>> W(length, std::vector<float>(length, 0.0));
        for (size_t i = 0; i < length; i++)
        {
            Vertex* vert = mesh->verts[i];
            if (vert->isItInLongestBoundary)
            {
                W[i][i] = 1.0;
            }
            else
            {
                int numberOfNeighbors = vert->vertList.size();
                float weightsSum = 0.0;
                for (size_t j = 0; j < numberOfNeighbors; j++)
                {
                    if (i == j)
                        continue;
                    int vertIdx = vert->vertList[j];
                    float weight = returnWeight(mesh, i, vertIdx, ParameterizationMethod::HARMONIC);
                    W[i][vertIdx] = weight;
                    weightsSum += weight;
                }
                W[i][i] = -weightsSum;
            }
        }
        return W;
    }
        
    std::vector<std::vector<float>> createWMean(Mesh* mesh)
    {
        int length = mesh->verts.size();
        std::vector<std::vector<float>> W(length, std::vector<float>(length, 0.0));
        for (size_t i = 0; i < length; i++)
        {
            Vertex* vert = mesh->verts[i];
            if (vert->isItInLongestBoundary)
            {
                W[i][i] = 1.0;
            }
            else
            {
                int numberOfNeighbors = vert->vertList.size();
                float weightsSum = 0.0;
                for (size_t j = 0; j < numberOfNeighbors; j++)
                {
                    if (i == j)
                        continue;
                    int vertIdx = vert->vertList[j];
                    float weight = returnWeight(mesh, i, vertIdx, ParameterizationMethod::MEAN);
                    W[i][vertIdx] = weight;
                    weightsSum += weight;
                }
                W[i][i] = -weightsSum;
            }
        }
        return W;
    }
    
    std::vector<std::vector<float>> createb(const Mesh* mesh, int coordinate, std::vector<Edge*> boundaryEdges)
    {
        auto perimeter = 0.0;
        std::map<int, float> boundaryEdgeIndexToAngle = {};

        for (auto edge : boundaryEdges)
        {
            Vertex* v1 = mesh->verts[edge->v1i];
            Vertex* v2 = mesh->verts[edge->v2i];
            float v1_coords[3] = { v1->coords[0], v1->coords[1], v1->coords[2] };
            float v2_coords[3] = { v2->coords[0], v2->coords[1], v2->coords[2] };
            float distance = VectorMath::distanceBetweenVectors(v1_coords, v2_coords);
            boundaryEdgeIndexToAngle[edge->edge_idx] = 0;
            perimeter += distance;
        }

        for (auto edge : boundaryEdges)
        {
            Vertex* v1 = mesh->verts[edge->v1i];
            Vertex* v2 = mesh->verts[edge->v2i];
            float v1_coords[3] = { v1->coords[0], v1->coords[1], v1->coords[2] };
            float v2_coords[3] = { v2->coords[0], v2->coords[1], v2->coords[2] };
            float distance = VectorMath::distanceBetweenVectors(v1_coords, v2_coords);
            float angle = 2 * M_PI * distance / perimeter;
            boundaryEdgeIndexToAngle[edge->edge_idx] = angle;
        }

        const int length = mesh->verts.size();
        std::vector<std::vector<float>> b(length, std::vector<float>(1, 0.0));
        //printVectorOfVectors(b);

        const int boundaryEdgesLength = boundaryEdges.size();
        if (boundaryEdgesLength == 0)
            return b;
        Edge* firstEdge = boundaryEdges[0];
        int firstEdgeIndex = firstEdge->edge_idx;
        Vertex* firstV1 = mesh->verts[firstEdge->v1i];
        Vertex* firstV2 = mesh->verts[firstEdge->v2i];
        int prevV1idx = firstV1->idx;
        int prevV2idx = firstV2->idx;
        float angle = boundaryEdgeIndexToAngle[firstEdgeIndex];
        
        const float r = 1.0;
        b[firstV1->idx][0] = coordinate == 0 ? r * cos(0) : r * sin(0);
        b[firstV2->idx][0] = coordinate == 0 ? r * cos(angle) : r * sin(angle);
        for (size_t i = 1; i < boundaryEdgesLength - 1; i++)
        {
            Edge* edge = boundaryEdges[i];
            int edgeIndex = edge->edge_idx;
            angle += boundaryEdgeIndexToAngle[edgeIndex];

            Vertex* V1 = mesh->verts[edge->v1i];
            Vertex* V2 = mesh->verts[edge->v2i];
            int currentV1idx = V1->idx;
            int currentV2idx = V2->idx;
            if (currentV1idx == prevV1idx || currentV1idx == prevV2idx)
            {
                b[V2->idx][0] = coordinate == 0 ? r * cos(angle) : r * sin(angle);
            }
            else if (currentV2idx == prevV1idx || currentV2idx == prevV2idx)
            {
                b[V1->idx][0] = coordinate == 0 ? r * cos(angle) : r * sin(angle);
            }
            prevV1idx = V1->idx;
            prevV2idx = V2->idx;
        }
        //printVectorOfVectors(b);

        return b;
    }
    
    Eigen::VectorXf calculateXDense(std::vector<std::vector<float>> W, std::vector<std::vector<float>> b)
    {
        //Eigen::MatrixXf matrixW = createEigenMatrix(W);
        Eigen::MatrixXf matrixW = createEigenMatrix(W).transpose();

        Eigen::VectorXf eigenVector(W.size());

        // Fill the Eigen vector with data from std::vector<std::vector<float>>
        for (int i = 0; i < W.size(); ++i) {
            eigenVector(i) = b[i][0];
        }

        //std::cout << std::endl << eigenVector << std::endl;
        Eigen::VectorXf x = matrixW.colPivHouseholderQr().solve(eigenVector);

        //std::cout << x << std::endl;
        //printVectorOfVectors(W);
        //std::ofstream fout;
        //fout.open("matrixW.txt");
        //std::cout << matrixW;
        //fout.close();

        return x;
    }

    void writeOFF(const Eigen::VectorXf& verticesx, const Eigen::VectorXf& verticesy, const std::string& filename) {
        std::ofstream ofs(filename);
        if (!ofs.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        // Write the header for OFF file
        ofs << "OFF" << std::endl;
        ofs << verticesx.size()<< " 0 0" << std::endl;

        // Write the vertices
        for (int i = 0; i < verticesx.size(); i++) {
            ofs << verticesx(i) << " " << verticesy(i) << " 0" << std::endl;
        }

        ofs.close();
        std::cout << "OFF file saved: " << filename << std::endl;
    }

    void manipulateFirstNLines(const Eigen::VectorXf& verticesx, const Eigen::VectorXf& verticesy, const std::string& readfilename, const std::string& writefilename) {
        std::ifstream ifs(readfilename);
        if (!ifs.is_open()) {
            std::cerr << "Error opening file: " << readfilename << std::endl;
            return;
        }

        // Read the entire file into a vector of strings
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(ifs, line)) {
            lines.push_back(line);
        }

        ifs.close();
        //lines[0] = "OFF\n";


        std::ofstream ofscopy(writefilename);
        if (!ofscopy.is_open()) {
            std::cerr << "Error opening file for writing: " << writefilename << std::endl;
            return;
        }

        for (const std::string& line : lines) {
            ofscopy << line << std::endl;
        }

        ofscopy.close();



        // Manipulate the first n lines
        for (int i = 2; i < verticesx.size() + 2 && i < lines.size(); ++i) {
            // Manipulate line i here
            // For example, you can print the line:
            //std::cout << "Original line " << i << ": " << lines[i] << std::endl;
            // Or you can modify the line, e.g., lines[i] = "new content";
            lines[i] = std::to_string(verticesx(i-2)) + " " + std::to_string(verticesy(i-2)) + " 0";
        }

        // Write the modified lines back to the file
        std::ofstream ofs(writefilename);
        if (!ofs.is_open()) {
            std::cerr << "Error opening file for writing: " << writefilename << std::endl;
            return;
        }

        for (const std::string& line : lines) {
            ofs << line << std::endl;
        }

        ofs.close();

        std::cout << "Manipulation complete for the first " << std::to_string(verticesx.size()) << " lines of " << writefilename << std::endl;
    }

    void printMatrixToAFile(std::vector<std::vector<float>> W)
    {
        std::ofstream fout;
        fout.open("matrixW.txt");
        for (int i = 0; i < W.size(); ++i)
        {
            for (size_t j = 0; j < W.size(); j++)
            {
                fout << std::fixed << std::setprecision(4) << W[i][j] << " ";
            }
            fout << std::endl;
        }
        fout.close();
    }

    Eigen::MatrixXf createEigenMatrix(const std::vector<std::vector<float>>& data) {
        // Get the dimensions of the data
        int rows = static_cast<float>(data.size());
        int cols = rows > 0 ? static_cast<float>(data[0].size()) : 0;

        // Create a dynamically allocated array to hold the data
        float* dataArray = new float[rows * cols];

        // Copy data from vector<vector<int>> to the dynamically allocated array
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                dataArray[i * cols + j] = data[i][j];
            }
        }

        // Create Eigen MatrixXx using the copied data
        Eigen::Map<Eigen::MatrixXf> eigenMatrix(dataArray, rows, cols);

        // Optionally, you can copy the data if you want to keep it after deleting the original array
        Eigen::MatrixXf copiedEigenMatrix = eigenMatrix;
        //std::cout << copiedEigenMatrix;
        // Free the dynamically allocated memory
        delete[] dataArray;

        //std::cout << std::endl << std::endl << copiedEigenMatrix << std::endl;

        // Return the Eigen MatrixXx

        return copiedEigenMatrix;
    }
    float cot(float angle) {
        float retVal = sinf(angle) == 0 ? 0 : cosf(angle) / sinf(angle);
        return retVal;
    }

    float returnWeight(Mesh* mesh, const int vertIdxi, const int vertIdxj, const ParameterizationMethod method)
    {
        float weight = 0.0;
        switch (method)
        {
        case ParameterizationMethod::UNIFORM:
        {
            weight = 1.0;
        }
        break;
        case ParameterizationMethod::HARMONIC:
        {
            if (vertIdxi == vertIdxj)
                return 0.0f;
            const Vertex* verti = mesh->verts[vertIdxi];
            const Vertex* vertj = mesh->verts[vertIdxj];
            const int neighEdgeNum = verti->edgeList.size();
            for (size_t eachNeighEdge = 0; eachNeighEdge < neighEdgeNum; eachNeighEdge++)
            {
                const int edgeIdx = verti->edgeList[eachNeighEdge];
                Edge* edge = mesh->edges[edgeIdx];
                const int endPoint1Idx = edge->v1i;
                const int endPoint2Idx = edge->v2i;
                if (endPoint1Idx == vertIdxi && endPoint2Idx == vertIdxj
                    || endPoint2Idx == vertIdxi && endPoint1Idx == vertIdxj)
                {
                    weight = calculateCotForEdge(mesh, edgeIdx);
                    return weight;
                }
            }
            return weight;
        }
            break;
        case ParameterizationMethod::MEAN:
        {
            if (vertIdxi == vertIdxj)
                return 0.0f;
            const Vertex* verti = mesh->verts[vertIdxi];
            const Vertex* vertj = mesh->verts[vertIdxj];
            const int neighEdgeNum = verti->edgeList.size();
            for (size_t eachNeighEdge = 0; eachNeighEdge < neighEdgeNum; eachNeighEdge++)
            {
                const int edgeIdx = verti->edgeList[eachNeighEdge];
                Edge* edge = mesh->edges[edgeIdx];
                const int endPoint1Idx = edge->v1i;
                const int endPoint2Idx = edge->v2i;
                if (endPoint1Idx == vertIdxi && endPoint2Idx == vertIdxj
                    || endPoint2Idx == vertIdxi && endPoint1Idx == vertIdxj)
                {
                    weight = calculateTanHalfForEdge(mesh, edgeIdx, vertIdxi);
                    return weight;
                }
            }
            return weight;
        }
            break;
        default:
            weight = 0.0f;
            break;
        }
        return weight;
    }

    float findAngle(Mesh* mesh, const int edgeIdx, const int vertIdx)
    {
        const float* edgeVert1coords = mesh->verts[mesh->edges[edgeIdx]->v1i]->coords;
        const float* edgeVert2coords = mesh->verts[mesh->edges[edgeIdx]->v2i]->coords;
        
        const float* otherVertcoords = mesh->verts[vertIdx]->coords;

        float vector1[3] = {};
        float vector2[3] = {};
        VectorMath::vector(vector1, edgeVert1coords, otherVertcoords);
        VectorMath::vector(vector2, edgeVert2coords, otherVertcoords);

        const float angle = VectorMath::calculateAngleBetweenVectors(vector1, vector2);
        return angle;
    }


    float findAngleForTan(Mesh* mesh, const int triIdx, const int vertIdxOfi)
    {
        const float* myVertcoords = mesh->verts[vertIdxOfi]->coords;

        triVertsIds triVertIdxs = TriangleMeshMath::getVertexIdsOfTriangle(mesh, triIdx);
        int selectedVertNumInTri = 0;
        for (size_t i = 0; i < triVertIdxs.size(); i++)
        {
            if (vertIdxOfi == triVertIdxs[i])
                selectedVertNumInTri = i;
        }

        triVertsCoords tri3VertsCoords = TriangleMeshMath::getCoordsOfTriangle(mesh, triVertIdxs);
        triOtherVertsCoords tri2VertsCoords = TriangleMeshMath::getOtherCoordsOfTriangle(tri3VertsCoords, selectedVertNumInTri);

        const float* otherVert1coords = tri2VertsCoords[0].data();
        const float* otherVert2coords = tri2VertsCoords[1].data();

        float vector1[3] = {};
        float vector2[3] = {};
        VectorMath::vector(vector1, otherVert1coords, myVertcoords);
        VectorMath::vector(vector2, otherVert2coords, myVertcoords);

        const float angle = VectorMath::calculateAngleBetweenVectors(vector1, vector2);
        return angle;
    }


    float calculateCotForEdge(Mesh* mesh, const int edgeIdx)
    {
        float totalCot = 0.0;
        Edge* edge = mesh->edges[edgeIdx];
        const int neighborTriangleNum = edge->existedTriangeNumber;
        //mesh->computeLength(edgeIdx);
        //const float edgeLength = edge->length;
        for (size_t i = 0; i < neighborTriangleNum; i++)
        {
            int triIdx = edge->triList[i];
            int vertIdx = TriangleMeshMath::getOtherVertexIdOfTriangle(mesh, triIdx, edgeIdx);
            float angle = findAngle(mesh, edgeIdx, vertIdx);
            totalCot += cot(angle);
        }
        return totalCot/static_cast<float>(neighborTriangleNum);
    }

    float calculateTanHalfForEdge(Mesh* mesh, const int edgeIdx, int vertIdxi)
    {
        float totalTanHalf = 0.0;
        Edge* edge = mesh->edges[edgeIdx];
        const int neighborTriangleNum = edge->existedTriangeNumber;
        
        for (size_t i = 0; i < neighborTriangleNum; i++)
        {
            int triIdx = edge->triList[i];
            int vertIdx = TriangleMeshMath::getOtherVertexIdOfTriangle(mesh, triIdx, edgeIdx);

            float angle = findAngleForTan(mesh, triIdx, vertIdxi);
            //totalTanHalf += cot(angle); //Creates improved (more closed mouth)

            totalTanHalf += tan(angle / 2.0f);
        }

        mesh->computeLength(edgeIdx);
        const float edgeLength = edge->length;
        return totalTanHalf / (static_cast<float>(neighborTriangleNum) * edgeLength);
        //return totalTanHalf / static_cast<float>(neighborTriangleNum);

    }

}
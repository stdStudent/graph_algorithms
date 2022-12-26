#ifndef LAB_3_GRAPH_H
#define LAB_3_GRAPH_H

#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <fstream>
#include <limits>
#include <list>
#include <algorithm>

using namespace std;

typedef vector<int> vint;
typedef vector<vint> vvint;

const int INF = std::numeric_limits<int>::max();

class Graph {
    struct Edge {
        int weight;
        int first;
        int second;

        Edge(int w, int f, int s) {
            weight = w;
            first = f;
            second = s;
        }
    };

    vector<vector<int>> AdjMatrix;
    vector<vector<int>> AdjMatrixINF;
    vector<vector<int>> AdjList;
    vector<Edge> Edges;
    int n_verts = 0;
    int n_edges = 0;

    void convertToAdjList() {
        vector<vector<int>> tmp(/*AdjMatrix.size()*/ n_verts);

        for (int i = 0; i < /*AdjMatrix.size()*/ n_verts; i++)
            for (int j = 0; j < /*AdjMatrix[i].size()*/ n_verts; j++)
                if (AdjMatrix[i][j] == 1)
                    tmp[i].push_back(j);

        AdjList = tmp;
    }

    void convertToAdjMatrixINF() {
        vector<vector<int>> tmp(n_verts, vector<int>(n_verts,0));

        for (int i = 0; i < n_verts; ++i) {
            for (int j = 0; j < n_verts; ++j) {
                if (i != j && AdjMatrix[i][j] == 0)
                    tmp[i][j] = INF;
                else
                    tmp[i][j] = AdjMatrix[i][j];
            }
        }

        AdjMatrixINF = tmp;
    }

    void convertToListEdges() {
        for (int i = 0; i < n_verts; ++i) {
            for (int j = 0; j < n_verts; ++j) {
                if (AdjMatrix[i][j] != 0)
                    Edges.emplace_back(AdjMatrix[i][j], i, j);
            }
        }
    }

    void writeMatrixToFile(vector<vector<int>> matrix, string file_path) {
        ofstream f;
        f.open(file_path);
        if (f.is_open()) {
            for (int i = 0; i < n_verts; i++) {
                for (int j = 0; j < n_verts; j++) {
                    f << matrix[i][j] << "\t";
                }
                f << "\n";
            }
            f.close();
        } else
            cout << "error";
    }

    void writePointsToFile(const vector<pair<int,int>>& points, const string& file)
    {
        vector <pair<int,int>>:: const_iterator i;

        if (!file.empty())
        {
            ofstream out(file.c_str());
            if(out.fail()) {
                out.close();
            }

            for(i=points.begin(); i != points.end();++i ) {
                out << i->first << " " << i->second << "\n";
            }

            out.close();
        }
    }

    void writeVectorToFile(vector<int> matrix, string file_path) {
        ofstream f;
        f.open(file_path);
        if (f.is_open()) {
            for (int i = 0; i < n_verts; i++)
                f << matrix[i] << ' ';

            f.close();
        } else
            cout << "error";
    }

public:
    [[nodiscard]] const vector<vector<int>> &getAdjMatrix() const {
        return AdjMatrix;
    }

    [[nodiscard]] const vector<vector<int>> &getAdjList() const {
        return AdjList;
    }

    explicit Graph(const string& filename) {
        fstream test_file;
        test_file.open(filename, ios::in);

        int n = 0, e = 0;
        test_file >> n;
        n_verts = n;
        vector<vector<int>> tmp(n,vector<int>(n,0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                test_file >> tmp[i][j];
                if (tmp[i][j] != 0)
                    ++e;
            }
        }

        test_file.close();
        AdjMatrix = tmp;
        convertToAdjList();
        convertToAdjMatrixINF();

        n_edges = e;
        convertToListEdges();
    }

    void printAdjMatrix() {
        for (int i = 0; i < /*AdjMatrix.size()*/ n_verts; i++)
        {
            for (int j = 0; j < /*AdjMatrix[i].size()*/ n_verts; j++)
                cout << AdjMatrix[i][j] << ' ';
            cout << endl;
        }
    }

    void printAdjList() {
        for (int i = 0; i < AdjList.size(); i++) {
            if (AdjList[i].empty())
                continue;

            cout << i << ':';

            for (int j = 0; j < AdjList[i].size(); j++)
                if (j == AdjList[i].size() - 1) {
                    cout << " -> " << AdjList[i][j] << endl;
                    break;
                } else
                    cout << " -> " << AdjList[i][j];
        }
    }

    void floydWarshall() {
        vector<vector<int>> matrix_res = AdjMatrixINF;

        auto start = chrono::steady_clock::now();
        for (int k = 0; k < n_verts; k++) {
            for (int i = 0; i < n_verts; i++) {
                for (int j = 0; j < n_verts; j++) {
                    if (matrix_res[i][k] == INF || matrix_res[k][j] == INF)
                        continue;
                    if (matrix_res[i][k] + matrix_res[k][j] < matrix_res[i][j])
                        matrix_res[i][j] = matrix_res[i][k] + matrix_res[k][j];
                }
            }
        }
        auto stop = std::chrono::steady_clock::now();
        chrono::duration<double> diff = stop - start;
        cout << "Time: " << diff.count() << " seconds." << endl;

        writeMatrixToFile(matrix_res, "floyd.txt");
    }

    int EdmondKarp(int from = -1, int to = -1) {
        if (from < 0) from = 0;
        if (to < 0) to = n_verts-1;

        vvint C = AdjMatrix; //	flow capacity
        int n = n_verts;

        vvint D(n, vint(n)); // residual network

        auto start = chrono::steady_clock::now();
        while(true) {
            vint P(n, -1); // Edges of the graph of the shortest paths
            vint q(n);

            q[0] = from;
            P[from] = 0;

            // (2)
            int h = 0, t = 1;
            while (h < t) { // Run through the graph vertexes
                int cur = q[h++];
                for (int v = 0; v < n; v++)
                    // P[v] == -1  <=>  not visited
                    if (P[v] == -1 && C[cur][v] - D[cur][v] > 0) {
                        q[t++] = v;
                        P[v] = cur;
                    }
            }

            // (2)
            if (P[to] == -1)
                break;

            // (3.a)
            int c_min = INF, cur = to;
            while (cur != from) {
                int prev = P[cur];
                c_min = min(c_min, C[prev][cur] - D[prev][cur]);
                cur = prev;
            }

            // (3.b, 3.C)
            cur = to;
            while (cur != from) {
                int prev = P[cur];
                D[prev][cur] += c_min;
                D[cur][prev] -= c_min;
                cur = prev;
            }

        }

        int phi_flow_size = 0;
        for (int i=0; i<n; i++)
            if (C[from][i])
                phi_flow_size += D[from][i];

        auto stop = std::chrono::steady_clock::now();
        chrono::duration<double> diff = stop - start;
        cout << "Time: " << diff.count() << " seconds." << endl;

        writeMatrixToFile(D, "Edmond_Karp.txt");
        return phi_flow_size;
    }

    int Kruskal() {
        sort(Edges.begin(), Edges.end(), [](Edge a, Edge b){return a.weight < b.weight;});

        vector<int> comp(n_verts);
        vector<pair<int, int>> result;
        for (int i = 0; i < n_verts; ++i)
            comp[i] = i;
        int ans = 0;


        auto start = chrono::steady_clock::now();
        for (auto& edge: Edges)
        {
            int weight = edge.weight;
            int begin  = edge.first;
            int end    = edge.second;

            if (comp[begin] != comp[end])
            {
                result.emplace_back(begin, end);
                ans += weight;
                int a = comp[begin];
                int b = comp[end];
                for (int i = 0; i < n_verts; ++i)
                    if (comp[i] == b)
                        comp[i] = a;
            }
        }
        auto stop = std::chrono::steady_clock::now();
        chrono::duration<double> diff = stop - start;
        cout << "Time: " << diff.count() << " seconds." << endl;

        if (result.size() != (n_verts - 1)) {
            cout << "The graph is not connected!\n";
            return -1;
        }

        writePointsToFile(result, "Kruskal.txt");
        return ans;
    }

    int Prim() {
        int index = 0;

        std::vector<bool> visited(n_verts, false);
        visited[0] = true;

        int ans = 0;
        vector<int> p(n_verts, INF);
        p[0] = 0;

        vector<int> result(n_verts);
        result[0] = -1;


        auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < n_verts; ++i) {
            int min_weight = INF;
            for (int j = 0; j < n_verts; ++j) {
                if (visited[j])
                    continue;
                if (p[j] < min_weight) {
                    min_weight = p[j];
                    index = j;
                }
            }

            if (min_weight < INF)
                ans += min_weight;
            visited[index] = true;

            for (int z = 0; z < n_verts; ++z) {
                if (AdjMatrixINF[index][z] == 0 || visited[z])
                    continue;
                p[z] = min(AdjMatrixINF[index][z], p[z]);
                result[z] = index;
            }
        }

        auto stop = std::chrono::steady_clock::now();
        chrono::duration<double> diff = stop - start;
        cout << "Time: " << diff.count() << " seconds." << endl;

        if (find(visited.begin(), visited.end(), false) != visited.end()) {
            cout << "The graph is not connected!\n";
            return -1;
        }

        writeVectorToFile(result, "Prim.txt");
        return ans;
    }

};


#endif //LAB_3_GRAPH_H

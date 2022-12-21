#ifndef LAB_3_GRAPH_H
#define LAB_3_GRAPH_H

#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <fstream>
#include <limits>
#include <list>

using namespace std;

typedef vector<int> vint;
typedef vector<vint> vvint;

const int INF = std::numeric_limits<int>::max();

class Graph {
    vector<vector<int>> AdjMatrix;
    vector<vector<int>> AdjMatrixINF;
    vector<vector<int>> AdjList;
    list<int> Edges;
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
                    Edges.push_back(AdjMatrix[i][j]);
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

        //n_edges = e;
        //convertToListEdges();
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

        vvint c = AdjMatrix;
        int n = n_verts;

        vvint D(n, vint(n)); // residual network

        auto start = chrono::steady_clock::now();
        while(true) {
            vint P(n, -1); // Edges of the graph of the shortest paths
            vint q(n);

            q[0] = from;
            P[from] = 0;

            // (3.c)
            int h = 0, t = 1;
            while (h < t) {
                int cur = q[h++];
                for (int v = 0; v < n; v++)
                    if (P[v] == -1 && c[cur][v] - D[cur][v] > 0) {
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
                c_min = min(c_min, c[prev][cur] - D[prev][cur]);
                cur = prev;
            }

            // (3.b)
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
            if (c[from][i])
                phi_flow_size += D[from][i];

        auto stop = std::chrono::steady_clock::now();
        chrono::duration<double> diff = stop - start;
        cout << "Time: " << diff.count() << " seconds." << endl;

        writeMatrixToFile(D, "Edmond_Karp.txt");
        return phi_flow_size;
    }

    bool bfs(int s, int t, vector<vector<int>> D, vector<int>& p) {
        vector<bool> visited(n_verts, false);
        queue <int> q;
        q.push(s);
        visited[s] = true;
        for (unsigned int i = 0; i < p.size(); ++i)
            p[i] = -1;
        p[s] = 0;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int v = 0; v < n_verts; v++)
                if (visited[v] == false && AdjMatrix[u][v] - D[u][v] > 0) {
                    q.push(v);
                    p[v] = u;
                    visited[v] = true;
                }
        }
        return visited[t] == true;
    }

    int ef(int s, int t) {
        int flow  = 0;
        vvint D(n_verts, vector<int>(n_verts, 0));
        vvint c = AdjMatrix;
        vector<int> p(n_verts);

        auto start = chrono::steady_clock::now();
        while (bfs(s, t, D, p)) {
            int cmin = INF;
            for (int v = t; v != s; v = p[v]) {
                int u = p[v];
                cmin = min(cmin, c[u][v] - D[u][v]);
            }

            for (int v = t; v != s; v = p[v]) {
                int u = p[v];
                D[u][v] += cmin;
                D[v][u] -= cmin;
            }

            flow += cmin;
        }

        auto stop = std::chrono::steady_clock::now();
        chrono::duration<double> diff = stop - start;
        cout << "Time: " << diff.count() << " seconds." << endl;

        writeMatrixToFile(D, "lab3_result.txt");

        return flow;
    }

};


#endif //LAB_3_GRAPH_H

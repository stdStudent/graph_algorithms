#include <iostream>

#include "graph.h"

using namespace std;

int main() {
    Graph G("SST/test5.txt");
    //G.printAdjMatrix();
    //G.printAdjList();
    //G.floydWarshall();

    //cout << G.EdmondKarp(0, 6);
    cout << G.Kruskal();
    cout << '\n' << G.Prim();

    return 0;
}



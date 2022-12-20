#include <iostream>

#include "graph.h"

using namespace std;

int main() {
    Graph G("6.txt");
    G.printAdjMatrix();
    G.printAdjList();
    //G.floydWarshall();

    cout << G.EdmondKarp(0, 111);
    // 0, 100 =
    //cout << G.maxFlow(0, 100);
    cout << '\n' << G.ef(0, 111);

    return 0;
}



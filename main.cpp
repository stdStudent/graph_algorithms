#include <iostream>

#include "graph.h"

using namespace std;

int main() {
    Graph G("6.txt");
    //G.printAdjMatrix();
    //G.printAdjList();
    //G.floydWarshall();

    cout << G.EdmondKarp(1, 499);
    // 0, 100 =
    //cout << G.maxFlow(0, 100);
    cout << '\n' << G.ef(1, 499);

    return 0;
}



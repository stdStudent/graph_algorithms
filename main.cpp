#include <iostream>

#include "graph.h"

using namespace std;

int main() {
    Graph G("6.txt");
    //G.printAdjMatrix();
    //G.printAdjList();
    //G.floydWarshall();

    cout << G.EdmondKarp(0,499);

    return 0;
}



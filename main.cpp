#include <iostream>

#include "graph.h"

using namespace std;

int main() {
    Graph G("BK/test4.txt");
    //G.printAdjMatrix();
    //G.printAdjList();
    //G.floydWarshall();

    //cout << G.EdmondKarp(0, 6);
    //cout << G.Kruskal();
    //cout << '\n' << G.Prim();
    auto [first, second] = G.BronKerbosch();
    cout << first << ' ' << second << '\n';

    //G.BronKerbosch();

//    Graph GG()
    return 0;
}



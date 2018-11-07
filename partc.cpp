/*
NAME: partc.cpp
COMPILATION: g++ partc.cpp -o PA2C
DESCRIPTION: reads in  an MST (distance) matrix of an undirected weighted graph
and outputs the DFS
AUTHOR: Amberlyn Schjoll
REFERENCES: This code uses https://www.geeksforgeeks.org/depth-first-search-or-dfs-for-a-graph/ for the algorithm
*/

#include<iostream>
#include<list>
#include<vector>
#include <stdlib.h>

using namespace std;

class Graph
{
    int V;    // No. of vertices

    list<int> *adj;

    void DFSUtil(int v, bool visited[]);
public:
    Graph(int V);   // Constructor

    void addEdge(int v, int w);

    void DFS(int v);
};

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w); // Add w to v’s list.
}

void Graph::DFSUtil(int v, bool visited[])
{
    visited[v] = true;
    cout << v+1 << " ";

    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            DFSUtil(*i, visited);
}

void Graph::DFS(int v)
{
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;

    DFSUtil(v, visited);
}

int main()
{
int num_verts, weight, start, edge; 
//#of vertices, edge weight, starting vertex, edge
vector <int> rowVec;
vector <int> colVec;

cin >> weight >> num_verts;

for (int i = 0; i < num_verts; i++){
    for (int j = 0; j < num_verts; j++){
        cin >> edge;
        if (edge != 0)
	{
            rowVec.push_back(i);
            colVec.push_back(j);
        }
    }
}

cin >> start;

Graph g(num_verts);
for (int i = 0; i < 2*num_verts - 2; i++)
{
    g.addEdge(rowVec[i], colVec[i]);
}

g.DFS(start - 1);
cout << endl;

return 0;
}

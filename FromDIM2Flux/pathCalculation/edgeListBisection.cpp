// Known inefficiencies:
// 1. Graph is represented by an adjacency matrix
// 2. Priority queue is implemented by an array
// -------------------------------------------------------------
//              History
// -------------------------------------------------------------
// written by: shruthi viswanath
/*
Note about additions/modifications from Metzner, Schutte and Van den Eijnden, 2009 
1. The paper talks about running the bisected edge-elimination algorithm on reduced graphs in each iteration. 
But here I don't reduce graphs. It is not clear that the additional time required to set up reduced graphs, and sort edges at each iteration (making sure that we don't count the already existing vertices on the path will compensate for the few extra BFS's we will need to do, each time, to find a bottleneck). This *might* be worth adding in the future if we really need it. 

2. The algorithm in the paper fails if there are cycles, so I had to modify the BFS search for connectedness to make sure that we don't traverse vertices that have already been encountered on the global maximum weight path (by storing already encountered vertices in a hash table)

3. The algorithm does not work if there are degenerate edges, so I make all edges explicitly non-degenerate by sorting them first and making them all uniform-spaced. e.g. 1,2,3...numEdges

4. For large graphs, we need double-precision in order for this algorithm to work correctly. 

5. I am deleting all the bottleneck edges from the global edgeList, this might make runtime a wee-bit shorter. I am also deleting those edges from the graph representation (g.A)

*/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <list>
#include <queue>
#include <vector>
#include <algorithm>
#include <ext/hash_set>
#include <sys/time.h> 

using namespace std;
using namespace __gnu_cxx;

#define FILENAME_LEN 150

#define ERR -1 // return code from functions for error message
#define FINE 0 // return code for successful completion

#define min(x,y) (x<y)?x:y

int numNodes;	// the number of nodes
int start;	// starting node for path calculation
int end;	// final node for path calculation

static long myclock()
{
    struct timeval tv; 
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + tv.tv_usec;
}


double getRuntime(long* end, long* start)
{
    return (*end - *start);
}

typedef struct {
int u;
int v;
double weight;
} edge; 

bool edgeSorter(edge  a, edge  b)
{ 
	return a.weight<b.weight;
}

vector<edge> edgeList; // global list of edges used 

hash_set<int> verticesOnPath; // global list of vertices on the path. This is stored in order to avoid cycles. 
// -------------------------------------------------------------
//              Object Definitions
// -------------------------------------------------------------

class Graph {				// define the Graph class
public:
	int n;				// includes an 'n' integer
	double **A;			// and a pointer to a pointer called 'A', adjacency matrix

	Graph(int size = 2);		// constructor of the Graph class
	~Graph();                  // destructor of the Graph class

	void addEdge(int x, int y, double w);	// add a new edge to Graph, A[x][y]=w
	void remEdge(int x, int y);	// remove an edge in Graph

	list<edge> maxfluxPath(int s, int t); // get the maxflux path between s and t by recursively getting subpaths
	bool connected(int s, int t, double weight); // check if s and t are connected through edges, all of which have a weight higher than given weight
	edge findBottleneck(int s,int t);  // find the bottleneck edge between s and t

}; 

Graph::Graph(int size) {		// constructor of Graph of size 
	int i, j;				// need two integers, i and j
	if (size < 2) n = 2;		// if the size is less than two elements, it's now two
	else n = size;			// otherwise, the size is size n
	
	A = new double*[n];			// A is a new array of size n

	for (i = 0; i < n; ++i)		// Build the first column of the A graph
		A[i] = new double[n];
	for (i = 0; i < n; ++i)		// Loop over the size
		for (j = 0; j < n; ++j)		// both dimensions
			A[i][j] = 0.0;		// initialize values to 0
	
} 

Graph::~Graph()  {			// destructor of Graph
	for (int i = 0; i < n; ++i)		// loop through rows of the graph
		delete [] A[i];
	delete [] A;			// delete graph
} 

void Graph::addEdge(int x, int y, double w) {		// add a new edge to the graph
	A[x][y] = w;		// by replacing the initialized 0 with w
} 

void Graph::remEdge(int x, int y) {		// remove edge from the graph 
	A[x][y] = 0.0;		 
} 

bool Graph::connected(int s, int t, double weight) 
{	// Use BFS algorithm to check connectivity between node s and node t.
	// The condition here is that connectivity is checked on only on the edges with weight greater than given "weight". 

	queue<int> Q;			 
	bool visited[n];	// declare the bool type pointer array declaring if a node has been visited
	int i,k;

	bool conn=false;        // by default not connected. At the end, this flag will tell whether the two nodes are connected or not. (with edges greater than given weight)

	for (i = 0; i < n; ++i)			// loop through the number of nodes
		visited[i] = false;		// declare the visitation state of each node to zero

	Q.push(s);				// create the queue for node x

	if(s == t)
	{	conn = true;
		return conn;			// if this is the final node, return
	}

	visited[s] = true;			// node s has now been visited

	while (!Q.empty())                      // while the Queue is not empty
	{			
		k = Q.front();		        // return the first element of the queue as integer k
		Q.pop();
	
		if(k == t)                      // if k is the desired node
		{			
			conn=true;
			return conn;		// if this is the final node, return
		}

		for (i = 0; i < n; ++i)		// for the rest of the nodes, i, neighbors of k. 
		{	
			
			if (A[k][i] >= weight && !visited[i]) 	// if node i is a neighbor of current node k, with edge weight greater than required weight, and has not been visited so far, then add it to the queue and visit it
			{	
				if(verticesOnPath.find(i)==verticesOnPath.end() || i==t) // don't explore edges to and from vertices that are already on the path. With the exception of the vertex t.
				{  
					Q.push(i);			// place i in the queue  
					visited[i] = true;		// and report that i has been visited
				}
			}
		}
	}
	
	return conn;
}

edge Graph::findBottleneck(int s,int t) 
{
	// Find bottleneck between the vertices s and t in current graph
	edge bottleneck;

	// Check if the last edge is a direct edge between s and t, then just return that edge as the bottleneck edge.
	if(edgeList.back().u==s and edgeList.back().v==t)
	{
		bottleneck=edgeList.back();
		edgeList.erase( edgeList.begin() + edgeList.size()-1 ); // erase the last edge
		A[bottleneck.u][bottleneck.v]=0.0;	// also remove from the graph
		
		return bottleneck;
	}

	int l = 0, r = edgeList.size()-1,m;

	while(r-l>1)
	{
		m = int(floor((r+l)/2)); 
		
		if(connected(s,t,edgeList[m].weight)) // check if the graph is connected between s and t with all connecting edges of weight > weight of the m'th edge
			l = m;
		else
			r = m;
	}

	bottleneck=edgeList[l];
	edgeList.erase( edgeList.begin() + l ); // need to return an iterator to erase
	A[bottleneck.u][bottleneck.v]=0.0;	// also remove from the graph

	// add vertices in the current list 
	verticesOnPath.insert(bottleneck.u);
	verticesOnPath.insert(bottleneck.v);

	return bottleneck; 
}

list<edge> Graph::maxfluxPath(int s, int t)
{
	edge bottleneck; // bottleneck edge between s and t (u,v)

	list<edge> pa,pb; // pa is the left part of the path from (s...u) and pb is the right part from (v...t), which is found recursively	

	if(s==t) // base case
	{	return(pa);  // return empty path, i.e. no edges
	}
	
	// Find the bottleneck between s and t in the current graph
	bottleneck=findBottleneck(s,t);    // find bottleneck between s and t
	
	// Decompose the current graph into two subgraphs, left of bottleneck and right of bottleneck
	if(s!=bottleneck.u)
	{
		pa=maxfluxPath(s,bottleneck.u);                 // Recursively obtain the path on the left subgraph

	}

	if(bottleneck.v!=t)
	{	pb=maxfluxPath(bottleneck.v,t);                 // Recursively obtain the path on the right subgraph

	}
	
	pa.push_back(bottleneck); // add bottleneck to the end of the list got from the first part of the path

	pa.insert(pa.end(),pb.begin(),pb.end()); // concatenate the second list to end of first one

	return(pa);		
}



int main(int argc, char *argv[])
{

	// ------------------------------------------------------
	//              Get input
	// ------------------------------------------------------
	char wtFilename[FILENAME_LEN];

	strcpy(wtFilename,argv[1]);
	numNodes=atoi(argv[2]); // total number of nodes in the graph
	start=atoi(argv[3]);  // vertex indices start from 1
	end=atoi(argv[4]); // vertex indices start from 1

	// translate all vertex numbers i, in the input, to vertex number i-1 in the program
	start = start -1;
	end = end -1;

	// Initialize graphs
	Graph g(numNodes);				// build the adjacency graph
	
	// return value for functions
	int retvalue;
	
	// ------------------------------------------------------
	//             Determine edges and weights in the graph
	// ------------------------------------------------------
	FILE *infile;
	int i = 0, j=0;
	infile = fopen(wtFilename,"r");

	// file opened properly?
	if(infile==NULL)
	{	printf("Could not properly load weight file!\n");
	return(ERR);
	}

	for(i=0;i<numNodes;i++)
		for(j=0;j<numNodes;j++)
		{
			fscanf(infile,"%lf",&g.A[i][j]);
		}

	fclose(infile);

	// ------------------------------------------------------
	//              Run modified Dijkstra's algorithm
	// ------------------------------------------------------
	// Get the timing information
	//long startTime=myclock();

	// Get the sorted list of edges.
	edge inp;

	for(i=0;i<numNodes;i++)
	{	for(j=0;j<numNodes;j++)
		{
			if(g.A[i][j] > 0.0)
			{
				inp.u=i;
				inp.v=j;
				inp.weight=g.A[i][j];
				edgeList.push_back(inp);
			}	
		}
	}

	sort(edgeList.begin(),edgeList.end(),edgeSorter);

	// This removes degeneracy of edges in the graph
	
	for(i=0;i<edgeList.size();i++) 
	{
		edgeList[i].weight=i+1; //uniformly spacing all edge weights from 1.....edgeList.size()
		
		g.A[edgeList[i].u][edgeList[i].v]=edgeList[i].weight;
	}	


	list<edge> m=g.maxfluxPath(start,end);	

	// Get the end time
	//long endTime=myclock();
	//cout << setprecision(15) << getRuntime(&endTime, &startTime)/1000.0 <<  endl ;	
	
	// Print output in Pajek format
		
	cout << "*Vertices " << numNodes << endl;
	
	for(i=0;i<numNodes;i++)
		cout << i+1  << "\"   c Blue" << endl;
	
	cout << "*Arcs" << endl; 
		
	for(list<edge>::iterator it=m.begin();it!=m.end();++it)
		cout << (*it).u+1 << " " << (*it).v+1 	<< " " << (*it).weight << "\"  c Black" << endl;
	
	
	return(0);

}

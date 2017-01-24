// maxflux.cpp
// Program will calculate the MaxFlux path through a network

// -------------------------------------------------------------
//              History
// -------------------------------------------------------------
// written by: steven m. kreuzer and shruthi viswanath 
#include<cstdlib>
#include <cstdio>
#include<cstring>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <vector>
#include <algorithm>
#include <queue>
#include <sys/time.h>

using namespace std;

#define FILENAME_LEN 150

#define ERR -1 // return code from functions for error message
#define FINE 0 // return code for successful completion
        
int numNodes;	// the number of nodes

typedef struct {
int u;
int v;
float weight;
} edge; 


vector<edge> maxfluxPathList;
// list of edges in the max flux path

bool edgeSorter(edge  a, edge  b)
{ 
	return a.weight<b.weight;
}



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

// -------------------------------------------------------------
//              Object Definitions
// -------------------------------------------------------------
  
class Graph {				// define the Graph class
    public:	
        int n;				// includes an 'n' integer
        float **A;			// and a pointer to a pointer called 'A'

        Graph(int size = 2);		// constructor of the Graph class
	Graph(const Graph *g);
        ~Graph();
	int conn;			// destructor of the Graph class
        void addEdge(int x, int y,float w);	// add a new edge to Graph
        void remEdge(int x, int y);	// rem a new edge to Graph
        bool BFS(int , int);		// perform the breadth-first search
}; 
 
Graph::Graph(int size) {		// constructor of Graph of size 
    int i, j;				// need two integers, i and j
    if (size < 2) n = 2;		// if the size is less than two elements, it's now two
    else n = size;			// otherwise, the size is size n
    A = new float*[n];			// A is a new array of size n
    for (i = 0; i < n; ++i)		// Build the first column of the A graph
        A[i] = new float[n];
    for (i = 0; i < n; ++i)		// Loop over the size
        for (j = 0; j < n; ++j)		// both dimensions
            A[i][j] = 0.0;		// initialize values to 0
} 

Graph::Graph(const Graph *gRef) {	
    int i,j;	 
    n=gRef->n;
    A = new float*[n];			// A is a new array of size n
    for (i = 0; i < n; ++i)		// Build the first column of the A graph
        A[i] = new float[n];
    for (i = 0; i < n; ++i)		// Loop over the size
        for (j = 0; j < n; ++j)		// both dimensions
            A[i][j] = gRef->A[i][j];		// initialize values to 0			 
 
} 
 
Graph::~Graph() {			// destructor of Graph
    for (int i = 0; i < n; ++i)		// loop through rows of the graph
    delete [] A[i];
    delete [] A;			// delete graph
} 
 
 
void Graph::addEdge(int x, int y,float w) {		// add a new edge to the graph
    A[x][y] = w;		// by replacing the initialized 0 with a 1
	//cout << A[x-1][y-1] << endl;
} 


void Graph::remEdge(int x, int y) {		// add a new edge to the graph
    A[x][y] = 0.0;		// by replacing the initialized 0 with a 1
} 


bool Graph::BFS(int x, int required) {		// General BFS algorithm to check connectivity between node x and node required

	queue<int> Q;					// instantiate a Queue named Q
	bool visited[n];		// declare the bool type pointer array declaring if a node has been visited
	int i,k;

	bool conn=false;

	for (i = 0; i < n; ++i)			// loop through the number of nodes
		visited[i] = false;				// declare the visitation state of each node to zero

	Q.push(x);				// create the queue for node x

	if(x == required)
	{	conn = true;
		return conn;			// if this is the final node, return
	}

	visited[x] = true;				// node x has now been visited

	while (!Q.empty()) {			// while the Queue is not empty
		k = Q.front();			// return the first element of the queue as integer k
		Q.pop();

		if(k == required){			// if k is the desired node
			conn=true;
			return conn;			// if this is the final node, return
		}

		for (i = 0; i < n; ++i)		// for the rest of the nodes, i
			if (A[k][i] > 0 && !visited[i]) {	// if node k is connected but i is not visited
				Q.push(i);				// place i in the queue for k
				visited[i] = true;			// and report that i has been visited
			}
	}
	
	return conn;
}
 
int main(int argc, char*argv[]) 
{
       // ------------------------------------------------------
       //              Get input
       // ------------------------------------------------------
	char wtFilename[FILENAME_LEN];

	strcpy(wtFilename,argv[1]); // weights file
	numNodes=atoi(argv[2]); // number of nodes in the network
	int start=atoi(argv[3]); // vertex indices numbered from 1
	int end=atoi(argv[4]); // vertex indices numbered from 1 
	

	// Initialize graphs
	Graph g(numNodes);	// build the adjacency graph
	
	// return value for functions
	int retvalue;
	
	// make vertex numbering zero based
	start=start-1;
	end=end-1;

	// ------------------------------------------------------
	//             Determine edges and weights in the input graph
	// ------------------------------------------------------
	FILE *infile;
	int i = 0, j=0;
	infile = fopen(wtFilename,"r");

	// file opened properly 
	if(infile==NULL)
	{	printf("Could not properly load weight file!\n");
		return(ERR);
	}

	// Add edges and weights to appropriate weight matrix element. 

	for(i=0;i<numNodes;i++)
		for(j=0;j<numNodes;j++)
		{
			fscanf(infile,"%f",&g.A[i][j]);
				
								
		}

	fclose(infile);

        // ------------------------------------------------------
        //          Run maxflux algorithm. 
        // ------------------------------------------------------

	// Get the timing information
	//long startTime=myclock();

	// Add non-zero edges to the edge list
	edge inp;
	
	vector<edge> edgeList; 
	// global list containing edges in sorted order. 
	// The minimum weight edge is deleted from here at each iteration, i.e. the minimum weight edge of the first path is deleted first, followed by minimum weight edge from second path etc. 

	for(i=0;i<numNodes;i++)
	{	for(j=0;j<numNodes;j++)
		{ 	if(g.A[i][j]>0.0)
			{
				inp.u=i;
				inp.v=j;
				inp.weight=g.A[i][j];
				edgeList.push_back(inp);				
			}

		}
	}

	// Sort edges first. this is a global list used in calculation of every path
	sort(edgeList.begin(),edgeList.end(),edgeSorter);

	if(g.BFS(start,end)) // if the graph is not connected at this stage, then the answer path will be 0's. 
	{	
		// new graph to check connectivity of current path. Once we have enough edges such that the start and end vertices are connected in the path, we can stop
		Graph pathGraph(numNodes);

		vector<edge>::iterator it=edgeList.begin();	// Go through the edge list

	
		while(!pathGraph.BFS(start,end)) // while the path graph between start and end vertices is missing edges
		{
			// temporarily remove current edge from graph and see if graph is still connected
			g.remEdge((*it).u,(*it).v); 

			if(!g.BFS(start,end))  // removal of this edge makes start and end vertices disconnected
			{     
				maxfluxPathList.push_back((*it)); // add this to the path list
				pathGraph.addEdge((*it).u,(*it).v,(*it).weight); // add this critical edge to the final path
				g.addEdge((*it).u,(*it).v,(*it).weight); // need to add back this edge to the graph for subsequent edge explorations! 

			}
			
			++it; // move to next edge in list
		}	
		
	}

	/*
	// Get the end time
	long endTime=myclock();
	cout << setprecision(15) << getRuntime(&endTime, &startTime)/1000.0 << endl ; // time in milliseconds
	*/ 		
	
	// Print output in Pajek format
	cout << "*Vertices " << numNodes << endl;
	
	for(i=0;i<numNodes;i++)
		cout << i+1  << "\"   c Blue" << endl;
	
	cout << "*Arcs" << endl; 
	for(vector<edge>::iterator it=maxfluxPathList.begin();it!=maxfluxPathList.end();++it)
		cout << (*it).u+1 << " " << (*it).v+1 	<< " " << (*it).weight << " c Black" << endl;
	

	return(0);

}

// dijkstraMaxflux.cpp
// Program that uses the modified Dijkstra algorithm to calculate the MaxFlux path through a network.
// Known inefficiencies:
// 1. Graph is represented by an adjacency matrix
// 2. Priority queue is implemented by an array
// -------------------------------------------------------------
//              History
// -------------------------------------------------------------
// written by: steven m. kreuzer and shruthi viswanath
// date written: 02/06/2013
// Please dont be judgemental about our code! 

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <list>
#include <queue>
#include <sys/time.h> 

using namespace std;

#define FILENAME_LEN 150

#define ERR -1 // return code from functions for error message
#define FINE 0 // return code for successful completion

#define INFINITY 9999999
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
float weight;
} edge; 


// -------------------------------------------------------------
//              Object Definitions
// -------------------------------------------------------------

class Graph {				// define the Graph class
public:
	int n;				// includes an 'n' integer
	float **A;			// and a pointer to a pointer called 'A', adjacency matrix

	edge** maxflux;

	Graph(int size = 2);		// constructor of the Graph class
	~Graph();                  // destructor of the Graph class

	void addEdge(int x, int y, float w);	// add a new edge to Graph, A[x][y]=w
	void remEdge(int x, int y);	// remove an edge in Graph


	void dijkstraMaxCapacity(int s); // calculate the bottleneck edges from s to all vertices

	bool dijkstraRun(int vertex); // check whether we ran dijkstraMaxCapacity fom vertex

	list<edge> maxfluxPath(int s, int t); // get the maxflux path between s and t by recursively getting subpaths

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

	// Similarly build the maxflux edges array
	maxflux = new edge*[n];			 

	for (i = 0; i < n; ++i)		 
		maxflux[i] = new edge[n];

	for (i = 0; i < n; ++i)		// Loop over the size
		for (j = 0; j < n; ++j)		// both dimensions
		{
			maxflux[j][j].u=0;
			maxflux[i][j].v=0;
			maxflux[i][j].weight=-1.0;
		}
	
} 

Graph::~Graph() {			// destructor of Graph
	for (int i = 0; i < n; ++i)		// loop through rows of the graph
		delete [] A[i];
	delete [] A;			// delete graph

	for (int i = 0; i < n; ++i)		 
		delete [] maxflux[i];
	delete [] maxflux;			 
} 

void Graph::addEdge(int x, int y, float w) {		// add a new edge to the graph
	A[x][y] = w;		// by replacing the initialized 0 with a 1
} 

void Graph::remEdge(int x, int y) {		// add a new edge to the graph
	A[x][y] = 0.0;		// by replacing the initialized 0 with a 1
} 


void Graph::dijkstraMaxCapacity(int s)
{
	//-------------------------------------------------------------
	//     Run a single iteration of Dijkstra's max flux path algorithm
	// -------------------------------------------------------------

	int i, currVertex, neibr;

	float maxValue, minValue;	// maximum flux estimate for each vertex
	bool visited[n];		// keeps track of whether a vertex is visited or not in the current function call 

	// Initialize the flux and predecessor arrays to default values.
	for(i=0;i<n;i++)
	{	visited[i]=false;
	}

	maxflux[s][s].weight=INFINITY; // start the algorithm at s
	currVertex=s;

	// Actual modified Dijkstra for max flux is inside this while loop
	while(!visited[currVertex])
	{			// while the Queue is not empty

		for(neibr=0;neibr<n;neibr++)
		{
			if(A[currVertex][neibr]>0.0) // (currVertex,neibr) is an edge in the graph
			{   minValue=min(maxflux[s][currVertex].weight,A[currVertex][neibr]);  // min weight edge on the path <s....currVertex,neibr>

			    if(maxflux[s][neibr].weight < minValue)     // Use d(v)=max(d(v),min(d(u),c(u,v)))
			    { 	
				maxflux[s][neibr].weight= minValue;

				if(A[currVertex][neibr]>minValue) // the current edge being examined is not the minimum weight edge on the path <s....currVertex,neibr>
				{
					maxflux[s][neibr].u=maxflux[s][currVertex].u;
					maxflux[s][neibr].v=maxflux[s][currVertex].v;
					maxflux[s][neibr].weight=maxflux[s][currVertex].weight;					
				}
				else // the current edge (currVertex,neibr) is the bottleneck edge on the the path <s....currVertex,neibr>
				{
					maxflux[s][neibr].u=currVertex;
					maxflux[s][neibr].v=neibr;
					maxflux[s][neibr].weight=A[currVertex][neibr];	

				} 
			    }
			}
		}

		visited[currVertex]=true; // mark current vertex as visited

		// Extract next vertex with maximum flux
		maxValue=0.0;
		for(i=0;i<n;i++) // Do a linear search over the vertex array
		{
			if(maxflux[s][i].weight>maxValue && !visited[i])
			{	maxValue=maxflux[s][i].weight;
			        currVertex=i;
			}
		}
	}

		
}

bool Graph::dijkstraRun(int vertex)
{
	if(maxflux[vertex][vertex].weight==INFINITY)
		return true;
	return false;
}	

list<edge> Graph::maxfluxPath(int s, int t)
{
	list<edge> pa,pb; 
	edge bdg; // bottleneck edge between s and t 

	if(s==t) // base case
	{
		return(pa);  
	}

	// run max capacity path from s if not already done so
	if(!dijkstraRun(s))
		dijkstraMaxCapacity(s);    // now all the bottleneck edges for paths from s to all vertices are stored

	bdg=maxflux[s][t]; // bottleneck edge

	pa=maxfluxPath(s,bdg.u);
	pb=maxfluxPath(bdg.v,t);
	
	pa.push_back(bdg); // add bdg to the end of the list got from the first part of the path
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
			fscanf(infile,"%f",&g.A[i][j]);
		}

	fclose(infile);
	
	
	// ------------------------------------------------------
	//              Run modified Dijkstra's algorithm
	// ------------------------------------------------------
	// Get the timing information
	//long startTime=myclock();

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

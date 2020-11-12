#define _CRT_SECURE_NO_WARNINGS
#define WIN32
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <omp.h> 
#include<time.h>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <climits>
#include<map>
using namespace std;
//Global variables to store num of threads 
int NUM_THREADS = 1;

//Hash边
vector<map<int, int>> Adj;

typedef unsigned int vid_t;
typedef unsigned int eid_t;
typedef struct {
	long n;
	long m;

	vid_t* adj;
	eid_t* num_edges;
	eid_t* eid;
} graph_t;


int vid_compare(const void* a, const void* b) {
	return (*(vid_t*)a - *(vid_t*)b);
}

long figuremax(long a, long b) {
	return a < b ? b : a;
}

static double timer() {
	//struct timeval tp;
	//gettimeofday(&tp, NULL);
	//return ((double)(tp.tv_sec) + tp.tv_usec * 1e-6);
	return 0;
}

int load_graph_from_file(char* filename, graph_t* g) {

	FILE* infp = fopen(filename, "r");
	long m = 0;
	vid_t u, v, t;

#pragma endregion

	if (infp == NULL) {
		fprintf(stderr, "Error: could not open inputh file: %s.\n Exiting ...\n", filename);
		exit(1);
	}

	fprintf(stdout, "Reading input file: %s\n", filename);

	double t0 = timer();

	//Read N and M
	fscanf(infp, "%ld %ld\n", &(g->n), &(g->m));
	printf("N: %ld, M: %ld \n", g->n, g->m);


	//Allocate space
	g->num_edges = (eid_t*)malloc((g->n + 1) * sizeof(eid_t));
	assert(g->num_edges != NULL);

#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		g->num_edges[i] = 0;
	}

	Adj.resize(g->n + 1);
#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		Adj[i].clear();
	}

	map<int, int>::iterator it1, it2;
	while (fscanf(infp, "%u %u \n", &u, &v) != EOF) {
		m++;
		//避免重复			
		if (Adj[u].find(v) == Adj[u].end()) {
			Adj[u][v] = 0;
			Adj[v][u] = 0;
			g->num_edges[u]++;
			g->num_edges[v]++;
		}
	}
	fclose(infp);

	if (m != g->m) {
		printf("Reading error: file does not contain %ld edges.\n", g->m);
		free(g->num_edges);
		exit(1);
	}

	m = 0;

	//需要n+1个点
	eid_t* temp_num_edges = (eid_t*)malloc((g->n + 1) * sizeof(eid_t));
	assert(temp_num_edges != NULL);

	temp_num_edges[0] = 0;

	for (long i = 0; i < g->n; i++) {
		m += g->num_edges[i];
		temp_num_edges[i + 1] = m;
	}

	//g->m is twice number of edges
	g->m = m;

	//Allocate space for adj
	g->adj = (vid_t*)malloc(m * sizeof(vid_t));
	assert(g->adj != NULL);

#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (long i = 0; i < g->n + 1; i++)
			g->num_edges[i] = temp_num_edges[i];

#pragma omp for schedule(static)
		for (long i = 0; i < m; i++)
			g->adj[i] = 0;
	}


	infp = fopen(filename, "r");
	if (infp == NULL) {
		fprintf(stderr, "Error: could not open input file: %s.\n Exiting ...\n", filename);
		exit(1);
	}

	//Read N and M
	fscanf(infp, "%ld %ld\n", &(g->n), &m);

#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		Adj[i].clear();
	}
	//Read the edges
	while (fscanf(infp, "%u %u\n", &u, &v) != EOF) {
		if (Adj[u].find(v) == Adj[u].end()) {
			Adj[u][v] = 0;
			Adj[v][u] = 0;
			g->adj[temp_num_edges[u]] = v;
			temp_num_edges[u]++;
			g->adj[temp_num_edges[v]] = u;
			temp_num_edges[v]++;
		}
	}

	fclose(infp);

	//Sort the adjacency lists
	for (long i = 0; i < g->n; i++) {
		qsort(g->adj + g->num_edges[i], g->num_edges[i + 1] - g->num_edges[i], sizeof(vid_t), vid_compare);
	}

	fprintf(stdout, "Reading input file took time: %.2lf sec \n", timer() - t0);
	free(temp_num_edges);
#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		Adj[i].clear();
	}
	return 0;
}


int load_graph_from_file2(char* filename, graph_t* g) {

	FILE* infp = fopen(filename, "r");
	long m = 0;
	vid_t u, v, t;

#pragma endregion

	if (infp == NULL) {
		fprintf(stderr, "Error: could not open inputh file: %s.\n Exiting ...\n", filename);
		exit(1);
	}

	fprintf(stdout, "Reading input file: %s\n", filename);

	double t0 = timer();

	//Read N and M
	fscanf(infp, "%ld %ld\n", &(g->n), &(g->m));
	printf("N: %ld, M: %ld \n", g->n, g->m);


	//Allocate space
	g->num_edges = (eid_t*)malloc((g->n + 1) * sizeof(eid_t));
	assert(g->num_edges != NULL);

#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		g->num_edges[i] = 0;
	}

	Adj.resize(g->n + 1);
#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		Adj[i].clear();
	}

	map<int, int>::iterator it1, it2;
	while (fscanf(infp, "%u %u \n", &u, &v) != EOF) {
		m++;
		//避免重复			
		if (Adj[u].find(v) == Adj[u].end()) {
			Adj[u][v] = 0;
			Adj[v][u] = 0;
			g->num_edges[u]++;
			g->num_edges[v]++;
		}
	}
	fclose(infp);

	if (m != g->m) {
		printf("Reading error: file does not contain %ld edges.\n", g->m);
		free(g->num_edges);
		exit(1);
	}

	m = 0;

	//需要n+1个点
	eid_t* temp_num_edges = (eid_t*)malloc((g->n + 1) * sizeof(eid_t));
	assert(temp_num_edges != NULL);

	temp_num_edges[0] = 0;

	for (long i = 0; i < g->n; i++) {
		m += g->num_edges[i];
		temp_num_edges[i + 1] = m;
	}

	//g->m is twice number of edges
	g->m = m;

	//Allocate space for adj
	g->adj = (vid_t*)malloc(m * sizeof(vid_t));
	assert(g->adj != NULL);

#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (long i = 0; i < g->n + 1; i++)
			g->num_edges[i] = temp_num_edges[i];

#pragma omp for schedule(static)
		for (long i = 0; i < m; i++)
			g->adj[i] = 0;
	}


	infp = fopen(filename, "r");
	if (infp == NULL) {
		fprintf(stderr, "Error: could not open input file: %s.\n Exiting ...\n", filename);
		exit(1);
	}

	//Read N and M
	fscanf(infp, "%ld %ld\n", &(g->n), &m);

#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		Adj[i].clear();
	}
	//Read the edges
	while (fscanf(infp, "%u %u\n", &u, &v) != EOF) {
		if (Adj[u].find(v) == Adj[u].end()) {
			Adj[u][v] = 0;
			Adj[v][u] = 0;
			g->adj[temp_num_edges[u]] = v;
			temp_num_edges[u]++;
			g->adj[temp_num_edges[v]] = u;
			temp_num_edges[v]++;
		}
	}

	fclose(infp);

	//Sort the adjacency lists
	for (long i = 0; i < g->n; i++) {
		qsort(g->adj + g->num_edges[i], g->num_edges[i + 1] - g->num_edges[i], sizeof(vid_t), vid_compare);
	}

	fprintf(stdout, "Reading input file took time: %.2lf sec \n", timer() - t0);
	free(temp_num_edges);
#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		Adj[i].clear();
	}
	return 0;
}


#pragma region 计算图文件节点

// 双边
void FigureEdgeAndNode2(char* filename) {
	FILE* infp = fopen(filename, "r");
	long m = 0;
	long _max = 0;
	vid_t u, v;
	// 检测节点
	while (fscanf(infp, "%u %u\n", &u, &v) != EOF) {
		m++;
		/*	g->num_edges[u]++;
			g->num_edges[v]++;*/
		_max = figuremax(_max, figuremax(u, v));
	}
	cout << "Max Node:" << _max << " Edges:" << m << endl;
}

// 三边
void FigureEdgeAndNode3(char* filename) {
	FILE* infp = fopen(filename, "r");
	long m = 0;
	long _max = 0;
	vid_t u, v, t;
	// 检测节点
	while (fscanf(infp, "%u %u %u\n", &u, &v, &t) != EOF) {
		m++;
		/*	g->num_edges[u]++;
			g->num_edges[v]++;*/
		_max = figuremax(_max, figuremax(u, v));
	}
	cout << "Max Node:" << _max << " Edges:" << m << endl;
}
#pragma endregion


int main() {


	char path[100] = "H:/竞赛相关/2020ccf-kmax/k-truss/filetest/example.txt";
	//char path[100] = "H:/竞赛相关/2020ccf-kmax/k-truss/k-truss/data/p2p-Gnutella31.txt";
	char path2[102] = "H:/竞赛相关/2020ccf-kmax/ktruss-data/s18.e16.rmat.edgelist.tsv";
	graph_t g;
	load_graph_from_file(path, &g);
	FILE* fp;
	fp = fopen("H:/竞赛相关/2020ccf-kmax/ktruss-data/s18.e16.rmat.edgelist.tsv", "r");
	int a, b, c;
	int ans = 0;
	while (fscanf(fp, "%d%d%d", &a, &b, &c) != EOF)
	{
		//cout << a << "-" << b << "-" << c << endl;
		ans++;
	}
	fclose(fp);
	cout << ans << endl;

	//for (int i = 0; i < 1000; i++)
	//{
	//	fscanf(fp, "%d%d%d", &a, &b, &c);
	//	cout << a << "-" << b << "-" << c << endl;
	//}
	system("pause");
}
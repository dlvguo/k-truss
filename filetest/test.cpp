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
using namespace std;
//Global variables to store num of threads 
int NUM_THREADS = 1;

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

static double timer() {
	//struct timeval tp;
	//gettimeofday(&tp, NULL);
	//return ((double)(tp.tv_sec) + tp.tv_usec * 1e-6);
	return 0;
}

int load_graph_from_file(char* filename, graph_t* g) {

	FILE* infp = fopen(filename, "r");
	if (infp == NULL) {
		fprintf(stderr, "Error: could not open inputh file: %s.\n Exiting ...\n", filename);
		exit(1);
	}

	fprintf(stdout, "Reading input file: %s\n", filename);

	double t0 = timer();

	//Read N and M
	fscanf(infp, "%ld %ld\n", &(g->n), &(g->m));
	long m = 0;
	vid_t u, v;
	vid_t _max = 0;

	//求到最大的max
	while (fscanf(infp, "%u %u\n", &u, &v) != EOF) {
		_max = max(_max, max(u, v));
	}

	fclose(infp);

	//计算max边长度
	printf("%d_max\n", _max);
	printf("N: %ld, M: %ld \n", g->n, g->m);


	infp = fopen(filename, "r");
	fscanf(infp, "%ld %ld\n", &(g->n), &(g->m));
	g->n = _max + 1;

	// Allocate space 分配空间
	g->num_edges = (eid_t*)malloc((g->n + 1) * sizeof(eid_t));
	assert(g->num_edges != NULL);

#pragma omp parallel for 
	for (long i = 0; i < g->n + 1; i++) {
		g->num_edges[i] = 0;
	}

	while (fscanf(infp, "%u %u\n", &u, &v) != EOF) {
		m++;
		g->num_edges[u]++;
		g->num_edges[v]++;
	}

	fclose(infp);

	// 避免边计算错误
	if (m != g->m) {
		printf("Reading error: file does not contain %ld edges.\n", g->m);
		free(g->num_edges);
		exit(1);
	}

	m = 0;

	eid_t* temp_num_edges = (eid_t*)malloc((g->n + 1) * sizeof(eid_t));
	assert(temp_num_edges != NULL);

	temp_num_edges[0] = 0;

	for (long i = 0; i < g->n + 1; i++) {
		g->num_edges[i] /= 2;
	}

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
	g->n = _max + 1;
	//Read the edges
	while (fscanf(infp, "%u %u\n", &u, &v) != EOF) {
		g->adj[temp_num_edges[u]] = v;
		temp_num_edges[u]++;
		g->adj[temp_num_edges[v]] = u;
		temp_num_edges[v]++;
	}

	/*for (long i = 0; i < g->n + 1; i++) {
		g->adj[i] /= 2;
	}*/


	fclose(infp);

	//Sort the adjacency lists
	for (long i = 0; i < g->n; i++) {
		qsort(g->adj + g->num_edges[i], g->num_edges[i + 1] - g->num_edges[i], sizeof(vid_t), vid_compare);
	}

	fprintf(stdout, "Reading input file took time: %.2lf sec \n", timer() - t0);
	free(temp_num_edges);
	cout << m << endl;
	return 0;
}


int main() {
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
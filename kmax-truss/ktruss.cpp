#define WIN32
#define  _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <omp.h> 
#include<time.h>
// #include <sys/time.h>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <climits>
#include<vector>
using namespace std;

typedef unsigned int UI;

// ϵͳ�߳��� Ĭ��Ϊ1
int ThreadNum = 1;


// ͼ�ṹ
typedef struct {
	// �ڵ�����
	long nodeNum;
	// �ܱ���*2 
	long m;
	// ��csr�������ɵı߱���
	UI* adj;
	// ������
	UI* num_edges;
	// ��ID
	UI* eid;
} Graph;

// �����߳���
void setting_thread() {
#pragma omp parallel
	{
#pragma omp master
		ThreadNum = omp_get_num_threads();
	}

	/*printf("NUM_PROCS:     %d \n", omp_get_num_procs());
	printf("NUM_THREADS:   %d \n", NUM_THREADS);*/
}

#pragma region Read File
//���� ���ֵ
UI figuremax(UI a, UI b) {
	return a < b ? b : a;
}

// ���ݱ�ID˳������
int vid_compare(const void* a, const void* b) {
	return (*(UI*)a - *(UI*)b);
}

// ��ͼ
void load_graph(Graph* g, char* path) {
	cout << path << endl;;
	FILE* infp = fopen(path, "r");
	if (infp == NULL) {
		fprintf(stderr, "Error: could not open inputh file: %s.\n Exiting ...\n", path);
		exit(1);
	}

	UI u, v, t;
	long m = 0, _max = 0;
	// fprintf(stdout, "Reading input file: %s\n", path);

	//TODO Ӧ�ÿ��Ը������ֶ�ȡ��ʽ �������ߺ���󶥵� 
	while (fscanf(infp, "%u %u %u\n", &u, &v, &t) != EOF) {
		m++;
		_max = figuremax(_max, figuremax(u, v));
	}
	fclose(infp);

	g->nodeNum = _max + 1;
	g->m = m;
	// cout << "N:" << g->n << " E:" << g->m << endl;

	//��ʼ����
	g->num_edges = (UI*)malloc((g->nodeNum + 1) * sizeof(UI));
	assert(g->num_edges != NULL);
#pragma omp parallel for 
	for (long i = 0; i < g->nodeNum + 1; i++) {
		g->num_edges[i] = 0;
	}


	// Adj.resize(g->n + 1);
	//#pragma omp parallel for 
	//	for (long i = 0; i < g->n + 1; i++) {
	//		Adj[i].clear();
	//	}


	infp = fopen(path, "r");
	//����
	while (fscanf(infp, "%u %u %u\n", &u, &v, &t) != EOF) {
		////�����ظ�			
		//if (Adj[u].find(v) == Adj[u].end()) {
		//	Adj[u][v] = 0;
		//	Adj[v][u] = 0;
		//	g->num_edges[u]++;
		//	g->num_edges[v]++;
		//}
		if (u > v) {
			g->num_edges[u]++;
			g->num_edges[v]++;
		}
	}

	fclose(infp);
	// ������һ�㲻��Ҫ
	if (m != g->m) {
		printf("Reading error: file does not contain %ld edges.\n", g->m);
		free(g->num_edges);
		exit(1);
	}

	m = 0;

	//��Ҫn+1���� ����һ��
	UI* temp_num_edges = (UI*)malloc((g->nodeNum + 1) * sizeof(UI));
	assert(temp_num_edges != NULL);

	temp_num_edges[0] = 0;

	//���Ż�
	for (long i = 0; i < g->nodeNum; i++) {
		m += g->num_edges[i];
		temp_num_edges[i + 1] = m;
	}

	//g->m is twice number of edges
	g->m = m;

	//Allocate space for adj
	g->adj = (UI*)malloc(m * sizeof(UI));
	assert(g->adj != NULL);

#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (long i = 0; i < g->nodeNum + 1; i++)
			g->num_edges[i] = temp_num_edges[i];

#pragma omp for schedule(static)
		for (long i = 0; i < m; i++)
			g->adj[i] = 0;
	}


	infp = fopen(path, "r");
	if (infp == NULL) {
		fprintf(stderr, "Error: could not open input file: %s.\n Exiting ...\n", path);
		exit(1);
	}

	//#pragma omp parallel for 
	//	for (long i = 0; i < g->n + 1; i++) {
	//		Adj[i].clear();
	//	}
		//Read the edges
	while (fscanf(infp, "%u %u %u\n", &u, &v, &t) != EOF) {
		/*if (Adj[u].find(v) == Adj[u].end()) {
			Adj[u][v] = 0;
			Adj[v][u] = 0;*/
		if (u > v) {
			g->adj[temp_num_edges[u]] = v;
			temp_num_edges[u]++;
			g->adj[temp_num_edges[v]] = u;
			temp_num_edges[v]++;
		}
		//}
	}

	fclose(infp);

	// ���ݱ�����
	for (long i = 0; i < g->nodeNum; i++) {
		qsort(g->adj + g->num_edges[i], g->num_edges[i + 1] - g->num_edges[i], sizeof(UI), vid_compare);
	}

	// fprintf(stdout, "Reading input file took time: %.2lf sec \n", timer() - t0);
	free(temp_num_edges);
}
#pragma endregion


int main(int argc, char* argv[]) {

	// fprintf(stderr, "Start Successfully \n");
	// TODO ע������߳����޸�
	if (argc < 2) {
		fprintf(stderr, "%s <Graph file>\n", argv[0]);
		exit(1);
	}
	// printf("FileFullName-%s\n", argv[1]);
	//char path[50] = "./example.txt";
	// ��ȡ�߳���
	setting_thread();
	Graph g;

	//��ȡͼ�ļ�
	load_graph(&g, argv[2]);

	//��ʱ
	//double t0 = timer();

	/************   Compute k - truss *****************************************/
	//��ȡ���б�
	/*Edge* edgeIdToEdge = (Edge*)malloc((g.m / 2) * sizeof(Edge));
	assert(edgeIdToEdge != NULL);*/

	//Populate the edge list array
	// getEidAndEdgeList(&g, edgeIdToEdge);

	int* EdgeSupport = (int*)calloc(g.m / 2, sizeof(int));
	assert(EdgeSupport != NULL);

#pragma omp parallel for 
	for (long i = 0; i < g.m / 2; i++) {
		EdgeSupport[i] = 0;
	}

	// PKT_intersection(&g, EdgeSupport, edgeIdToEdge);

	// ������
	// display_stats(EdgeSupport, g.m / 2);

	// fprintf(stdout, "Figure took time: %.2lf sec \n", timer() - t0);

	//Free memory
	// free_graph(&g);

	/*if (edgeIdToEdge != NULL)
		free(edgeIdToEdge);

	if (EdgeSupport != NULL)
		free(EdgeSupport);*/

	return 0;
}


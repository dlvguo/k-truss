#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <omp.h> 
#include<time.h>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <climits>
#include<vector>
#include<algorithm>
using namespace std;

typedef unsigned int UI;

// ϵͳ�߳��� Ĭ��Ϊ1
int ThreadNum = 1;

// ͼ�ṹ
struct Graph {
	// �ڵ�����
	long nodeNums;
	// �ܱ���*2 
	long m;
	// ��csr�������ɵı߱���
	UI* adj;
	// ������
	UI* num_edges;
	// ��ID
	UI* eid;
};

// ���ļ�������ͼ��
vector<UI> vts;

// �߽ṹ��
typedef struct Edge {
	UI u, v;
}El;

// �����߳���
void setting_thread() {
#pragma omp parallel
	{
#pragma omp master
		ThreadNum = omp_get_num_threads();
	}
}

#pragma region Read File
//���� ���ֵ
UI figuremax(UI a, UI b) {
	return a < b ? b : a;
}

// ��ͼ
void load_graph(Graph* g, char* path) {
	FILE* infp = fopen(path, "r");
	if (infp == NULL) {
		fprintf(stderr, "�Ҳ����ļ�: %s.\n Exiting ...\n", path);
		exit(1);
	}
	UI u, v, t, _max = 0;
	long m = 0;
	while (fscanf(infp, "%u %u %u\n", &u, &v, &t) != EOF) {
		//ֻ��ȡһ��
		if (u > v) {
			_max = figuremax(_max, u);
			vts.push_back(u);
			vts.push_back(v);
		}
	}
	fclose(infp);
	g->nodeNums = _max + 1;
	//��ʼ����
	g->num_edges = (UI*)malloc((g->nodeNums + 1) * sizeof(UI));

#pragma omp parallel for 
	for (long i = 0; i < g->nodeNums + 1; i++) {
		g->num_edges[i] = 0;
	}

	for (UI i = 0; i < vts.size(); i += 2)
	{
		g->num_edges[vts[i]]++;
		g->num_edges[vts[i + 1]]++;
	}

	//��Ҫn+1���� ����һ��
	UI* temp_num_edges = (UI*)malloc((g->nodeNums + 1) * sizeof(UI));
	temp_num_edges[0] = 0;

	//csr������� ���Ż�
	for (long i = 0; i < g->nodeNums; i++) {
		m += g->num_edges[i];
		temp_num_edges[i + 1] = m;
	}

	//g�µ�m��edges������
	g->m = m;
	g->adj = (UI*)malloc(m * sizeof(UI));

#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (long i = 0; i < g->nodeNums + 1; i++)
			g->num_edges[i] = temp_num_edges[i];

#pragma omp for schedule(static)
		for (long i = 0; i < m; i++)
			g->adj[i] = 0;
	}

	for (UI i = 0; i < vts.size(); i += 2)
	{
		g->adj[temp_num_edges[vts[i]]] = vts[i + 1];
		temp_num_edges[vts[i]]++;

		g->adj[temp_num_edges[vts[i + 1]]] = vts[i];
		temp_num_edges[vts[i + 1]]++;
	}
	vts.clear();
	free(temp_num_edges);
}


//��ȡ��ID�ͱ��б�Populate eid and edge list
void getEidAndEdgeList(Graph* g, Edge* els) {

	//����Eid��� ��СΪ�ܱ���*2 
	g->eid = (UI*)malloc(g->m * sizeof(UI));

	//ÿһ���ߵĳ�ʼ�ڵ����
	UI* num_edges_copy = (UI*)malloc(g->nodeNums * sizeof(UI));

	for (UI i = 0; i < g->nodeNums; i++) {
		num_edges_copy[i] = g->num_edges[i];
	}

	long edgeId = 0;

	//����ͼ����Ӧ��ֻ��Ҫһ��Գ� �������б߲�����ID
	for (UI u = 0; u < g->nodeNums; u++) {
		for (UI j = g->num_edges[u]; j < g->num_edges[u + 1]; j++) {
			UI v = g->adj[j];
			if (u < v) {
				Edge e;
				e.u = u;
				e.v = v;

				g->eid[j] = edgeId;
				num_edges_copy[u]++;

				if (g->adj[num_edges_copy[v]] == u) {
					g->eid[num_edges_copy[v]] = edgeId;
					num_edges_copy[v]++;
				}
				els[edgeId] = e;
				edgeId++;
			}
		}
	}
}

#pragma endregion

#pragma region ���м���edgeSupport��k-truss

void figureSubLevel_intersection(Graph* g, UI* curr, bool* InCurr, long currTail, int* edgeSupport,
	int level, UI* next, bool* InNext, long* nextTail, bool* processed, Edge* els) {

	//��������
	const long BUFFER_SIZE_BYTES = 4096;
	const long BUFFER_SIZE = BUFFER_SIZE_BYTES / sizeof(UI);

	UI buff[BUFFER_SIZE];
	long index = 0;

#pragma omp for schedule(dynamic,4)
	for (long i = 0; i < currTail; i++) {

		//�����
		UI e1 = curr[i];

		Edge edge = els[e1];

		UI u = edge.u;
		UI v = edge.v;

		UI uStart = g->num_edges[u], uEnd = g->num_edges[u + 1];
		UI vStart = g->num_edges[v], vEnd = g->num_edges[v + 1];

		UI numElements = (uEnd - uStart) + (vEnd - vStart);
		UI j_index = uStart, k_index = vStart;

		for (UI innerIdx = 0; innerIdx < numElements; innerIdx++) {
			if (j_index >= uEnd || k_index >= vEnd) {
				break;
			}
			else if (g->adj[j_index] == g->adj[k_index]) {
				UI e2 = g->eid[k_index];  //<v,w>
				UI e3 = g->eid[j_index];  //<u,w>
				//����� e1, e2, e3 ���һ��������
				if ((!processed[e2]) && (!processed[e3])) {
					if (edgeSupport[e2] > level && edgeSupport[e3] > level) {
						int supE2 = __sync_fetch_and_sub(&edgeSupport[e2], 1);
						if (supE2 == (level + 1)) {
							buff[index] = e2;
							InNext[e2] = true;
							index++;
						}
						if (supE2 <= level) {
							__sync_fetch_and_add(&edgeSupport[e2], 1);
						}
						if (index >= BUFFER_SIZE) {
							long tempIdx = __sync_fetch_and_add(nextTail, BUFFER_SIZE);
							for (long bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
								next[tempIdx + bufIdx] = buff[bufIdx];
							index = 0;
						}
						int supE3 = __sync_fetch_and_sub(&edgeSupport[e3], 1);
						if (supE3 == (level + 1)) {
							buff[index] = e3;
							InNext[e3] = true;
							index++;
						}
						if (supE3 <= level) {
							__sync_fetch_and_add(&edgeSupport[e3], 1);
						}
						if (index >= BUFFER_SIZE) {
							long tempIdx = __sync_fetch_and_add(nextTail, BUFFER_SIZE);
							for (long bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
								next[tempIdx + bufIdx] = buff[bufIdx];
							index = 0;
						}
					}
					else if (edgeSupport[e2] > level) {
						if (e1 < e3 && InCurr[e3]) {
							int supE2 = __sync_fetch_and_sub(&edgeSupport[e2], 1);
							if (supE2 == (level + 1)) {
								buff[index] = e2;
								InNext[e2] = true;
								index++;
							}
							if (supE2 <= level) {
								__sync_fetch_and_add(&edgeSupport[e2], 1);
							}
							if (index >= BUFFER_SIZE) {
								long tempIdx = __sync_fetch_and_add(nextTail, BUFFER_SIZE);
								for (long bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
									next[tempIdx + bufIdx] = buff[bufIdx];
								index = 0;
							}
						}
						if (!InCurr[e3]) {
							int supE2 = __sync_fetch_and_sub(&edgeSupport[e2], 1);
							if (supE2 == (level + 1)) {
								buff[index] = e2;
								InNext[e2] = true;
								index++;
							}
							if (supE2 <= level) {
								__sync_fetch_and_add(&edgeSupport[e2], 1);
							}
							if (index >= BUFFER_SIZE) {
								long tempIdx = __sync_fetch_and_add(nextTail, BUFFER_SIZE);
								for (long bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
									next[tempIdx + bufIdx] = buff[bufIdx];
								index = 0;
							}
						}
					}
					else if (edgeSupport[e3] > level) {
						if (e1 < e2 && InCurr[e2]) {
							int supE3 = __sync_fetch_and_sub(&edgeSupport[e3], 1);
							if (supE3 == (level + 1)) {
								buff[index] = e3;
								InNext[e3] = true;
								index++;
							}
							if (supE3 <= level) {
								__sync_fetch_and_add(&edgeSupport[e3], 1);
							}
							if (index >= BUFFER_SIZE) {
								long tempIdx = __sync_fetch_and_add(nextTail, BUFFER_SIZE);
								for (long bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
									next[tempIdx + bufIdx] = buff[bufIdx];
								index = 0;
							}
						}
						if (!InCurr[e2]) {
							int supE3 = __sync_fetch_and_sub(&edgeSupport[e3], 1);
							if (supE3 == (level + 1)) {
								buff[index] = e3;
								InNext[e3] = true;
								index++;
							}
							if (supE3 <= level) {
								__sync_fetch_and_add(&edgeSupport[e3], 1);
							}
							if (index >= BUFFER_SIZE) {
								long tempIdx = __sync_fetch_and_add(nextTail, BUFFER_SIZE);
								for (long bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
									next[tempIdx + bufIdx] = buff[bufIdx];
								index = 0;
							}
						}
					}
				}
				j_index++;
				k_index++;
			}
			else if (g->adj[j_index] < g->adj[k_index]) {
				j_index++;
			}
			else if (g->adj[k_index] < g->adj[j_index]) {
				k_index++;
			}
		}
	}
	if (index > 0) {
		long tempIdx = __sync_fetch_and_add(nextTail, index);;
		for (long bufIdx = 0; bufIdx < index; bufIdx++)
			next[tempIdx + bufIdx] = buff[bufIdx];
	}

#pragma omp barrier

#pragma omp for schedule(static)
	for (long i = 0; i < currTail; i++) {
		UI e = curr[i];
		processed[e] = true;
		InCurr[e] = false;
	}
#pragma omp barrier
}

// ɨ���Ϊlevel�ı�
void scanLevel(long edgeNums, int* edgeSupport, int level, UI* curr, long* currTail, bool* inCurr) {
	// �������� ���Ǹ���
	const long BUFFER_SIZE_BYTES = 4096;
	const long BUFFER_SIZE = BUFFER_SIZE_BYTES / sizeof(UI);

	UI buff[BUFFER_SIZE];
	long index = 0;

#pragma omp for schedule(static) 
	for (long i = 0; i < edgeNums; i++) {
		if (edgeSupport[i] == level) {
			// ����level�ı�
			buff[index] = i;
			inCurr[i] = true;
			index++;

			if (index >= BUFFER_SIZE) {
				// ʹ��ԭ�Ӳ��� ���û����� ����ռ�ù�����
				long tempIdx = __sync_fetch_and_add(currTail, BUFFER_SIZE);
				for (long j = 0; j < BUFFER_SIZE; j++) {
					curr[tempIdx + j] = buff[j];
				}
				index = 0;
			}
		}
	}

	if (index > 0) {
		long tempIdx = __sync_fetch_and_add(currTail, index);
		for (long j = 0; j < index; j++) {
			curr[tempIdx + j] = buff[j];
		}
	}
#pragma omp barrier
}

// ����ktruss
void figurektruss(Graph* g, int* edgeSupport, Edge* els) {
	// ������
	long edgeNums = g->m / 2;
	long nodeNums = g->nodeNums;
	// ��ǰIndex ����һ��Index
	long currTail = 0, nextTail = 0;
	// �жϱ��Ƿ���
	bool* processed = (bool*)malloc(edgeNums * sizeof(bool));
	// ��ǰ������Ҫ���ʵĽڵ�
	UI* curr = (UI*)malloc(edgeNums * sizeof(UI));
	// �жϵ�ǰ�ڵ��Ƿ����
	bool* InCurr = (bool*)malloc(edgeNums * sizeof(bool));
	// �´η�����Ҫ����������
	UI* next = (UI*)malloc(edgeNums * sizeof(UI));
	// ���̷߳�ֹ�´η��ʵ��ظ��Ľڵ�
	bool* InNext = (bool*)malloc(edgeNums * sizeof(bool));
	// ��ʼ����ı�
	UI* startEdge = (UI*)malloc(nodeNums * sizeof(UI));

#pragma region ����ģ��
#pragma omp parallel 
	{
		int tid = omp_get_thread_num();

		// �ж��ھ��Ƿ񱻷���
		UI* X = (UI*)malloc(g->nodeNums * sizeof(UI));
		for (UI i = 0; i < g->nodeNums; i++) {
			X[i] = 0;
		}

#pragma omp for schedule(static) 
		for (UI e = 0; e < edgeNums; e++) {
			processed[e] = false;
			InCurr[e] = false;
			InNext[e] = false;
		}

#pragma omp for schedule(static) 
		for (UI i = 0; i < nodeNums; i++) {
			UI j = g->num_edges[i];
			UI endIndex = g->num_edges[i + 1];

			while (j < endIndex) {
				if (g->adj[j] > i)
					break;
				j++;
			}
			// �ҵ���С�����ĵ�
			startEdge[i] = j;
		}
#pragma omp for schedule(dynamic,10) 
		for (UI u = 0; u < nodeNums; u++) {

			for (UI j = startEdge[u]; j < g->num_edges[u + 1]; j++) {
				UI w = g->adj[j];
				X[w] = j + 1;
			}
			for (UI j = g->num_edges[u]; j < startEdge[u]; j++) {
				UI v = g->adj[j];

				for (UI k = g->num_edges[v + 1] - 1; k >= startEdge[v]; k--) {
					UI w = g->adj[k];
					// ������
					if (w <= u) {
						break;
					}
					if (X[w]) {
						//�ж�����һ��������
						//��ID: <u,w> : g->eid[ X[w] -1] 
						//<u,w> : g->eid[ X[w] -1] 
						//<v,u> : g->eid[ j ]  
						//<v,w> : g->eid[ k ]		
						UI e1 = g->eid[X[w] - 1], e2 = g->eid[j], e3 = g->eid[k];
						__sync_fetch_and_add(&edgeSupport[e1], 1);
						__sync_fetch_and_add(&edgeSupport[e2], 1);
						__sync_fetch_and_add(&edgeSupport[e3], 1);
					}
				}
			}
			// ��0
			for (UI j = startEdge[u]; j < g->num_edges[u + 1]; j++) {
				UI w = g->adj[j];
				X[w] = 0;
			}
		}
#pragma omp barrier

		//����truss
		int level = 0;
		long todo = edgeNums;

		while (todo > 0) {
			// ɨ��Level
			scanLevel(edgeNums, edgeSupport, level, curr, &currTail, InCurr);
			while (currTail > 0) {
				todo = todo - currTail;
				// �ֽ��ӽṹ
				figureSubLevel_intersection(g, curr, InCurr, currTail, edgeSupport, level, next, InNext, &nextTail, processed, els);
				if (tid == 0) {
					UI* tempCurr = curr;
					curr = next;
					next = tempCurr;
					bool* tempInCurr = InCurr;
					InCurr = InNext;
					InNext = tempInCurr;
					currTail = nextTail;
					nextTail = 0;
				}
#pragma omp barrier	
			}
			level = level + 1;
#pragma omp barrier	

		}
		free(X);
	}

#pragma endregion
	//�ͷ��ڴ�
	free(next);
	free(InNext);
	free(curr);
	free(InCurr);
	free(processed);
	free(startEdge);
}

#pragma endregion


// ���Kmax
void display_kmaxtruss(int* edgeSupport, long edgeNums) {
	int maxSup = 0;
	for (long i = 0; i < edgeNums; i++) {

		if (maxSup < edgeSupport[i]) {
			maxSup = edgeSupport[i];
		}
	}

	long numEdgesWithMaxSup = 0;

	for (long i = 0; i < edgeNums; i++) {
		if (edgeSupport[i] == maxSup) {
			numEdgesWithMaxSup++;
		}
	}
	printf("kmax = %d, Edges in kmax-truss = %ld.\n", maxSup + 2, numEdgesWithMaxSup);
}

int main(int argc, char* argv[]) {

	if (argc < 2) {
		fprintf(stderr, "%s �Ҳ�������·��\n", argv[0]);
		exit(1);
	}
	// ��ȡ�߳���
	setting_thread();
	if (argc > 4) {
		int tempThreadNum = atoi(argv[4]);
		ThreadNum = tempThreadNum < ThreadNum ? tempThreadNum : ThreadNum;
	}
	Graph g;

	//��ȡͼ�ļ�
	load_graph(&g, argv[2]);

	//��ȡ���б�
	Edge* els = (Edge*)malloc((g.m / 2) * sizeof(Edge));

	//��ȡ��list
	getEidAndEdgeList(&g, els);

	//��֧�ֶ�
	int* edgeSupport = (int*)calloc(g.m / 2, sizeof(int));

#pragma omp parallel for 
	for (long i = 0; i < g.m / 2; i++) {
		edgeSupport[i] = 0;
	}

	//����ktruss
	figurektruss(&g, edgeSupport, els);

	// ���kmax
	display_kmaxtruss(edgeSupport, g.m / 2);

	//Free memory
	// free_graph(&g);

	if (els != NULL)
		free(els);

	if (edgeSupport != NULL)
		free(edgeSupport);

	return 0;
}


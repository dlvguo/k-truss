#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <stack>
#include <time.h>
#include <fstream>
#include <stdlib.h>
using namespace std;

struct triNode
{
	int st, mid, ed;
	double PrSt, PrEd;
	triNode(int s, int m, int e, double ps, double pe) {
		st = s; mid = m; ed = e;
		PrSt = ps; PrEd = pe;
	}
};

struct boundData
{
	int up_bound;
	int low_bound;
};

struct NeiNode {
	int id;
	NeiNode* next;
	// NeiNeiNode * firstNeiNei;
	int valid = 1;
	int sup;
	double PrSup = 0;
	int SS = 0;
	int CS = 0;
	int tuss = 0;
	double Pr = 0.1;
	//int eeMap = 0;
	vector<triNode> tris;
	NeiNode() { Pr = 0.5; }
};

struct GraphNode
{
	NeiNode* firstNei;
	int id;
	int degree = 0;
};

struct edge
{
	int st, ed;
	int st2, ed2;
	int sup;
	edge(int s, int e, int p) { st = s; ed = e; sup = p; }
	edge() {}
};

struct cmp {
	bool operator ()(const edge a, const edge b) {
		if (a.sup != b.sup)
			return a.sup < b.sup;
		else if (a.st != b.st)
			return a.st < b.st;
		return a.ed < b.ed;
	}
};


struct cmp2 {
	bool operator ()(const edge a, const edge b) {
		if (a.st != b.st)
			return a.st < b.st;
		return a.ed < b.ed;
	}
};

class Graph {
public:
	Graph();
	void buildGraph(int edge_num, int node_num);
	int getEdgeNum() { return edgeNum; }
	int getNodeNum() { return nodeNum; }
	void initSup();
	void addEdge(int stId, int edId);
	void addEdge(int stId, int edId, double pr);
	void output();
	void egdeBinSort();
	void greed();
	void distribute();
	void writeFile(char* filename);
	void initPrSup();
	void PrGreed();
	void PrDistribute();
	void dynamicInsert(int stId, int edId);
	void dynamicDelete(int stId, int edId);
	void supInitDelete(int stId, int edId);
	void centerInsert(int stId, int edId, int model);
	void centerDelete(int stId, int edId, int model);
	void centerMultInsert(vector<int> stIds, vector<int > edIds);
	void centerMultDelete(vector<int> stIds, vector<int > edIds);
	void cover();//用sup的值来覆盖truss
	void enableAllNodes();
	void initSuperSup();
	void initConstrainSup();
	void SingleNodeInsert(set<pair<int, int> > edgs, int nodeId);
	void MultNodeInsert(vector<set<pair<int, int> > > edgs, vector<int > nodeIds);
	int changedEdgeNum = 0;
	int changedNodeNum = 0;
	int totalSteps = 0;
	double totalTime = 0;
	void outputDynamicInfo(int computeType) {
		cout << "Compute mode:" << computeType << endl;
		cout << "Total Affected nodes number:" << changeNodes.size() << endl;
		cout << "The depth of broadcast:" << maxDepth << endl;
		cout << "Total steps:" << totalSteps << endl;
		cout << "Total Message Number:" << totalMsgNumber << endl;
	}
	int startCntSteps = 0;
	string logfile = "log.txt";
	void log(int cycleNumber, int number, double time) {
		FILE* fp;//定义一个文件指针
		fp = fopen(logfile.c_str(), "a");//以只写的方式打开文件，前面的参数是文件路径，后面的参数是表示只写
		if (fp) {
			fprintf(fp, "Parameter Message : %s %d %d %d %d \n", filename, method, graphType, computeType, enumber);
			fprintf(fp, "Number of cycles: %d \n", cycleNumber);
			fprintf(fp, "Max inserted edges number of all cycles: %d \n", number);
			fprintf(fp, "Time: %f mm\n", time / CLOCKS_PER_SEC);
			fprintf(fp, "--------------------------------\n\n");
		}
		else cout << "Error occurred when open the file!" << endl;
		fclose(fp);
	}
	char* filename;//文件路径
	int method = -1;//求解方法 0-贪心 1-分布式
	int graphType = -1;//图类型：0-static graph  1-temporal graph
	int computeType = -1;//0-不插入不删除 1-插入：逐条初始化并分布式计算 2-插入：批次插入初始化并分布式计算 3-插入：批次插入，初始化为Sup+2并分布式计算 4-删除：批次删除并分布式计算 5-删除：初始化min{sup+2，truss(e)}并分布式计算
	int enumber = -1; //插入/删除 数量级
	void recoard(char* file, int m, int g, int c, int n) {
		filename = file;
		method = m;
		graphType = g;
		computeType = c;
		enumber = n;
	}


private:
	int edgeNum;
	int nodeNum;
	double gamma = 0.1;
	GraphNode* nodeList;
	map<int, int> findId;
	//map<pair<int,int>,edge> findInfo; 
	set<edge, cmp> edgeSet;
	int Find(int id);
	int cnt = 0;
	vector<pair<int, int> >* bin;
	map<pair<int, int>, int> changeEdge;//动态插入时边改变的数
	set<pair<int, int> > cntEdge;//需要改变的边
	set<int > changeNodes;//涉及到的点
	map<pair<int, int>, int> visit;//动态插入时边改变的数
	int maxDepth = -1;
	int totalMsgNumber = 0;
	map<pair<int, int>, int> Vset;
	map<pair<int, int>, int> Xset;
	map<pair<int, int>, int > Sset;
	stack<pair<int, int> > Stack;
	clock_t startTime;
	clock_t endTime;

	void setEdgeMess(int stId, int edId, int k, int initModle);
	void remEdge(int stId, int edId);
	void PrRemEdge(int stId, int edId);
	double computePrSup(int st, int ed, double p);
	void deleteEdge(int st, int ed);
	void outputSet();
	int computeSup(int st, int ed);
	boundData bound(int st, int ed);
	void upAdjust(int st, int ed);//动态添边  记录递归深度
	void upAdjust_2(int st, int ed, int deep);//动态添边扩散
	void downAdjust(int st, int ed);//动态删边
	int getEdgeMess(int stId, int edId, int model);
	bool inCntEdge(int st, int ed);
	int computeSuperSup(int st, int ed, int truss);
	int computeConstrainSup(int st, int ed, int truss);
	void Eliminate(pair<int, int> e, int t);
	int computeLB(int st, int ed);
	int computeUB(int st, int ed);
	void centerAdjust(set<edge, cmp> PES);
	int MapToIndex(int id);
	bool exist_edge(int st, int ed);
	bool noCross(vector<edge > eset, int st, int ed);
	bool noNodeCorss(vector<set<pair<int, int> > > edgs, set<pair<int, int> > edg);
	set<edge, cmp> oneNodePES(set<pair<int, int> > edgs, int nodeId);
	void DeleteCenterAdjust(set<edge, cmp> PES);

};

/*
struct NeiNeiNode
{
	int id;
	int sup;
	NeiNeiNode *next;
};
*/


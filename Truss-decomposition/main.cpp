#include "Graph.h"
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <map>
#define random(x) rand()%x;
using namespace std;

int seed = 100;

Graph G;
int node_num, edge_num;//边数，点数

string temp;//文件输出临时存储字符串

//在进行随机插入or删除的时候，需要进行抽样
int changeNum;//抽样数目
map<int, int > rdm; //用于随机抽样,map去重
vector<int > shhuffle_rdm;//用于随机抽样,进行shuffle随机排列
vector<int> v;//保存抽样出的数据的下标
vector<int > stId;//要插入/删除边的起点
vector<int > edId;//要插入/删除边的终点
set<pair<int, int> > clearSet;
set<int > clearNodeSet;

char* filename;//文件路径
int method = -1;//求解方法 0-贪心 1-分布式
int graphType = -1;//图类型：0-static graph  1-temporal graph
int computeType = -1;//0-不插入不删除 1-插入：逐条初始化并分布式计算 2-插入：批次插入初始化并分布式计算 3-插入：批次插入，初始化为Sup+2并分布式计算 4-删除：批次删除并分布式计算 5-删除：初始化min{sup+2，truss(e)}并分布式计算
int number = -1; //插入/删除 数量级
char write;  //'w'-将truss值写入"output.txt"文件  

int Pow(int k) {
	int res = 1;
	while (k--) res = res * 10;
	return res;
}

string getName(char* filename) {
	int len = strlen(filename);
	int st = -1, ed = -1;
	int fg1 = 0, fg2 = 0;
	for (int i = len - 1; i >= 0; i--) {
		if (filename[i] == '.' && !fg1) {
			fg1 = 1;
			ed = i - 1;
		}
		if ((filename[i] == '/' || filename[i] == '\\') && !fg2) {
			fg2 = 1;
			st = i + 1;
		}
	}
	if (st == -1) st = 0;
	string res = "";
	for (int i = st; i <= ed; i++) res += filename[i];
	return res;
}

int main(int argc, char* argv[]) {
	//strcat(filename,argv[1]);
	filename = argv[1];
	method = argv[2][0] - '0';
	graphType = argv[3][0] - '0';
	computeType = argv[4][0] - '0';
	number = argv[5][0] - '0';
	write = argv[6][0];
	seed = atoi(argv[7]);
	G.recoard(filename, method, graphType, computeType, number);

	ifstream in(filename);
	if (!in.is_open()) cout << "fail to open file!\n" << endl;

	if (argv[2][0] == '2' || argv[2][0] == '3') {//概率图
		getline(in, temp);
		istringstream strs(temp);
		int node_num;
		int node;
		strs >> node_num >> node;
		G.buildGraph(0, node_num);
	}

	//读入文件头，确定点数和边数，建图
	while (getline(in, temp)) {
		istringstream str(temp);
		if (temp[0] == '#') {
			if (temp[2] == 'N' && temp[3] == 'o' && temp[4] == 'd') {
				string t1, t2, t3;
				str >> t1 >> t2 >> node_num >> t3 >> edge_num;
				//cout<<temp<<endl;
				//cout<<t1<<" "<<t2<<" "<<node_num<<" "<<t3<<" "<<edge_num<<"  ss"<<endl;
			}
			else continue;
		}
		else break;
	}
	do {
		istringstream str(temp);
		if (temp[0] != '#') {
			int st, ed;
			str >> st >> ed;
			int a = min(st, ed);
			int b = max(st, ed);
			clearSet.insert(make_pair(a, b));
			clearNodeSet.insert(a);
			clearNodeSet.insert(b);
			//cout<<a<<" "<<b<<" "<<clearSet.size()<<endl;
		}
	} while (getline(in, temp));

	edge_num = clearSet.size();
	G.buildGraph(edge_num, node_num);

	changeNum = number;
	changeNum = Pow(changeNum);
	int RandomUpBound = 0;
	if (graphType == 3) RandomUpBound = node_num;
	else RandomUpBound = edge_num;

	srand(seed);
	if (graphType == 0 || graphType == 3) {
		for (int i = 0; i < RandomUpBound; i++) shhuffle_rdm.push_back(i);
		random_shuffle(shhuffle_rdm.begin(), shhuffle_rdm.end());
		for (int i = 0; i < changeNum; i++) v.push_back(shhuffle_rdm[i]);
		if (!v.empty()) sort(v.begin(), v.end());
	}
	else if (graphType == 1) {
		for (int i = 0; i < changeNum; i++) {
			v.push_back(RandomUpBound - changeNum + i);
		}
	}
	else if (graphType == 2) {
		v.push_back(rand() % RandomUpBound);
	}

	//for (int i = 0; i < v.size(); i++) 		cout<<v[i]<<endl;

	if (graphType == 3) {
		set<int>::iterator  it;
		int vcnt = 0;
		int cnt = 0;
		for (it = clearNodeSet.begin(); it != clearNodeSet.end(); it++) {
			if (v[vcnt] == cnt)
			{
				v[vcnt++] = (*it);
			}
			cnt++;
		}
		vector<set<pair<int, int> > > edgeFromNodeToBeInserted;
		for (int i = 0; i < v.size(); i++) {
			set<pair<int, int> > tempset;
			edgeFromNodeToBeInserted.push_back(tempset);
		}
		set<pair<int, int> >::iterator  iter = clearSet.begin();
		//cout<<1<<endl;
		while (iter != clearSet.end()) {
			int st = (*iter).first;
			int ed = (*iter).second;
			int flg = 0;
			for (int i = 0; i < v.size(); i++) {
				if (st == v[i] || ed == v[i]) {
					edgeFromNodeToBeInserted[i].insert(make_pair((*iter).first, (*iter).second));
					flg = 1;
					break;
				}
			}
			if (flg) {
				set<pair<int, int> >::iterator  tempIt = iter;
				iter++;
				clearSet.erase(tempIt);
			}
			else iter++;
		}
		//cout<<2<<endl;
		for (iter = clearSet.begin(); iter != clearSet.end(); iter++) {
			int st = (*iter).first;
			int ed = (*iter).second;
			G.addEdge(st, ed);
		}

		if (method == 0) {
			G.initSup();
			G.cover();
			G.greed();
		}
		else if (method == 1) {
			G.initSup();
			G.cover();
			G.distribute();
		}
		else cout << "Wrong modle." << endl;

		G.enableAllNodes();

		switch (computeType)
		{
		case 0:
			for (int i = 0; i < edgeFromNodeToBeInserted.size(); i++) {
				G.SingleNodeInsert(edgeFromNodeToBeInserted[i], v[i]);
			}
			G.log(changeNum, changeNum, G.totalTime);
			break;

		case 1:
			G.MultNodeInsert(edgeFromNodeToBeInserted, v);
			break;

		default:
			break;
		}

	}
	else {

		int rNum = 0;
		int pNum = 0;
		set<pair<int, int> >::iterator  iter;
		for (iter = clearSet.begin(); iter != clearSet.end(); iter++) {
			int needInsert = 0;
			int needDelete = 0;
			if (((computeType <= 3 && computeType >= 1) || computeType == 6 || computeType == 7 || computeType == 8 || computeType == 9) && !v.empty() && v[pNum] == rNum && pNum < changeNum) {
				//cout<<rNum<<" "<<pNum<<endl;
				pNum++;
				needInsert = 1;
				if (computeType == 8 || computeType == 9) {
					needDelete = 1;
					needInsert = 0;
				}
			}
			int st = (*iter).first;
			int ed = (*iter).second;
			//if(needInsert) cout<<st<<" "<<ed<<endl;;

			if (!needInsert) {
				G.addEdge(st, ed);
			}
			if (needInsert || computeType == 4 || computeType == 5 || needDelete) {//如果需要插入或者删除
				stId.push_back(min(st, ed));
				edId.push_back(max(st, ed));
			}
			rNum++;
			//cout<<rNum<<endl;
		}

		if (method == 0) {
			G.initSup();
			G.cover();
			G.greed();
		}
		else if (method == 1) {
			G.initSup();
			G.cover();
			G.distribute();
		}
		else cout << "Wrong modle." << endl;

		G.enableAllNodes();
		G.startCntSteps = 1;
		switch (computeType)
		{
		case 0:
			break;
		case 1:
			for (int i = 0; i < changeNum; i++) {
				G.dynamicInsert(stId[i], edId[i]);
				G.distribute();
			}
			break;
		case 2:
			for (int i = 0; i < changeNum; i++) G.dynamicInsert(stId[i], edId[i]);
			G.distribute();
			break;
		case 3:
			for (int i = 0; i < changeNum; i++) G.addEdge(stId[i], edId[i], 0);
			G.initSup();
			G.cover();
			G.distribute();
			break;
		case 4:
			for (int i = 0; i < changeNum; i++) {
				G.dynamicDelete(stId[i], edId[i]);
				G.distribute();
			}
			break;
		case 5:
			for (int i = 0; i < changeNum; i++) {
				G.supInitDelete(stId[i], edId[i]);
			}
			G.distribute();
			break;
		case 6:
			for (int i = 0; i < changeNum; i++) {
				G.centerInsert(stId[i], edId[i], 0);
			}
			G.log(changeNum, changeNum, G.totalTime);
			break;
		case 7:
			G.centerMultInsert(stId, edId);
			break;
		case 8:
			for (int i = 0; i < changeNum; i++) {
				G.centerDelete(stId[i], edId[i], 0);
			}
			G.log(changeNum, changeNum, G.totalTime);
			break;
		case 9:
			//cout<<stId.size()<<" "<<edId.size()<<endl;
			G.centerMultDelete(stId, edId);
			break;
		default:
			if (graphType != 2)
				cout << "Wrong modle!" << endl;
			break;
		}
		if (computeType != 0 && graphType != 2 && computeType != 6 && computeType != 7 && computeType != 8 && computeType != 9) G.outputDynamicInfo(computeType);

	}

	if (write == 'w') {
		string str = getName(filename);
		str = argv[5][0] + str;
		//cout<<str<<endl;
		G.writeFile((char*)str.c_str());
	}
}
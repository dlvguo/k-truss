#include "Graph.h"
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stack>
#include <algorithm>
#include <stdlib.h>

Graph::Graph() {}
void Graph::buildGraph(int edge_num, int node_num) {
	edgeNum = 0;
	nodeNum = node_num;
	nodeList = new GraphNode[node_num + 10]();
	bin = new vector<pair<int, int> >[node_num - 1]();
}

int Graph::Find(int id) {
	map<int, int>::iterator iter;
	iter = findId.find(id);
	if (iter != findId.end())
		return iter->second;
	else
		return -1;
}

int Graph::MapToIndex(int stId) {
	int findSt = Find(stId);
	if (findSt < 0) {
		findSt = cnt;
		nodeList[cnt].id = stId;
		findId[stId] = cnt++;
	}
	return findSt;
}

void Graph::addEdge(int stId, int edId) {
	int findSt = MapToIndex(stId);
	int findEd = MapToIndex(edId);
	NeiNode* p = nodeList[findSt].firstNei;//无重边
	while (p) {
		if (findEd == p->id && p->valid) {
			cout << "Wrong model:" << stId << " " << edId << endl;
			return;
		}
		p = p->next;
	}
	NeiNode* e = new NeiNode;
	e->id = findEd;
	e->next = nodeList[findSt].firstNei;
	e->sup = 0;
	nodeList[findSt].firstNei = e;
	nodeList[findSt].degree++;

	e = new NeiNode;
	e->id = findSt;
	e->next = nodeList[findEd].firstNei;
	e->sup = 0;
	nodeList[findEd].firstNei = e;
	nodeList[findEd].degree++;
	visit[make_pair(min(findSt, findEd), max(findSt, findEd))] = 0;

	edgeNum++;
}

void Graph::addEdge(int stId, int edId, double pr) {
	//edgeNum ++;
	if (stId == edId) {//无自环
		return;
	}
	int findSt = Find(stId);
	if (findSt < 0) {
		findSt = cnt;
		nodeList[cnt].id = stId;
		findId[stId] = cnt++;
	}
	int findEd = Find(edId);
	if (findEd < 0) {
		findEd = cnt;
		nodeList[cnt].id = edId;
		findId[edId] = cnt++;
	}
	NeiNode* p = nodeList[findSt].firstNei;
	while (p) {
		if (findEd == p->id && p->valid) return;
		p = p->next;
	}

	NeiNode* e = new NeiNode;
	e->id = findEd;
	e->next = nodeList[findSt].firstNei;
	e->sup = 0;
	e->Pr = pr;
	nodeList[findSt].firstNei = e;
	nodeList[findSt].degree++;

	e = new NeiNode;
	e->id = findSt;
	e->next = nodeList[findEd].firstNei;
	e->sup = 0;
	e->Pr = pr;
	nodeList[findEd].firstNei = e;
	nodeList[findEd].degree++;
	changeEdge[make_pair(min(findSt, findEd), max(findSt, findEd))] = 0;
	visit[make_pair(min(findSt, findEd), max(findSt, findEd))] = 0;
}

void Graph::initSup() {
	//cout<<"Init support start"<<endl;
	for (int k = 0; k < nodeNum; k++) {
		NeiNode* nei = nodeList[k].firstNei;
		while (nei) {
			if (nei->valid) {
				nei->sup = computeSup(k, nei->id);
				//cout<<nei->sup<<endl;
			}
			nei = nei->next;
		}
	}
	//cout<<"Init support complish"<<endl;
}

void Graph::initSuperSup() {
	//cout<<"Init SS start"<<endl;
	for (int k = 0; k < nodeNum; k++) {
		NeiNode* nei = nodeList[k].firstNei;
		while (nei) {
			if (nei->valid) {
				nei->SS = computeSuperSup(k, nei->id, nei->tuss);
				//cout<<nei->sup<<endl;
			}
			nei = nei->next;
		}
	}
	//cout<<"Init SS complish"<<endl;
}
void Graph::initConstrainSup() {
	//cout<<"Init CS start"<<endl;
	for (int k = 0; k < nodeNum; k++) {
		NeiNode* nei = nodeList[k].firstNei;
		while (nei) {
			if (nei->valid) {
				nei->CS = computeConstrainSup(k, nei->id, nei->tuss);
				//cout<<nei->sup<<endl;
			}
			nei = nei->next;
		}
	}
	//cout<<"Init CS complish"<<endl;
}

int Graph::computeSup(int st, int ed) {
	int res = 0;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) res++;
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	return res;
}
int Graph::computeSuperSup(int st, int ed, int truss) {
	int res = 0;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id && j->tuss >= truss && i->tuss >= truss) res++;
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	return res;
}

int Graph::computeConstrainSup(int st, int ed, int truss) {
	int res = 0;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						int minn = min(i->tuss, j->tuss);
						if (minn > truss)
							res++;
						else if (minn == truss) {
							int flg = 1;
							if (i->tuss == truss) {
								if (i->SS <= truss - 2) flg = 0;
							}
							if (j->tuss == truss) {
								if (j->SS <= truss - 2) flg = 0;
							}
							if (flg == 1) res++;
						}
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	return res;
}

void Graph::cover() {
	for (int k = 0; k < nodeNum; k++) {
		NeiNode* nei = nodeList[k].firstNei;
		while (nei) {
			if (nei->valid) {
				nei->tuss = nei->sup + 2;
				if (k < nei->id)
					bin[nei->sup].push_back(make_pair(k, nei->id));
			}
			nei = nei->next;
		}
	}
}

//将边的truss设置为k
void Graph::setEdgeMess(int stId, int edId, int k, int model) {
	NeiNode* p = nodeList[stId].firstNei;
	while (p)
	{
		if (p->id == edId) {
			switch (model)
			{
			case 1:
				p->tuss = k;
				break;
			case 2:
				p->sup = k;
				break;
			case 3:
				p->SS = k;
				break;
			case 4:
				p->CS = k;
				break;
			default:
				break;
			}
			break;
		}
		p = p->next;
	}
	p = nodeList[edId].firstNei;
	while (p)
	{
		if (p->id == stId) {
			switch (model)
			{
			case 1:
				p->tuss = k;
				break;
			case 2:
				p->sup = k;
				break;
			case 3:
				p->SS = k;
				break;
			case 4:
				p->CS = k;
				break;
			default:
				break;
			}
			break;
		}
		p = p->next;
	}
	//cout<<"Set edge ("<<stId<<","<<edId<<") Truss in:"<<k<<endl;
}

int Graph::getEdgeMess(int stId, int edId, int model) {
	NeiNode* p = nodeList[stId].firstNei;
	int res;
	while (p)
	{
		if (p->id == edId) {
			switch (model)
			{
			case 1:
				res = p->tuss;
				break;
			case 2:
				res = p->sup;
				break;
			case 3:
				res = p->SS;
				break;
			case 4:
				res = p->CS;
				break;
			default:
				break;
			}
			break;
		}
		p = p->next;
	}
	return res;
}

void Graph::remEdge(int stId, int edId) {
	int st = min(stId, edId);
	int ed = max(stId, edId);
	NeiNode* i = nodeList[st].firstNei;
	while (i)
	{
		if (i->valid && i->id != ed)
		{
			NeiNode* j = nodeList[ed].firstNei;
			while (j)
			{
				if (j->valid && j->id != st)
				{
					if (j->id == i->id) {
						i->sup--;
						//  cout<<"sup of "<<nodeList[st].id<<" "<<nodeList[i->id].id<<" down"<<endl;
						NeiNode* ip = nodeList[i->id].firstNei;
						while (ip) {
							if (ip->id == st) {
								ip->sup--;
								//  cout<<"sup of "<<nodeList[i->id].id<<" "<<nodeList[ip->id].id<<" down"<<endl;
								for (int k = 0; k < bin[i->sup + 1].size(); k++) {
									if (bin[i->sup + 1][k].first == min(i->id, st) && bin[i->sup + 1][k].second == max(i->id, st))
									{
										//  cout<<"erase "<<i->sup+1<<" in "<<nodeList[bin[i->sup+1][k].first].id <<" "<<nodeList[bin[i->sup+1][k].second].id<<endl;
										bin[i->sup + 1].erase(bin[i->sup + 1].begin() + k);
										bin[i->sup].push_back(make_pair(min(i->id, st), max(i->id, st)));
										break;
									}

								}
								break;
							}
							ip = ip->next;
						}
						j->sup--;
						//cout<<"sup of "<<nodeList[ed].id<<" "<<nodeList[j->id].id<<" down"<<endl;
						ip = nodeList[j->id].firstNei;
						while (ip) {
							if (ip->id == ed) {
								ip->sup--;
								// cout<<"sup of "<<nodeList[j->id].id<<" "<<nodeList[ip->id].id<<" down"<<endl;
								for (int k = 0; k < bin[j->sup + 1].size(); k++) {
									if (bin[j->sup + 1][k].first == min(j->id, ed) && bin[j->sup + 1][k].second == max(j->id, ed))
									{
										// cout<<"erase "<<j->sup+1<<" in "<<nodeList[bin[j->sup+1][k].first].id <<" "<<nodeList[bin[j->sup+1][k].second].id<<endl;
										bin[j->sup + 1].erase(bin[j->sup + 1].begin() + k);
										bin[j->sup].push_back(make_pair(min(j->id, ed), max(j->id, ed)));
										break;
									}

								}
								break;
							}
							ip = ip->next;
						}
						break;
					}
				}
				else if (j->id == st) j->valid = 0;
				j = j->next;
			}
		}
		else if (i->id == ed) i->valid = 0;
		i = i->next;
	}
	//edgeNum--;
}

void Graph::greed() {
	int k = 2;
	while (1)
	{
		int minSup = 0;
		while (bin[minSup].size() == 0) {
			minSup++;
			if (minSup > nodeNum - 2) return;
		}

		//cout<<"minSup:"<<minSup<<endl;
		while (minSup > k - 2) k++;
		int stId = bin[minSup][0].first;
		int edId = bin[minSup][0].second;
		bin[minSup].erase(bin[minSup].begin() + 0);
		//cout<<nodeList[stId].id<<" "<<nodeList[edId].id<<endl;
		setEdgeMess(stId, edId, k, 1);
		remEdge(stId, edId);
		//output();
	}
}

void Graph::distribute() {
	cout << "Start DIS_Algrithm" << endl;
	int maxIter = 99999999;
	int cnt_iter = 0;
	while (maxIter--) {
		int changeEdgeNum = 0;
		int changeNodeNum = 0;
		for (int k = 0; k < nodeNum; k++) {
			NeiNode* nei = nodeList[k].firstNei;
			int oldChangeEdgeNum = changeEdgeNum;
			//  cout<<"Node:"<<k<<endl;
			while (nei) {
				if (nei->valid) {
					vector<int> newTuss;
					int maxTuss = -1;
					NeiNode* i = nodeList[k].firstNei;
					while (i)
					{
						if (i->valid && i->id != nei->id) {
							NeiNode* j = nodeList[nei->id].firstNei;
							while (j)
							{
								if (j->valid && j->id != k) {
									if (j->id == i->id) {
										newTuss.push_back(min(i->tuss, j->tuss));
										maxTuss = max(maxTuss, min(i->tuss, j->tuss));
									}
								}
								j = j->next;
							}

						}
						i = i->next;
					}
					// cout<<"nei(k-2 >=k):"<<nei->id<<endl;
					while (maxTuss > 1) {
						int cnt = 0;
						for (int q = 0; q < newTuss.size(); q++) {
							//       cout<<newTuss[q]<<" ";
							if (newTuss[q] >= maxTuss)
								cnt++;
						}
						if (cnt >= maxTuss - 2) break;
						maxTuss--;
					}
					//   cout<<endl;
					if (maxTuss > 1 && maxTuss < nei->tuss) {
						nei->tuss = maxTuss;
						changeEdgeNum++;
					}
				}
				nei = nei->next;
			}
			if (changeEdgeNum - oldChangeEdgeNum > 0)changeNodeNum++;
		}
		if (startCntSteps) {
			totalMsgNumber += changeNodeNum;
			totalSteps += 1;
		}
		//cout<<"第"<<++cnt_iter<<"轮： 更新点数-"<<changeNodeNum<<"   更新边数-"<<changeEdgeNum<<endl;
		cout << "Step " << setw(4) << ++cnt_iter << ":  Changed Point Num = " << setw(4) << changeNodeNum << "  Changed Edge(Message) Num:" << setw(4) << changeEdgeNum << endl;
		if (changeEdgeNum == 0) {
			break;
		}
	}
	cout << endl;
}

void Graph::output() {
	cout << "Node:" << nodeNum << " Edge:" << edgeNum << endl;
	for (int i = 0; i < nodeNum; i++)
	{
		cout << "Point:" << nodeList[i].id << " Degree:" << nodeList[i].degree << "  Arc:" << endl;
		NeiNode* nei = nodeList[i].firstNei;
		while (nei)
		{
			cout << "Id:" << nodeList[nei->id].id << "  Valid:" << nei->valid << " Sup:" << nei->sup << " Tuss:" << nei->tuss << endl;
			nei = nei->next;
		}
	}

	cout << "Bin:" << endl;
	for (int i = 0; i < nodeNum - 1; i++)
	{
		cout << "i:" << i << "    >>>>>>  ";
		for (int j = 0; j < bin[i].size(); j++)
		{
			cout << "<" << nodeList[bin[i][j].first].id << "," << nodeList[bin[i][j].second].id << ">" << " ";
		}
		cout << endl;
	}


}

void Graph::writeFile(char* filename) {
	FILE* fp;//定义一个文件指针
	string tmp = "";
	for (int i = 0; i < strlen(filename); i++)
		tmp += filename[i];
	cout << filename << endl;
	fp = fopen(strcat(filename, ".truss"), "w");//以只写的方式打开文件，前面的参数是文件路径，后面的参数是表示只写
	fprintf(fp, "# Nodes: %d Edges: %d\n", nodeNum, edgeNum);
	fprintf(fp, "Node 1\tNode 2\tTruss\n");
	int totalNnum = 0;
	map<int, int > InPercent;
	for (int i = 0; i < nodeNum; i++)
	{
		NeiNode* nei = nodeList[i].firstNei;
		while (nei)
		{
			if (nodeList[i].id < nodeList[nei->id].id) {
				fprintf(fp, "%d\t%d\t%d\n", nodeList[i].id, nodeList[nei->id].id, nei->tuss);
				totalNnum++;
				InPercent[nei->tuss]++;
			}
			nei = nei->next;
		}
	}
	filename = (char*)tmp.c_str();
	fp = fopen(strcat(filename, ".percent"), "w");//以只写的方式打开文件，前面的参数是文件路径，后面的参数是表示只写
	fprintf(fp, "Total Affected nodes number:%d\n", changeNodes.size());
	fprintf(fp, "The depth of broadcast:%d\n", maxDepth);
	fprintf(fp, "Total Message Number:%d\n", totalMsgNumber);
	fprintf(fp, "Total steps:%d\n\n", totalSteps);
	fprintf(fp, "Truss\tPercent\n");
	map<int, int >::iterator it;
	for (it = InPercent.begin(); it != InPercent.end(); it++)
	{
		fprintf(fp, "%d\t%f\n", it->first, (double)(((double)it->second) / totalNnum));
	}
}


//初始化Pr SuP  每条边的sup值用dp计算  并且将sup覆盖truss
void Graph::initPrSup() {
	for (int k = 0; k < nodeNum; k++) {
		NeiNode* nei = nodeList[k].firstNei;
		while (nei) {
			if (nei->valid) {
				nei->PrSup = computePrSup(k, nei->id, nei->Pr);
				nei->tuss = nei->PrSup + 2;
				if (k < nei->id) {
					edge e(k, nei->id, nei->PrSup);
					edgeSet.insert(e);
					//cout<<k<<" "<<nei->id<<" "<<nei->PrSup<<endl;
				}
				//nei->tuss = nei->PrSup+2;
			}
			nei = nei->next;
		}
	}
}

//dp及计算sup dp[i][j] = dp[i-1][j-1]p(st,j)p(ed,j) +  (1-p(st,j)p(ed,j))dp[i-1][j-1]
//前i个三角形里面有j个是存在的
double Graph::computePrSup(int st, int ed, double p) {
	NeiNode* i = nodeList[st].firstNei;

	// cout<<"<"<<st<<","<<ed<<">:";
	int ek = 0;
	vector<double> pr;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						pr.push_back(i->Pr * j->Pr);
						ek++;
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	// cout<<" = "<<PrSup<<endl;
	double* e = new double[ek + 1];
	double** f = new double* [ek + 1];
	for (int i = 0; i < ek + 1; i++)
	{
		f[i] = new double[ek + 1];
	}

	//g++不允许此方法计算
	//double e[ek + 1];
	//double f[ek + 1][ek + 1];
	f[0][0] = 1;
	// for(int k =1;k<=ek+1;k++) f[0][k] = 0;
	for (int l = 0; l <= ek; l++) {
		for (int k = l + 1; k <= ek; k++) {
			f[k][l] = 0;
		}
	}
	for (int l = 1; l <= ek; l++) {
		f[0][l] = (1 - pr[l - 1]) * f[0][l - 1];
	}
	for (int k = 1; k <= ek; k++) {
		for (int l = 1; l <= ek; l++) {
			f[k][l] = pr[l - 1] * f[k - 1][l - 1] + (1 - pr[l - 1]) * f[k][l - 1];
		}
	}
	e[0] = 1;
	for (int t = 1; t <= ek; t++) {
		e[t] = e[t - 1] - f[t][ek];//\sum pr(sup>=t) sup大于等于t的所有情况的概率 
		//cout<<"t:"<<t<<" e[t]:"<<e[t]<<" e[t-1]:"<<e[t-1]<<" p:"<<p<<" e[t]*p:"<<e[t]*p<<" f[t][ek]:"<<f[t][ek]<<" gamma:"<<gamma<<endl;
		if (e[t] * p < gamma) return t - 1;
	}

	return ek;
}

void Graph::PrGreed() {
	int k = 2;
	while (1)
	{
		if (edgeSet.size() == 0) break;
		int minSup = (*edgeSet.begin()).sup;

		//cout<<"minSup:"<<minSup<<endl;
		while (minSup > k - 2) k++;
		int stId = (*edgeSet.begin()).st;
		int edId = (*edgeSet.begin()).ed;
		edgeSet.erase(*edgeSet.begin());
		//cout<<nodeList[stId].id<<" - "<<nodeList[edId].id<<endl;
		setEdgeMess(stId, edId, k, 1);
		PrRemEdge(stId, edId);
	}
}

void Graph::PrRemEdge(int stId, int edId) {
	int st = min(stId, edId);
	int ed = max(stId, edId);
	deleteEdge(st, ed);
	NeiNode* i = nodeList[st].firstNei;
	while (i)
	{
		if (i->valid && i->id != ed)
		{
			NeiNode* j = nodeList[ed].firstNei;
			while (j)
			{
				if (j->valid && j->id != st)
				{
					if (j->id == i->id) {
						edge temp1(min(i->id, st), max(i->id, st), i->PrSup);
						edgeSet.erase(temp1);
						i->PrSup = computePrSup(st, i->id, i->Pr);
						edge temp2(min(i->id, st), max(i->id, st), i->PrSup);
						edgeSet.insert(temp2);

						edge temp3(min(j->id, ed), max(j->id, ed), j->PrSup);
						edgeSet.erase(temp3);
						j->PrSup = computePrSup(ed, j->id, j->Pr);
						edge temp4(min(j->id, ed), max(j->id, ed), j->PrSup);
						edgeSet.insert(temp4);
						//cout<<"sup of "<<nodeList[st].id<<" "<<nodeList[i->id].id<<" down to (i,j):"<<i->PrSup<<" "<<j->PrSup<<endl;
						break;
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}

}

void Graph::deleteEdge(int st, int ed) {
	//cout<<"st:"<<st<<" ed:"<<ed<<endl;
	NeiNode* i = nodeList[st].firstNei;
	while (i)
	{
		if (i->valid && i->id == ed)
		{
			i->valid = 0;
		}
		i = i->next;
	}
	NeiNode* j = nodeList[ed].firstNei;
	while (j)
	{
		if (j->valid && j->id == st)
		{
			j->valid = 0;
		}
		j = j->next;
	}
}

void Graph::outputSet() {
	set<edge>::iterator iter = edgeSet.begin();
	cout << "Set size:" << edgeSet.size() << endl;
	while (iter != edgeSet.end())
	{
		cout << (*iter).st << " " << (*iter).ed << endl;
		iter++;
	}
}

void Graph::PrDistribute() {
	int maxIter = 99999999;
	int cnt_iter = 0;
	while (maxIter--) {
		int changeEdgeNum = 0;
		int changeNodeNum = 0;
		for (int k = 0; k < nodeNum; k++) {
			NeiNode* nei = nodeList[k].firstNei;
			//  cout<<"Node:"<<k<<endl;
			int oldChangeEdgeNum = changeEdgeNum;
			while (nei) {
				if (nei->valid) {
					vector<int> newTuss;
					int maxTuss = -1;
					NeiNode* i = nodeList[k].firstNei;
					while (i)
					{
						if (i->valid && i->id != nei->id) {
							NeiNode* j = nodeList[nei->id].firstNei;
							while (j)
							{
								if (j->valid && j->id != k) {
									if (j->id == i->id) { //形成了一个三角形
										newTuss.push_back(min(i->tuss, j->tuss));//放入三角形临边较小的truss
										maxTuss = max(maxTuss, min(i->tuss, j->tuss));//记录数组中最大的maxTruss
									}
								}
								j = j->next;
							}

						}
						i = i->next;
					}
					// cout<<"nei(k-2 >=k):"<<nei->id<<endl;
					while (maxTuss > 1) {
						int cnt = 0;
						for (int q = 0; q < newTuss.size(); q++) {
							//       cout<<newTuss[q]<<" ";
							if (newTuss[q] >= maxTuss)
								cnt++;
						}
						if (cnt >= maxTuss - 2) break;
						maxTuss--;
					}
					//   cout<<endl;
					if (maxTuss > 1 && maxTuss < nei->tuss) {
						nei->tuss = maxTuss;
						changeEdgeNum++;
					}
				}
				nei = nei->next;
			}
			if (changeEdgeNum - oldChangeEdgeNum > 0)changeNodeNum++;
		}
		//cout<<"第"<<++cnt_iter<<"轮： 更新点数-"<<changeNodeNum<<"   更新边数-"<<changeEdgeNum<<endl;
		cout << "Step " << setw(4) << ++cnt_iter << ":  Changed Point Num = " << setw(4) << changeNodeNum << "  Changed Edge Num:" << setw(4) << changeEdgeNum << endl;
		if (changeEdgeNum == 0) {
			break;
		}
	}

}

bool Graph::inCntEdge(int st, int ed) {
	int stId = min(st, ed);
	int edId = max(st, ed);
	int it = visit[make_pair(min(st, ed), max(st, ed))];;
	if (it == 0) return false;
	return true;
}

void Graph::upAdjust_2(int st, int ed, int deep) {
	NeiNode* i = nodeList[st].firstNei;
	int t = i->tuss;
	visit[make_pair(min(st, ed), max(st, ed))] = 1;
	maxDepth = max(maxDepth, deep);
	//cout<<"fff"<< cntEdge.size()<<endl;
	//cout<<"ADD@2:"<<min(st,ed)<<" "<<max(st,ed)<<"  "<<endl;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						if (i->tuss == t && j->tuss >= t && !inCntEdge(st, i->id)) {
							//setEdgeMess(st,i->id,1+getTuss(st,i->id));
							//changeEdge[make_pair(min(st,i->id),max(st,i->id))] += 1 ;
							//cout<<"ADD@2:"<<min(st,i->id)<<" "<<max(st,i->id)<<"  "<<endl;
							cntEdge.insert(make_pair(min(st, i->id), max(st, i->id)));
							upAdjust_2(min(st, i->id), max(st, i->id), deep + 1);
							changeNodes.insert(st);
							changeNodes.insert(i->id);
						}
						if (j->tuss == t && i->tuss >= t && !inCntEdge(ed, j->id)) {
							//setEdgeMess(ed,j->id,1+getTuss(ed,j->id));
							//changeEdge[make_pair(min(ed,i->id),max(ed,i->id))] += 1 ;
							//cout<<"ADD@2:"<<min(ed,i->id)<<" "<<max(ed,i->id)<<"  "<<endl;
							cntEdge.insert(make_pair(min(ed, i->id), max(ed, i->id)));
							upAdjust_2(min(ed, i->id), max(ed, i->id), deep + 1);
							changeNodes.insert(ed);
							changeNodes.insert(j->id);
						}
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	//cout<<"add2 return"<<endl;
}

void Graph::upAdjust(int st, int ed) {

	NeiNode* i = nodeList[st].firstNei;
	int t = i->tuss;
	visit[make_pair(min(st, ed), max(st, ed))] = 1;;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						//cout<<i->tuss<<" "<<j->tuss<<" "<<t<<endl;
						if (i->tuss < t && !inCntEdge(st, i->id)) {
							//setEdgeMess(st,i->id,1+getTuss(st,i->id));
							//changeEdge[make_pair(min(st,i->id),max(st,i->id))] += 1;
							cntEdge.insert(make_pair(min(st, i->id), max(st, i->id)));
							//cout<<"ADD:"<<min(st,i->id)<<" "<<max(st,i->id)<<"  "<<changeEdge[make_pair(min(st,i->id),max(st,i->id))]<<endl;
							upAdjust_2(min(st, i->id), max(st, i->id), 1);
							changeNodes.insert(t);
							changeNodes.insert(i->id);
						}
						if (j->tuss < t && !inCntEdge(ed, j->id)) {
							//setEdgeMess(ed,j->id,1+getTuss(ed,j->id));
							//changeEdge[make_pair(min(ed,i->id),max(ed,i->id))] += 1 ;
							cntEdge.insert(make_pair(min(ed, i->id), max(ed, i->id)));
							//cout<<"ADD:"<<min(ed,i->id)<<" "<<max(ed,i->id)<<"  "<<changeEdge[make_pair(min(ed,i->id),max(ed,i->id))]<<endl;
							upAdjust_2(min(ed, i->id), max(ed, i->id), 1);
							changeNodes.insert(t);
							changeNodes.insert(j->id);
						}
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	//cout<<"add return"<<endl;

}

void Graph::dynamicInsert(int stId, int edId) {
	//cout<<"Insert:"<<stId<<" "<<edId<<"<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	int s = stId;
	int e = edId;
	addEdge(stId, edId, 0);
	stId = Find(stId);
	edId = Find(edId);
	int sp = computeSup(stId, edId);
	nodeList[stId].firstNei->sup = sp;
	nodeList[edId].firstNei->sup = sp;
	setEdgeMess(stId, edId, sp + 2, 1);
	boundData bd = bound(stId, edId);
	//cout<<"Dot:"<<s<<" "<<e<<" ->bound:"<<bd.up_bound<<" sup+2:" <<sp+2<<endl;

	setEdgeMess(stId, edId, min(sp + 2, bd.up_bound), 1);
	//cout<<nodeList[stId].id <<" "<<nodeList[edId].id<<" Up_bound:"<<min(sp+2,bd.up_bound)<<endl;
	upAdjust(stId, edId);
	//cout<<"Adjust finished"<<endl;
	/*
	set<pair<int,int> >::iterator it;
	for(it =cntEdge.begin();it!=cntEdge.end();it++){
		visit[*it] = 0;
	}
	cntEdge.clear();
	//cout<<"Insert finish:"<<s<<" "<<e<<"<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	*/
	set<pair<int, int> >::iterator it;
	for (it = cntEdge.begin(); it != cntEdge.end(); it++) {
		//cout<<(*it).first<<" "<<(*it).second<<"  "<<changeEdge[*it]<<" "<<getTuss((*it).first,(*it).second)<<" "<<getSup((*it).first,(*it).second)+2<<endl;
		int a = computeSup((*it).first, (*it).second) + 2;
		int b = getEdgeMess((*it).first, (*it).second, 1) + 1;//changeEdge[*it];
		//cout<<"xutao xitong d :"<<changeEdge[*it]<<endl;
		int newTuss = min(a, b);
		//cout<<a<<" "<<b<<endl;
//cout<<nodeList[(*it).first].id<<" "<<nodeList[(*it).second].id<<" newTruss:"<< newTuss<<"  getTruss and sup:"<<getTuss((*it).first,(*it).second)<<endl;
		setEdgeMess((*it).first, (*it).second, newTuss, 1);
		//changeEdge[*it] = 0;
		visit[*it] = 0;
	}
	cntEdge.clear();
	visit[make_pair(min(stId, edId), max(stId, edId))] = 0;
}

void Graph::downAdjust(int st, int ed) {
	NeiNode* i = nodeList[st].firstNei;
	int t = i->tuss;
	visit[make_pair(min(st, ed), max(st, ed))] = 1;;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						if (i->tuss <= t) {
							cntEdge.insert(make_pair(min(st, i->id), max(st, i->id)));
							//cout<<"ADD:"<<min(st,i->id)<<" "<<max(st,i->id)<<"  "<<changeEdge[make_pair(min(st,i->id),max(st,i->id))]<<endl;
							upAdjust_2(min(st, i->id), max(st, i->id), 1);
							changeNodes.insert(t);
							changeNodes.insert(i->id);
						}
						if (j->tuss <= t) {
							cntEdge.insert(make_pair(min(ed, i->id), max(ed, i->id)));
							//cout<<"ADD:"<<min(ed,i->id)<<" "<<max(ed,i->id)<<"  "<<changeEdge[make_pair(min(ed,i->id),max(ed,i->id))]<<endl;
							upAdjust_2(min(ed, i->id), max(ed, i->id), 1);
							changeNodes.insert(t);
							changeNodes.insert(j->id);
						}
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
}

void Graph::dynamicDelete(int stId, int edId) {
	//cout<<"Start delete:"<<stId<<" "<<edId<<endl;
	//addEdge(stId,edId);
	stId = Find(stId);
	edId = Find(edId);
	//boundData bd = bound(stId,edId);
	remEdge(stId, edId);
	downAdjust(stId, edId);
	//cout<<"Adjust finished"<<endl;
	set<pair<int, int> >::iterator it;
	for (it = cntEdge.begin(); it != cntEdge.end(); it++) {
		int a = computeSup((*it).first, (*it).second) + 2;
		int b = getEdgeMess((*it).first, (*it).second, 1);
		int newTuss = min(a, b);
		//cout<<(*it).first<<" "<<(*it).second<<" "<< newTuss<<" "<<a<<" "<<b<<endl;
		setEdgeMess((*it).first, (*it).second, newTuss, 1);
		visit[*it] = 0;
	}
	cntEdge.clear();
	//cout<<"Delete finished"<<endl;
	visit[make_pair(min(stId, edId), max(stId, edId))] = 0;
}

boundData Graph::bound(int st, int ed) {
	vector<int > v;
	boundData res;
	res.up_bound = 9999999;
	res.low_bound = 9999999;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						int a = i->tuss;// getTuss(i->id,st);
						int b = j->tuss;//getTuss(j->id,ed);
						int d = min(a, b);
						v.push_back(d);
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	int up = 0;
	int low = 0;
	int maxTuss = v.size() + 2;
	while (maxTuss > 2) {
		int cnt = 0;
		int cnt2 = 0;
		for (int q = 0; q < v.size(); q++) {
			//       cout<<newTuss[q]<<" ";
			if (v[q] >= maxTuss - 1)
				cnt++;
			if (v[q] >= maxTuss)
				cnt2++;
		}
		if (cnt >= maxTuss - 2 && up == 0) {
			res.up_bound = maxTuss;
			up = 1;
		}
		if (cnt2 >= maxTuss - 2 && low == 0) {
			res.low_bound = maxTuss;
			low = 1;
		}
		if (up && low) break;
		maxTuss--;
	}
	return res;

}

void Graph::supInitDelete(int stId, int edId) {
	stId = Find(stId);
	edId = Find(edId);
	int st = min(stId, edId);
	int ed = max(stId, edId);
	deleteEdge(st, ed);
}

int Graph::computeLB(int st, int ed) {
	vector<int > v;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						int a = i->tuss;// getTuss(i->id,st);
						int b = j->tuss;//getTuss(j->id,ed);
						int d = min(a, b);
						v.push_back(d);
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	//for(int i=0;i<v.size();i++) cout<<v[i]<<" ";
	//cout<<endl;
	int maxTuss = v.size() + 2;
	int LB = maxTuss;
	while (maxTuss > 2) {
		int cnt = 0;
		for (int q = 0; q < v.size(); q++) {
			//       cout<<newTuss[q]<<" ";
			if (v[q] >= maxTuss)
				cnt++;
		}
		if (cnt >= maxTuss - 2) {
			LB = maxTuss;
			break;
		}
		maxTuss--;
	}
	//cout<<LB<<endl;
	return LB;
}

int Graph::computeUB(int st, int ed) {
	vector<int > v;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						int a = i->tuss;// getTuss(i->id,st);
						int b = j->tuss;//getTuss(j->id,ed);
						int d = min(a, b);
						v.push_back(d);
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
	//for(int i=0;i<v.size();i++) cout<<v[i]<<" ";
	//cout<<endl;
	int maxTuss = v.size() + 2;
	int UB = maxTuss;
	while (maxTuss > 2) {
		int cnt = 0;
		for (int q = 0; q < v.size(); q++) {
			//       cout<<newTuss[q]<<" ";
			if (v[q] >= maxTuss - 1)
				cnt++;
		}
		if (cnt >= maxTuss - 2) {
			UB = maxTuss;
			break;
		}
		maxTuss--;
	}
	//cout<<UB<<endl;
	return UB;
}

void Graph::centerAdjust(set<edge, cmp> PES) {
	//set<pair<int,int> > resSet;

	set<edge >::iterator iter;
	/*
	cout<<"PES:"<<endl;
	for(iter = PES.begin();iter!=PES.end();iter++){
		cout<<(*iter).st<<" "<<(*iter).ed<<endl;
	}
	cout<<"PES_end:"<<endl;
	*/
	for (iter = PES.begin(); iter != PES.end(); iter++) {
		while (!Stack.empty()) Stack.pop();
		pair<int, int> e_root = make_pair((*iter).st, (*iter).ed);
		int rootTruss = getEdgeMess(e_root.first, e_root.second, 1);
		Sset[e_root] = getEdgeMess(e_root.first, e_root.second, 4);
		Stack.push(e_root);
		Vset[e_root] = 1;
		while (!(Stack.empty())) {
			pair<int, int> now_e = Stack.top();
			Stack.pop();
			//cout<<Sset[now_e]<<" "<<rootTruss-2<<endl;
			if (Sset[now_e] > rootTruss - 2) {
				//cout<<"Change:"<<now_e.first<<" "<<now_e.second<<endl;
				//resSet.insert(now_e);
				int a = now_e.first;
				int b = now_e.second;
				NeiNode* i = nodeList[a].firstNei;
				while (i) {
					if (i->valid && i->id != b) {
						NeiNode* j = nodeList[b].firstNei;
						while (j) {
							if (j->valid && j->id != a) {
								if (j->id == i->id) {
									//cout<<"Turss 1("<<a<<","<<i->id<<"):"<<i->tuss<<" "<<"Truss 2("<<b<<","<<i->id<<"):"<<j->tuss<<"    rootTruss:"<<rootTruss<<endl;
									if (min(i->tuss, j->tuss) >= rootTruss) {
										pair<int, int> new_e_1 = make_pair(min(i->id, a), max(i->id, a));
										if (i->tuss == rootTruss && i->SS > rootTruss - 2 && Vset[new_e_1] == 0 && j->tuss > rootTruss) {
											Stack.push(new_e_1);
											//cout<<"Add:"<<nodeList[min(i->id,a)].id<<" "<<nodeList[max(i->id,a)].id<<endl;
											Vset[new_e_1] = 1;
											Sset[new_e_1] = Sset[new_e_1] + getEdgeMess(new_e_1.first, new_e_1.second, 4);
										}
										pair<int, int> new_e_2 = make_pair(min(i->id, b), max(i->id, b));
										if (j->tuss == rootTruss && j->SS > rootTruss - 2 && Vset[new_e_2] == 0 && i->tuss > rootTruss) {
											Stack.push(new_e_2);
											//cout<<"Add:"<<nodeList[min(i->id,b)].id<<" "<<nodeList[max(i->id,b)].id<<endl;
											Vset[new_e_2] = 1;
											Sset[new_e_2] = Sset[new_e_2] + getEdgeMess(new_e_2.first, new_e_2.second, 4);
										}
										if (j->tuss == rootTruss && i->tuss == rootTruss && i->SS > rootTruss - 2 && j->SS > rootTruss - 2) {
											if (Vset[new_e_2] == 0) {
												Stack.push(new_e_2);
												//cout<<"Add:"<<nodeList[min(i->id,b)].id<<" "<<nodeList[max(i->id,b)].id<<endl;
												Vset[new_e_2] = 1;
												Sset[new_e_2] = Sset[new_e_2] + getEdgeMess(new_e_2.first, new_e_2.second, 4);
											}
											if (Vset[new_e_1] == 0) {
												Stack.push(new_e_1);
												//cout<<"Add:"<<nodeList[min(i->id,a)].id<<" "<<nodeList[max(i->id,a)].id<<endl;
												Vset[new_e_1] = 1;
												Sset[new_e_1] = Sset[new_e_1] + getEdgeMess(new_e_1.first, new_e_1.second, 4);
											}

										}
									}
								}
							}
							j = j->next;
						}
					}
					i = i->next;
				}
				//cout<<"--------------------------------------------"<<endl<<endl;
			}
			else {
				//cout<<"Qu zhu:"<<now_e.first<<" "<<now_e.second<<endl;
				if (Xset[now_e] == 0) {
					Xset[now_e] = 1;
					Eliminate(now_e, rootTruss);
				}
			}

		}
	}
	/*
	set<pair<int,int> >::iterator it;
	for(it = resSet.begin();it!=resSet.end();it++){
		cout<<nodeList[(*it).first].id<<" "<<nodeList[(*it).second].id<<" "<<getEdgeMess( (*it).first,(*it).second,1)<<endl;
		setEdgeMess( (*it).first,(*it).second,getEdgeMess( (*it).first,(*it).second,1)+1,1);
	}
	*/

	//cout<<"UP edge: "<<endl;
	map<pair<int, int>, int >::iterator vit;
	for (vit = Vset.begin(); vit != Vset.end(); vit++) {
		//cout<<nodeList[(*vit).first.first].id<<" "<<nodeList[(*vit).first.second].id<<" "<<(*vit).second<<" "<<Xset[(*vit).first]<<endl;
		if (Vset[(*vit).first] == 1 && Xset[(*vit).first] == 0) {
			//cout<<nodeList[(*vit).first.first].id<<" "<<nodeList[(*vit).first.second].id<<endl;
			setEdgeMess((*vit).first.first, (*vit).first.second, getEdgeMess((*vit).first.first, (*vit).first.second, 1) + 1, 1);
		}
	}

}

bool Graph::exist_edge(int st, int ed) {
	if (Find(st) < 0 || Find(ed) < 0) return false;
	int findSt = Find(st);
	int findEd = Find(ed);
	NeiNode* p = nodeList[findSt].firstNei;
	while (p) {
		if (findEd == p->id && p->valid) return true;
		p = p->next;
	}
	return false;
}




void Graph::centerInsert(int stId, int edId, int initModle) {
	//cout<<"Insert>>>>>>"<<endl;
	//cout<<stId<<" "<<edId<<" <<<<<  input "<<endl;
	if (exist_edge(stId, edId)) {
		cout << "Repeated edge." << endl;
		return;
	}

	addEdge(stId, edId);
	stId = Find(stId);
	edId = Find(edId);
	//cout<<stId<<" "<<edId<<" "<<nodeNum<<endl;
	int st = min(stId, edId);
	int ed = max(stId, edId);

	double tmpTime = 0;
	startTime = clock();
	int LB = computeLB(st, ed);
	//int UB = computeUB(st,ed);
	int UB = 0;
	int myTruss = LB;

	setEdgeMess(st, ed, LB, 1);//Truss
	int mySS, myCS;
	endTime = clock();
	tmpTime += 1000.0 * (endTime - startTime);
	startTime = clock();
	initSuperSup();
	initConstrainSup();
	mySS = getEdgeMess(st, ed, 3);
	myCS = getEdgeMess(st, ed, 4);

	//cout<<">>"<<mySS<<" "<<myCS<<endl;
	if (myCS > myTruss - 2) UB = LB + 1;
	else UB = LB;
	//cout<<LB<<" "<<UB<<endl;
	//setEdgeMess(st,ed,UB,1);
	//UB = computeUB(st,ed);

	startTime = clock();

	set<edge, cmp> PES;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						if (i->tuss <= UB) PES.insert(edge(min(i->id, st), max(i->id, st), i->tuss));
						if (j->tuss <= UB) PES.insert(edge(min(j->id, ed), max(j->id, ed), j->tuss));
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}

	Vset.clear();
	Xset.clear();
	Sset.clear();

	centerAdjust(PES);
	endTime = clock();
	tmpTime += 1000.0 * (endTime - startTime);
	totalTime += tmpTime;
	//log(1,1,time);

	//cout<<"Time of insert One edge:"<<time<<endl;
}

void Graph::Eliminate(pair<int, int> e, int t) {
	Xset[e] = 1;
	int k = t - 2;
	int a = e.first;
	int b = e.second;
	NeiNode* i = nodeList[a].firstNei;
	while (i) {
		if (i->valid && i->id != b) {
			NeiNode* j = nodeList[b].firstNei;
			while (j) {
				if (j->valid && j->id != a) {
					if (j->id == i->id) {

						if (min(i->tuss, j->tuss) > t) {
							/*
							pair<int,int> new_e_1 = make_pair(min(i->id,a),max(i->id,a));
							if(i->tuss == t ){
								Sset[new_e_1] -- ;
								if(Sset[new_e_1] == k && Xset[new_e_1]==0)
									Eliminate(new_e_1,t);
							}
							pair<int,int> new_e_2 = make_pair(min(j->id,b),max(j->id,b));
							if(j->tuss == t ){
								Sset[new_e_2] -- ;
								if(Sset[new_e_2] == k && Xset[new_e_2]==0)
									Eliminate(new_e_2,t);
							}
							*/
						}
						else if (j->tuss == t && i->tuss == j->tuss) {
							pair<int, int> new_e_1 = make_pair(min(i->id, a), max(i->id, a));
							pair<int, int> new_e_2 = make_pair(min(j->id, b), max(j->id, b));
							if (i->SS > k && j->SS > k) {
								Sset[new_e_1] --; Sset[new_e_2] --;
								if (Xset[new_e_1] == 0 && Sset[new_e_1] == k) {
									Eliminate(new_e_1, t);
								}
								if (Xset[new_e_2] == 0 && Sset[new_e_2] == k) {
									Eliminate(new_e_2, t);
								}
							}
						}
						else if (i->tuss > t && j->tuss == t) {
							pair<int, int> new_e_2 = make_pair(min(j->id, b), max(j->id, b));
							Sset[new_e_2] --;
							if (Xset[new_e_2] == 0 && Sset[new_e_2] == k) {
								Eliminate(new_e_2, t);
							}
						}
						else if (j->tuss > t && i->tuss == t) {
							pair<int, int> new_e_1 = make_pair(min(j->id, b), max(j->id, b));
							Sset[new_e_1] --;
							if (Xset[new_e_1] == 0 && Sset[new_e_1] == k) {
								Eliminate(new_e_1, t);
							}
						}

					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}
}

void Graph::enableAllNodes() {
	for (int k = 0; k < nodeNum; k++) {
		NeiNode* nei = nodeList[k].firstNei;
		while (nei) {
			nei->valid = 1;
			nei = nei->next;
		}
	}
}
bool GreaterSort(edge a, edge b) {
	if (a.st != b.st)
		return a.st < b.st;
	return a.ed < b.ed;
}
bool LessSort(edge a, edge b) {
	if (a.st != b.st)
		return a.st < b.st;
	return a.ed < b.ed;
}


bool Graph::noCross(vector<edge> eset, int st, int ed) {
	st = Find(st);
	ed = Find(ed);
	for (int k = 0; k < eset.size(); k++) {
		map<pair<int, int>, int> a;
		NeiNode* i = nodeList[st].firstNei;
		while (i) {
			if (i->valid && i->id != ed) {
				a[make_pair(min(i->id, st), max(i->id, st))] = 1;
			}
			i = i->next;
		}
		i = nodeList[ed].firstNei;
		while (i) {
			if (i->valid && i->id != st) {
				a[make_pair(min(i->id, ed), max(i->id, ed))] = 1;
			}
			i = i->next;
		}
		int st2 = Find(eset[k].st);
		int ed2 = Find(eset[k].ed);
		map<pair<int, int>, int> b;
		i = nodeList[st2].firstNei;
		while (i) {
			if (i->valid && i->id != ed2) {
				if (a[make_pair(min(i->id, st2), max(i->id, st2))] == 1) return 0;
			}
			i = i->next;
		}
		i = nodeList[ed2].firstNei;
		while (i) {
			if (i->valid && i->id != st2) {
				if (a[make_pair(min(i->id, ed2), max(i->id, ed2))] == 1) return 0;
			}
			i = i->next;
		}
	}
	return 1;
}

void Graph::centerMultInsert(vector<int> stIds, vector<int > edIds) {
	//initSuperSup();
	//initConstrainSup();
	double maxTime;
	vector<edge> insertEdge;
	for (int i = 0; i < stIds.size(); i++) {
		MapToIndex(stIds[i]);
		MapToIndex(edIds[i]);
		edge tmp(stIds[i], edIds[i], 2);
		insertEdge.push_back(tmp);
	}
	int cntStep = 0;
	int maxNumberOneCycle = 0;
	while (!insertEdge.empty()) {
		initSuperSup();
		initConstrainSup();

		double tmpTime = 0;
		startTime = clock();
		cout << "Step " << cntStep++ << "->" << endl;
		// set<edge,cmp>::iterator iter;
		int cntNuumberOnecycle = 0;
		vector<edge > nowInsertEdgeSet;
		addEdge(insertEdge[0].st, insertEdge[0].ed);
		int a = Find(insertEdge[0].st);
		int b = Find(insertEdge[0].ed);
		int st = min(a, b);
		int ed = max(a, b);
		int LB = computeLB(st, ed);
		setEdgeMess(st, ed, LB, 1);
		insertEdge[0].st2 = LB;
		nowInsertEdgeSet.push_back(insertEdge[0]);

		insertEdge.erase(insertEdge.begin() + 0);
		for (int k = 0; k < insertEdge.size(); k++) {
			if (noCross(nowInsertEdgeSet, insertEdge[k].st, insertEdge[k].ed)) {
				addEdge(insertEdge[k].st, insertEdge[k].ed);
				int a = Find(insertEdge[k].st);
				int b = Find(insertEdge[k].ed);
				int st = min(a, b);
				int ed = max(a, b);
				int LB = computeLB(st, ed);
				setEdgeMess(st, ed, LB, 1);
				insertEdge[k].st2 = LB;

				nowInsertEdgeSet.push_back(insertEdge[k]);
				//cout<<"fuck:"<<k<<" "<<insertEdge.size()<<" "<<insertEdge[k].st<<" "<<insertEdge[k].ed<<" "<<a<<" "<<b<<endl; 
				insertEdge.erase(insertEdge.begin() + k);
				k--;
				cntNuumberOnecycle++;

			}
			else continue;
		}
		endTime = clock();
		tmpTime += 1000.0 * (endTime - startTime);
		initSuperSup();
		initConstrainSup();
		startTime = clock();
		set<edge, cmp> PES;
		for (int k = 0; k < nowInsertEdgeSet.size(); k++) {
			int a = Find(nowInsertEdgeSet[k].st);
			int b = Find(nowInsertEdgeSet[k].ed);
			//cout<<nowInsertEdgeSet[k].st<<" "<<nowInsertEdgeSet[k].ed<<" "<<a<<" "<<b<<endl;
			int st = min(a, b);
			int ed = max(a, b);
			int LB = nowInsertEdgeSet[k].st2;
			//int UB = computeUB(st,ed);
			int UB = 0;
			int mySS, myCS;
			mySS = getEdgeMess(st, ed, 3);
			myCS = getEdgeMess(st, ed, 4);
			if (myCS > LB - 2) UB = LB + 1;
			else UB = LB;

			setEdgeMess(st, ed, LB, 1);
			if (LB == UB && UB == 2)continue;

			NeiNode* i = nodeList[st].firstNei;
			while (i) {
				if (i->valid && i->id != ed) {
					NeiNode* j = nodeList[ed].firstNei;
					while (j) {
						if (j->valid && j->id != st) {
							if (j->id == i->id) {
								if (i->tuss < UB) {
									PES.insert(edge(min(i->id, st), max(i->id, st), i->tuss));
									//cout<<"Insert CS:"<<getEdgeMess(min(i->id,st),max(i->id,st),4)<<endl;
								}
								if (j->tuss < UB) {
									PES.insert(edge(min(j->id, ed), max(j->id, ed), j->tuss));
									//cout<<"Insert CS:"<<getEdgeMess(min(j->id,ed),max(j->id,ed),4)<<endl;;
								}
							}
						}
						j = j->next;
					}
				}
				i = i->next;
			}
			//PES.insert(edge(min(ed,st),max(ed,st),getEdgeMess(min(ed,st),max(ed,st),1)));
		}

		Vset.clear();
		Xset.clear();
		Sset.clear();
		centerAdjust(PES);

		endTime = clock();
		tmpTime += 1000.0 * (endTime - startTime);
		maxTime = max(maxTime, tmpTime);
		maxNumberOneCycle = max(maxNumberOneCycle, cntNuumberOnecycle);
	}
	log(cntStep, maxNumberOneCycle, maxTime);
}




void Graph::SingleNodeInsert(set<pair<int, int> > edgs, int nodeId) {
	/*cout<<"Nodde::::::::::::::::::::>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<" "<<nodeId<<endl;
	set<pair<int,int> >::iterator it;
	for(it = edgs.begin();it!=edgs.end();it++){
		cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<(*it).first<<" "<<(*it).second<<endl;
	}*/
	double tmpTime = 0;
	startTime = clock();
	set<edge, cmp> PES = oneNodePES(edgs, nodeId);

	Vset.clear();
	Xset.clear();
	Sset.clear();
	endTime = clock();
	tmpTime += 1000.0 * (endTime - startTime);

	initSuperSup();
	initConstrainSup();

	startTime = clock();

	centerAdjust(PES);
	endTime = clock();
	tmpTime += 1000.0 * (endTime - startTime);
	totalTime += tmpTime;
}

set<edge, cmp> Graph::oneNodePES(set<pair<int, int> > edgs, int nodeId) {
	//cout<<nodeId<<endl;
	set<pair<int, int> >::iterator it;
	int maxLB = -1;
	for (it = edgs.begin(); it != edgs.end(); it++) {
		//cout<<(*it).first<<" "<<(*it).second<<endl; 
		addEdge((*it).first, (*it).second);
	}
	for (it = edgs.begin(); it != edgs.end(); it++) {
		//cout<<(*it).first<<" "<<(*it).second<<endl; 
		int lb = computeLB(Find((*it).first), Find((*it).second));
		maxLB = max(maxLB, lb);
	}
	//cout<<"MaxLB:"<<maxLB<<endl;
	set<edge, cmp> PES;
	for (it = edgs.begin(); it != edgs.end(); it++) {
		//cout<<(*it).first<<" "<<(*it).second<<endl; 
		int st = min(Find((*it).first), Find((*it).second));
		int ed = max(Find((*it).first), Find((*it).second));
		setEdgeMess(st, ed, maxLB, 1);
		PES.insert(edge(st, ed, maxLB));
	}
	return PES;
}

void Graph::MultNodeInsert(vector<set<pair<int, int> > > edgs, vector<int > nodeIds) {
	double maxTime;
	for (int i = 0; i < edgs.size(); i++) {
		set<pair<int, int> >::iterator it;
		for (it = edgs[i].begin(); it != edgs[i].end(); it++) {
			//cout<<(*it).first<<" "<<(*it).second<<endl; 
			MapToIndex((*it).first);
			MapToIndex((*it).second);
		}
	}
	int cntStep = 0;
	int maxNumberOneCycle = 0;
	while (!edgs.empty()) {
		double tmpTime = 0;
		//startTime = clock();
		cout << "Step " << cntStep++ << "->" << endl;
		// set<edge,cmp>::iterator iter;
		int cntNuumberOnecycle = 0;

		startTime = clock();
		vector<set<pair<int, int> > > tmpEdges;
		vector<int > tmpNodes;
		tmpEdges.push_back(edgs[0]);
		tmpNodes.push_back(nodeIds[0]);
		edgs.erase(edgs.begin() + 0);
		nodeIds.erase(nodeIds.begin() + 0);
		for (int i = 0; i < edgs.size(); i++) {
			if (noNodeCorss(tmpEdges, edgs[0])) {
				tmpEdges.push_back(edgs[0]);
				tmpNodes.push_back(nodeIds[0]);
				edgs.erase(edgs.begin() + i);
				nodeIds.erase(nodeIds.begin() + i);
				i--;
				cntNuumberOnecycle++;
			}
			else continue;
		}
		endTime = clock();
		tmpTime += 1000.0 * (endTime - startTime);
		initSuperSup();
		initConstrainSup();

		startTime = clock();
		set<edge, cmp> PES;
		for (int k = 0; k < tmpEdges.size(); k++) {
			set<pair<int, int> >::iterator it;
			int maxLB = -1;
			for (it = tmpEdges[k].begin(); it != tmpEdges[k].end(); it++) {
				//cout<<(*it).first<<" "<<(*it).second<<endl; 
				addEdge((*it).first, (*it).second);
			}
			for (it = tmpEdges[k].begin(); it != tmpEdges[k].end(); it++) {
				//cout<<(*it).first<<" "<<(*it).second<<endl; 
				int lb = computeLB(Find((*it).first), Find((*it).second));
				maxLB = max(maxLB, lb);
			}
			//cout<<"MaxLB:"<<maxLB<<endl;
			for (it = tmpEdges[k].begin(); it != tmpEdges[k].end(); it++) {
				//cout<<(*it).first<<" "<<(*it).second<<endl; 
				setEdgeMess(Find((*it).first), Find((*it).second), maxLB, 1);
				PES.insert(edge(Find((*it).first), Find((*it).second), maxLB));
			}
		}

		centerAdjust(PES);
		endTime = clock();
		tmpTime += 1000.0 * (endTime - startTime);
		maxTime = max(maxTime, tmpTime);
		maxNumberOneCycle = max(maxNumberOneCycle, cntNuumberOnecycle);
	}
	log(cntStep, maxNumberOneCycle, maxTime);
}

bool Graph::noNodeCorss(vector<set<pair<int, int> > > edgs, set<pair<int, int> > edg) {
	set<pair<int, int> >::iterator iter;
	for (iter = edg.begin(); iter != edg.end(); iter++) {
		for (int i = 0; i < edgs.size(); i++) {
			set<pair<int, int> >::iterator it;
			for (it = edgs[i].begin(); it != edgs[i].end(); it++) {
				//cout<<(*it).first<<" "<<(*it).second<<endl; 
				vector<edge> t;
				t.push_back(edge((*it).first, (*it).second, 0));
				if (!noCross(t, (*iter).first, (*iter).second)) return false;
			}
		}
	}

	return true;
}


void Graph::DeleteCenterAdjust(set<edge, cmp> PES) {
	set<edge >::iterator iter;
	for (iter = PES.begin(); iter != PES.end(); iter++) {
		while (!Stack.empty()) Stack.pop();
		pair<int, int> e_root = make_pair((*iter).st, (*iter).ed);
		int rootTruss = getEdgeMess(e_root.first, e_root.second, 1);
		Sset[e_root] = getEdgeMess(e_root.first, e_root.second, 3);
		Stack.push(e_root);
		Vset[e_root] = 1;
		while (!(Stack.empty())) {
			pair<int, int> now_e = Stack.top();
			Stack.pop();
			//cout<<Sset[now_e]<<" "<<rootTruss-2<<endl;
			if (Sset[now_e] < rootTruss - 2) {
				//cout<<"Change:"<<now_e.first<<" "<<now_e.second<<endl;
				Eliminate(now_e, rootTruss);
				int a = now_e.first;
				int b = now_e.second;
				NeiNode* i = nodeList[a].firstNei;
				while (i) {
					if (i->valid && i->id != b) {
						NeiNode* j = nodeList[b].firstNei;
						while (j) {
							if (j->valid && j->id != a) {
								if (j->id == i->id) {
									//cout<<"Turss 1("<<a<<","<<i->id<<"):"<<i->tuss<<" "<<"Truss 2("<<b<<","<<i->id<<"):"<<j->tuss<<"    rootTruss:"<<rootTruss<<endl;
									if (min(i->tuss, j->tuss) >= rootTruss) {
										pair<int, int> new_e_1 = make_pair(min(i->id, a), max(i->id, a));
										if (i->tuss == rootTruss && j->tuss > rootTruss && Vset[new_e_1] == 0) {
											Stack.push(new_e_1);
											//cout<<"Add:"<<nodeList[min(i->id,a)].id<<" "<<nodeList[max(i->id,a)].id<<endl;
											Vset[new_e_1] = 1;
											Sset[new_e_1] = Sset[new_e_1] + getEdgeMess(new_e_1.first, new_e_1.second, 3);
										}
										pair<int, int> new_e_2 = make_pair(min(i->id, b), max(i->id, b));
										if (j->tuss == rootTruss && i->tuss > rootTruss && Vset[new_e_2] == 0) {
											Stack.push(new_e_2);
											//cout<<"Add:"<<nodeList[min(i->id,b)].id<<" "<<nodeList[max(i->id,b)].id<<endl;
											Vset[new_e_2] = 1;
											Sset[new_e_2] = Sset[new_e_2] + getEdgeMess(new_e_2.first, new_e_2.second, 3);
										}
										if (j->tuss == rootTruss && i->tuss == rootTruss) {
											if (Vset[new_e_2] == 0) {
												Stack.push(new_e_2);
												//cout<<"Add:"<<nodeList[min(i->id,b)].id<<" "<<nodeList[max(i->id,b)].id<<endl;
												Vset[new_e_2] = 1;
												Sset[new_e_2] = Sset[new_e_2] + getEdgeMess(new_e_2.first, new_e_2.second, 3);
											}
											if (Vset[new_e_1] == 0) {
												Stack.push(new_e_1);
												//cout<<"Add:"<<nodeList[min(i->id,a)].id<<" "<<nodeList[max(i->id,a)].id<<endl;
												Vset[new_e_1] = 1;
												Sset[new_e_1] = Sset[new_e_1] + getEdgeMess(new_e_1.first, new_e_1.second, 3);
											}

										}
									}
								}
							}
							j = j->next;
						}
					}
					i = i->next;
				}
				//cout<<"--------------------------------------------"<<endl<<endl;
			}

		}
	}

	map<pair<int, int>, int >::iterator vit;
	for (vit = Vset.begin(); vit != Vset.end(); vit++) {
		//cout<<nodeList[(*vit).first.first].id<<" "<<nodeList[(*vit).first.second].id<<" "<<(*vit).second<<" "<<Xset[(*vit).first]<<endl;
		if (Vset[(*vit).first] == 1 && Xset[(*vit).first] == 0) {
			//cout<<nodeList[(*vit).first.first].id<<" "<<nodeList[(*vit).first.second].id<<endl;
			setEdgeMess((*vit).first.first, (*vit).first.second, getEdgeMess((*vit).first.first, (*vit).first.second, 1) - 1, 1);
		}
	}

}


void Graph::centerDelete(int stId, int edId, int model) {
	stId = Find(stId);
	edId = Find(edId);
	//cout<<stId<<" "<<edId<<" "<<nodeNum<<endl;
	int st = min(stId, edId);
	int ed = max(stId, edId);

	double tmpTime = 0;
	initSuperSup();

	startTime = clock();

	int myTruss = getEdgeMess(st, ed, 1);
	set<edge, cmp> PES;
	NeiNode* i = nodeList[st].firstNei;
	while (i) {
		if (i->valid && i->id != ed) {
			NeiNode* j = nodeList[ed].firstNei;
			while (j) {
				if (j->valid && j->id != st) {
					if (j->id == i->id) {
						if (i->tuss <= myTruss) PES.insert(edge(min(i->id, st), max(i->id, st), i->tuss));
						if (j->tuss <= myTruss) PES.insert(edge(min(j->id, ed), max(j->id, ed), j->tuss));
					}
				}
				j = j->next;
			}
		}
		i = i->next;
	}

	deleteEdge(st, ed);

	Vset.clear();
	Xset.clear();
	Sset.clear();

	DeleteCenterAdjust(PES);
	endTime = clock();
	tmpTime += 1000.0 * (endTime - startTime);
	totalTime += tmpTime;
	//log(1,1,time);

	//cout<<"Time of insert One edge:"<<time<<endl;
}

void Graph::centerMultDelete(vector<int> stIds, vector<int > edIds) {
	//initSuperSup();
	//initConstrainSup();
	double maxTime;
	vector<edge> insertEdge;
	for (int i = 0; i < stIds.size(); i++) {
		//if(!exist_edge(stIds[i],edIds[i])){
		//    cout<<"No exiSt"<<endl;
		//    continue;
		//}
		edge tmp(stIds[i], edIds[i], 2);
		insertEdge.push_back(tmp);
	}
	int cntStep = 0;
	int maxNumberOneCycle = 0;
	while (!insertEdge.empty()) {
		initSuperSup();

		double tmpTime = 0;
		startTime = clock();
		cout << "Step " << cntStep++ << "->" << endl;
		// set<edge,cmp>::iterator iter;
		int cntNuumberOnecycle = 0;
		vector<edge > nowInsertEdgeSet;
		int a = Find(insertEdge[0].st);
		int b = Find(insertEdge[0].ed);
		int st = min(a, b);
		int ed = max(a, b);
		deleteEdge(st, ed);
		nowInsertEdgeSet.push_back(insertEdge[0]);

		//int sss = 0;
		insertEdge.erase(insertEdge.begin() + 0);
		for (int k = 0; k < insertEdge.size(); k++) {
			if (noCross(nowInsertEdgeSet, insertEdge[k].st, insertEdge[k].ed)) {
				//cout<<sss<<endl;
				int a = Find(insertEdge[k].st);
				int b = Find(insertEdge[k].ed);
				int st = min(a, b);
				int ed = max(a, b);
				deleteEdge(st, ed);
				nowInsertEdgeSet.push_back(insertEdge[k]);
				//cout<<"fuck:"<<k<<" "<<insertEdge.size()<<" "<<insertEdge[k].st<<" "<<insertEdge[k].ed<<" "<<a<<" "<<b<<endl; 
				insertEdge.erase(insertEdge.begin() + k);
				k--;
				cntNuumberOnecycle++;
				//cout<<sss++<<endl;
			}
			else continue;
		}
		endTime = clock();
		tmpTime += 1000.0 * (endTime - startTime);
		initSuperSup();

		startTime = clock();
		set<edge, cmp> PES;
		//cout<<nowInsertEdgeSet.size()<<endl;
		for (int k = 0; k < nowInsertEdgeSet.size(); k++) {
			int a = Find(nowInsertEdgeSet[k].st);
			int b = Find(nowInsertEdgeSet[k].ed);
			//cout<<nowInsertEdgeSet[k].st<<" "<<nowInsertEdgeSet[k].ed<<" "<<a<<" "<<b<<endl;
			int st = min(a, b);
			int ed = max(a, b);

			int myTruss = getEdgeMess(st, ed, 1);
			NeiNode* i = nodeList[st].firstNei;
			while (i) {
				if (i->valid && i->id != ed) {
					NeiNode* j = nodeList[ed].firstNei;
					while (j) {
						if (j->valid && j->id != st) {
							if (j->id == i->id) {
								if (i->tuss <= myTruss) {
									PES.insert(edge(min(i->id, st), max(i->id, st), i->tuss));
									//cout<<"Insert CS:"<<getEdgeMess(min(i->id,st),max(i->id,st),4)<<endl;
								}
								if (j->tuss <= myTruss) {
									PES.insert(edge(min(j->id, ed), max(j->id, ed), j->tuss));
									//cout<<"Insert CS:"<<getEdgeMess(min(j->id,ed),max(j->id,ed),4)<<endl;;
								}
							}
						}
						j = j->next;
					}
				}
				i = i->next;
			}
			//PES.insert(edge(min(ed,st),max(ed,st),getEdgeMess(min(ed,st),max(ed,st),1)));
		}

		Vset.clear();
		Xset.clear();
		Sset.clear();
		DeleteCenterAdjust(PES);

		endTime = clock();
		tmpTime += 1000.0 * (endTime - startTime);
		maxTime = max(maxTime, tmpTime);
		maxNumberOneCycle = max(maxNumberOneCycle, cntNuumberOnecycle);
	}
	log(cntStep, maxNumberOneCycle, maxTime);
}

/*
	for(int i=0;i<stIds.size();i++){
		int a = Find(stIds[i]);
		int b = Find(edIds[i]);
		int st = min(a,b);
		int ed = max(a,b);
		int LB = computeLB(st,ed);
		int UB=0;
		//for(int i=0;i<v.size();i++) cout<<v[i]<<" ";
		//cout<<endl<<v.size()<<endl;
		setEdgeMess(st,ed,LB,1);

		int mySS = computeSuperSup(st,ed,LB);
		int myCS = computeConstrainSup(st,ed,LB);
		//cout<<">>"<<mySS<<" "<<myCS<<endl;
		if(myCS > LB-2) UB = LB+1;
		else UB = LB;
		int small = 999999999;
		NeiNode *ii = nodeList[st].firstNei;
		while(ii){
			if(ii->valid && ii->id != ed){
				NeiNode *j = nodeList[ed].firstNei;
				while(j){
					if(j->valid && j->id!= st){
						if(j->id == ii->id) {
							int a = ii->tuss;// getTuss(i->id,st);
							int b = j->tuss;//getTuss(j->id,ed);
							int d = min(a,b);
							small = min(small,d);
						}
					}
					j = j->next;
				}
			}
			ii = ii->next;
		}
		edge tmp(small,UB-1,0);
		tmp.st2 = stIds[i];
		tmp.ed2 = edIds[i];
		insertEdge.push_back(tmp);
		//cout<<st<<" + "<<ed<<" UB:"<<UB<<" LB:"<<LB<<endl;
	}

	sort(insertEdge.begin(),insertEdge.end(),LessSort);
	int cntStep = 0;
	int maxNumberOneCycle = 0;

	while(!insertEdge.empty()){
		startTime = clock();
		cout<<"Step "<<cntStep++<<" :"<<endl;
	   // set<edge,cmp>::iterator iter;
		int lastEnd = -1;
		int cntNuumberOnecycle = 0;
		for(int k=0;k<insertEdge.size();k++){
			if(lastEnd<insertEdge[k].st) {
				lastEnd = insertEdge[k].ed;
				addEdge(insertEdge[k].st2,insertEdge[k].ed2);
			}else continue;
		}
		initSuperSup();
		initSuperSup();
		for(int k=0;k<insertEdge.size();k++){
			if(lastEnd<insertEdge[k].st) {
				lastEnd = insertEdge[k].ed;
				cout<<"Insert:"<<insertEdge[k].st2<<" "<<insertEdge[k].ed2<<endl;
				int st = Find(insertEdge[k].st2);
				int ed = Find(insertEdge[k].ed2);
				int LB = computeLB(st,ed);
				int UB=0;
				setEdgeMess(st,ed,LB,1);
				int mySS,myCS;
				mySS = getEdgeMess(st,ed,3);
				myCS = getEdgeMess(st,ed,4);
				if(myCS > LB-2) UB = LB+1;
				else UB = LB;
				centerAdjust(PES);

				//centerInsert(insertEdge[k].st2,insertEdge[k].ed2,1);
				insertEdge.erase(insertEdge.begin()+k);
				k--;
				cntNuumberOnecycle ++;
			}else continue;
		}
		endTime = clock();
		maxTime = max(maxTime,1000.0*(endTime - startTime));
		maxNumberOneCycle = max(maxNumberOneCycle,cntNuumberOnecycle);
		cntStep++;
	}
	log(cntStep,maxNumberOneCycle,maxTime);
	cout<<"Number of cycles:"<<cntStep<<endl;
	cout<<"Time of longest cycle:"<<maxTime<<endl;
*/
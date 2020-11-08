#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include <fstream>
using namespace std;

int main() {
	FILE* fp;
	fp = fopen("H:/æ∫»¸œ‡πÿ/2020ccf-kmax/ktruss-data/s18.e16.rmat.edgelist.tsv", "r");
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
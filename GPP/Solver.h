#ifndef SOLVER_H
#define SOLVER_H

#include <ctime>
#include <vector>
#include <map>
#include <unordered_set>

#include "Graph.h"
#include "Config.h"
#include "Instance.h"

using namespace std;


struct TabuStruct {
	int partNum, maxDegree;

	int level, nodeNum;
	int curObj, bestObj;
	shared_ptr<Graph> pg;

	vector<int> curSol;
	vector<int> mvFreq;
	vector<vector<int>> tabuList;

	vector<int> bvNum4P;
	vector<int> partWgt;
	vector<int> maxGainIdx4P;

	vector<int> inWgt4V;
	vector<vector<int>> extWgt4V2P;
	vector<vector<unordered_set<int>>> gainNode4P;


	void init(int partnum, int maxdegree) {
		partNum = partnum;
		maxDegree = maxdegree;
		gainNode4P = vector<vector<unordered_set<int>>>(partNum, vector<unordered_set<int>>(2 * maxDegree + 1));
	}

	void reset(int lev, const shared_ptr<Graph> &g, const vector<int> &initSol) {
		level = lev;
		nodeNum = g->nodes.size();
		curObj = bestObj = 0;
		pg = g;

		curSol.assign(initSol.begin(), initSol.end());
		mvFreq.assign(nodeNum, 0);
		tabuList.assign(nodeNum, vector<int>(partNum, 0));
		bvNum4P.assign(partNum, 0);
		partWgt.assign(partNum, 0);
		maxGainIdx4P.assign(partNum, 0);
		inWgt4V.assign(nodeNum, 0);
		extWgt4V2P.assign(nodeNum, vector<int>(partNum, 0));
		for (auto &bktArr : gainNode4P) {
			for (auto &bkt : bktArr) {
				bkt.clear();
			}
		}

		load();
	}

	void load() {
		for (int v = 0; v < nodeNum; ++v) {
			partWgt[curSol[v]] += pg->nodes[v].wgt;
			auto &adjList(pg->nodes[v].adjList);
			for (auto &adj : adjList) {
				if (curSol[adj.vid] == curSol[v]) {
					inWgt4V[v] += adj.wgt;
				}
				else {
					extWgt4V2P[v][curSol[adj.vid]] += adj.wgt;
					curObj += adj.wgt;
				}
			}

			for (int k = 0; k < partNum; ++k) {
				if (extWgt4V2P[v][k] != 0 && curSol[v] != k) {
					++bvNum4P[k];
					int gainIndex = extWgt4V2P[v][k] - inWgt4V[v] + maxDegree;
					gainNode4P[k][gainIndex].insert(v);
					if (gainIndex > maxGainIdx4P[k]) { maxGainIdx4P[k] = gainIndex; }
				}
			}
		}
		curObj /= 2;
		bestObj = curObj;
	}

	void insert(int gain, int vertex, int part) {
		++bvNum4P[part];
		gainNode4P[part][gain + maxDegree].insert(vertex);
	}

	void erase(int gain, int vertex, int part) {
		gainNode4P[part][gain + maxDegree].erase(vertex);
		--bvNum4P[part];
	}
};

class Solver {
public:
	void input();
	void output(int runId);

	void init();
	void coarsen();
	void uncoarsen();

	void initialPartition();
	void its(int level, vector<int> &partition, const shared_ptr<Graph> &pg);
	int getFirstPart(TabuStruct &tss);
	int getSecondPart(TabuStruct &tss, int firstPart);
	int firstNode2Mv(int iter, TabuStruct &tss, int tarPart);
	int secondNode2Mv(int iter, TabuStruct &tss, int tarPart, int firstPart);
	void updateGain(int iter, TabuStruct &tss, int vertex, int tarPart);
	void mvOneNode(int iter, TabuStruct &tss, int vertex, int tarPart);
	void mvTwoNode(int iter, TabuStruct &tss, int firstNode, int firstPart);
	void perturbe(int iter, TabuStruct &tss);

	int getBalanceQuality(TabuStruct &tss);
	double getImbalance();
	double getImbalance(const vector<int> &partWgts);
	size_t check(vector<int> &sol, const shared_ptr<Graph> &pg);


	Instance ins;
	Config cfg;

	vector<int> curNodes;
	map<pair<int, int>, int> curEdges;
	vector<shared_ptr<Graph>> pGraphs;
	vector<vector<int>> nodeMap;

	TabuStruct tss;
	vector<int> partition;

	clock_t startTime;
	size_t seed;
};


#endif // !SOLVER_H

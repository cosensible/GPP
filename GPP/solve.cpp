#include "Solver.h"

#include <fstream>
#include <sstream>
#include <iostream>

#include <string>
#include <memory>
#include <numeric>
#include <algorithm>
#include <unordered_set>

using namespace std;

void Solver::input() {
	string infile = ins.insPath + ins.insName + ".graph";
	ifstream ifs(infile);
	ifs >> ins.nodeNum >> ins.edgeNum;
	ifs.get();
	ins.isolated = 0;

	string line;
	int adj = -1;
	for (int v = 0; v < ins.nodeNum; ++v) {
		getline(ifs, line);
		if (line.size() == 0) { ++ins.isolated; }
		istringstream iss(line);
		while (iss >> adj) {
			curEdges[{v, adj - 1}] = 1;
		}
	}
	curNodes.assign(ins.nodeNum, 1);
}

void Solver::init() {
	seed = 1608020520;
	//seed = time(0);
	srand(seed);
	startTime = clock();

	// 根据算例属性设置参数 cfg, ins
	tss.init(ins.partNum, ins.maxDegree);
}

void Solver::coarsen() {
	while (curNodes.size() > cfg.cst + ins.isolated) {
		// 构造当前图的邻接表
		pGraphs.emplace_back(new Graph(curNodes, curEdges));
		nodeMap.emplace_back(curNodes.size());

		unordered_set<int> unmatchSet; // 该层尚未匹配的节点集合
		for (int i = 0; i < curNodes.size(); ++i) { unmatchSet.insert(i); }
		int newNodeNum = 0; // 新节点个数
		while (!unmatchSet.empty()) {
			// 随机挑选未匹配节点 nid
			auto it = unmatchSet.begin();
			for (int i = rand() % unmatchSet.size(); i > 0; --i, ++it) {}
			int nid = *it;
			unmatchSet.erase(it);

			auto &adjList(pGraphs.back()->nodes[nid].adjList);
			for (auto it = adjList.begin(); it != adjList.end(); ++it) {
				if (unmatchSet.count(it->vid) != 0) { // 邻居未匹配
					unmatchSet.erase(it->vid);
					nodeMap.back()[it->vid] = newNodeNum;
					break;
				}
			}
			nodeMap.back()[nid] = newNodeNum++;
		}

		// 根据最大匹配（节点映射关系）压缩图
		vector<int> newNodes(newNodeNum, 0);
		map<pair<int, int>, int> newEdges;
		for (int i = 0; i < curNodes.size(); ++i) {
			newNodes[nodeMap.back()[i]] += curNodes[i];
		}
		for (auto e = curEdges.begin(); e != curEdges.end(); ++e) {
			int beg = nodeMap.back()[e->first.first], end = nodeMap.back()[e->first.second];
			if (beg == end) { continue; }
			newEdges[{beg, end}] += e->second;
		}
		swap(curNodes, newNodes);
		swap(curEdges, newEdges);
	}
	pGraphs.emplace_back(new Graph(curNodes, curEdges));
	//pGraphs.back()->printGraph();
}


void Solver::initialPartition() {
	int nodeNum = pGraphs.back()->nodes.size();
	vector<int> idx(nodeNum);
	iota(idx.begin(), idx.end(), 0);
	random_shuffle(idx.begin(), idx.end());

	//for (int i = idx.size() - 1; i > 0; --i) {
	//	int index = rand() % (i + 1);
	//	swap(idx[i], idx[index]);
	//}

	partition.assign(nodeNum, -1);
	for (int i = 0; i < nodeNum; ++i) {
		partition[idx[i]] = i % ins.partNum;
	}
}

void Solver::uncoarsen() {
	initialPartition();
	cout << "level=" << pGraphs.size() - 1 << "\t";
	//cout << "before its : " << check(partition, pGraphs.back()) << ", ";
	its(pGraphs.size() - 1, partition, pGraphs.back());
	//cout << "after its : " << check(partition, pGraphs.back()) << endl;

	for (int level = pGraphs.size() - 2; level >= 0; --level) {
		int nodeNum = pGraphs[level]->nodes.size();
		vector<int> partition1(nodeNum, -1);
		for (int v = 0; v < nodeNum; ++v) {
			partition1[v] = partition[nodeMap[level][v]];
		}
		cout << "level=" << level << "\t";
		//cout << "before its : " << check(partition1, pGraphs[level]) << ", ";
		its(level, partition1, pGraphs[level]);
		swap(partition, partition1);
		//cout << "after its : " << check(partition, pGraphs[level]) << endl;
	}
	cout << "ins=" << ins.insName << ", partnum=" << ins.partNum << ", imbalance=" << getImbalance()
		<< ", seed=" << seed << ", duration=" << (clock() - startTime) / static_cast<double>(CLOCKS_PER_SEC) << endl << endl;
}

int Solver::getFirstPart(TabuStruct &tss) {
	int maxPart = -1, maxWgt4P = 0;
	int randPart = rand() % ins.partNum;
	for (int k = 0; k < ins.partNum; ++k) {
		if (tss.partWgt[k] > maxWgt4P) {
			maxPart = k;
			maxWgt4P = tss.partWgt[k];
		}
	}
	while (randPart == maxPart) {
		randPart = rand() % ins.partNum;
	}
	return randPart;
}

int Solver::getSecondPart(TabuStruct &tss, int firstPart) {
	int maxPart = -1, maxWgt4P = 0;
	int randPart = rand() % ins.partNum;
	for (int k = 0; k < ins.partNum; ++k) {
		if (tss.partWgt[k] > maxWgt4P) {
			maxPart = k;
			maxWgt4P = tss.partWgt[k];
		}
	}
	while (randPart == maxPart || randPart == firstPart) {
		randPart = rand() % ins.partNum;
	}
	return randPart;
}

int Solver::firstNode2Mv(int iter, TabuStruct &tss, int tarPart) {
	int gainIndex = tss.maxGainIdx4P[tarPart];
	int minDiff = INT_MAX, minVal = INT_MAX;
	bool notFound = true;
	int vertex = -1;

	while (notFound) {
		auto &bkt(tss.gainNode4P[tarPart][gainIndex]);
		for (int adj : bkt) {
			if (tss.partWgt[tss.curSol[adj]] <= tss.partWgt[tarPart]) { continue; }
			int gain = tss.extWgt4V2P[adj][tarPart] - tss.inWgt4V[adj];
			if (tss.tabuList[adj][tarPart] <= iter || tss.curObj - gain < tss.bestObj) {
				int diff = abs(tss.partWgt[tarPart] + 2 * tss.pg->nodes[adj].wgt - tss.partWgt[tss.curSol[adj]]);
				int val = gain + tss.mvFreq[adj] / (iter + 1);
				if (minDiff > diff && minVal > val) {
					minDiff = diff;
					minVal = val;
					vertex = adj;
				}
			}
		}
		if (vertex == -1) { --gainIndex; }
		else { notFound = false; }
		if (gainIndex <= 20) break;
	}
	return vertex;
}

int Solver::secondNode2Mv(int iter, TabuStruct &tss, int tarPart, int firstPart) {
	int gainIndex = tss.maxGainIdx4P[tarPart];
	int minDiff = INT_MAX, minVal = INT_MAX;
	bool notFound = true;
	int vertex = -1;

	while (notFound) {
		auto &bkt(tss.gainNode4P[tarPart][gainIndex]);
		for (int adj : bkt) {
			if (tss.curSol[adj] == firstPart) { continue; }
			int gain = tss.extWgt4V2P[adj][tarPart] - tss.inWgt4V[adj];
			if (tss.tabuList[adj][tarPart] <= iter || tss.curObj - gain < tss.bestObj) {
				int diff = abs(tss.partWgt[tarPart] + 2 * tss.pg->nodes[adj].wgt - tss.partWgt[tss.curSol[adj]]);
				int val = gain + tss.mvFreq[adj] / (iter + 1);
				if (minDiff > diff && minVal > val) {
					minDiff = diff;
					minVal = val;
					vertex = adj;
				}
			}
		}
		if (vertex == -1) { --gainIndex; }
		else { notFound = false; }
		if (gainIndex <= 20) break;
	}
	return vertex;
}

void Solver::updateGain(int iter, TabuStruct &tss, int vertex, int tarPart) {
	int srcPart = tss.curSol[vertex];
	tss.tabuList[vertex][srcPart] = iter + tss.bvNum4P[srcPart] * 0.3 + rand() % 3;
	// 删除待移动节点到各分区的桶内 gain 值节点
	for (int k = 0; k < ins.partNum; ++k) {
		if (tss.extWgt4V2P[vertex][k] != 0 && tss.curSol[vertex] != k) {
			tss.erase(tss.extWgt4V2P[vertex][k] - tss.inWgt4V[vertex], vertex, k);
		}
	}
	// 更新节点到源分区和目标分区的外部边权和、节点的内部边权和
	tss.extWgt4V2P[vertex][srcPart] = tss.inWgt4V[vertex];
	tss.inWgt4V[vertex] = tss.extWgt4V2P[vertex][tarPart];
	tss.extWgt4V2P[vertex][tarPart] = 0;
	tss.curSol[vertex] = tarPart;
	// 因为节点的内部边权和被改变, 所以节点到所有分区的gain值被改变
	for (int k = 0; k < ins.partNum; ++k) {
		if (tss.extWgt4V2P[vertex][k] != 0 && tss.curSol[vertex] != k) {
			int gain = tss.extWgt4V2P[vertex][k] - tss.inWgt4V[vertex];
			if (gain + ins.maxDegree > tss.maxGainIdx4P[k]) {
				tss.maxGainIdx4P[k] = gain + ins.maxDegree;
			}
			tss.insert(gain, vertex, k);
		}
	}
	++tss.mvFreq[vertex];
}

void Solver::mvOneNode(int iter, TabuStruct &tss, int vertex, int tarPart) {
	if (vertex != -1 && tss.curSol[vertex] != tarPart) {
		int srcPart = tss.curSol[vertex];
		tss.partWgt[tarPart] += tss.pg->nodes[vertex].wgt;
		tss.partWgt[srcPart] -= tss.pg->nodes[vertex].wgt;

		// 更新待移动节点到各分区的 gain 值
		updateGain(iter, tss, vertex, tarPart);
		// 更新邻居节点到各分区的 gain 值
		auto &adjList = tss.pg->nodes[vertex].adjList;
		for (auto &adj : adjList) {
			for (int k = 0; k < ins.partNum; ++k) {
				if (tss.extWgt4V2P[adj.vid][k] != 0 && tss.curSol[adj.vid] != k) {
					tss.erase(tss.extWgt4V2P[adj.vid][k] - tss.inWgt4V[adj.vid], adj.vid, k);
				}
			}
			if (tss.curSol[adj.vid] == tarPart) {
				tss.inWgt4V[adj.vid] += adj.wgt;
				tss.extWgt4V2P[adj.vid][srcPart] -= adj.wgt;
				tss.curObj -= adj.wgt;
			}
			else if (tss.curSol[adj.vid] == srcPart) {
				tss.inWgt4V[adj.vid] -= adj.wgt;
				tss.extWgt4V2P[adj.vid][tarPart] += adj.wgt;
				tss.curObj += adj.wgt;
			}
			else {
				tss.extWgt4V2P[adj.vid][srcPart] -= adj.wgt;
				tss.extWgt4V2P[adj.vid][tarPart] += adj.wgt;
			}
			for (int k = 0; k < ins.partNum; ++k) {
				if (tss.extWgt4V2P[adj.vid][k] != 0 && tss.curSol[adj.vid] != k) {
					int gain = tss.extWgt4V2P[adj.vid][k] - tss.inWgt4V[adj.vid];
					if (gain + ins.maxDegree > tss.maxGainIdx4P[k]) {
						tss.maxGainIdx4P[k] = gain + ins.maxDegree;
					}
					tss.insert(gain, adj.vid, k);
				}
			}
		}
	}
}

void Solver::mvTwoNode(int iter, TabuStruct &tss, int firstNode, int firstPart) {
	int secondPart = getSecondPart(tss, firstPart);
	int secondNode = secondNode2Mv(iter, tss, secondPart, firstPart);
	if (secondNode != -1 && firstNode != -1) {
		mvOneNode(iter, tss, firstNode, firstPart);
		mvOneNode(iter, tss, secondNode, secondPart);
	}
}

int Solver::getBalanceQuality(TabuStruct &tss) {
	int bq = 0;
	int div = ins.nodeNum / ins.partNum;
	int vmp = ins.nodeNum % ins.partNum;
	for (int k = 0; k < ins.partNum; ++k) {
		if (vmp != 0 && tss.partWgt[k] != (div + 1) && tss.partWgt[k] != div)
			bq += abs(div - tss.partWgt[k]);
		else if (vmp == 0)
			bq += abs(div - tss.partWgt[k]);
	}
	return bq;
}

void Solver::perturbe(int iter, TabuStruct &tss) {
	int strength = tss.nodeNum * 0.02;
	for (int t = 0; t < strength; ++t) {
		int randPart = getFirstPart(tss);
		int vertex = rand() % tss.nodeNum;
		if (tss.level != 0) {
			while (tss.curSol[vertex] == randPart || tss.partWgt[tss.curSol[vertex]] < tss.partWgt[randPart]) {
				vertex = rand() % tss.nodeNum;
			}
		}
		else {
			while (tss.curSol[vertex] == randPart) {
				if (tss.partWgt[tss.curSol[vertex]] - tss.partWgt[randPart] >= 2 * tss.pg->nodes[vertex].wgt) {
					break;
				}
				vertex = rand() % tss.nodeNum;
			}
		}
		mvOneNode(iter, tss, vertex, randPart);
	}
}

void Solver::its(int level, vector<int> &partition, const shared_ptr<Graph> &pg) {
	tss.reset(level, pg, partition);
	//cout << "Before ITS, obj=" << tss.bestObj << ", ";

	const size_t maxIter = 10 * tss.nodeNum;
	const size_t smlOptCnt = 20, bigOptCnt = 0.1*tss.nodeNum;
	size_t localOpt1 = 0, localOpt2 = 0, noImpv = 0;
	size_t bestBalance = INT_MAX;

	for (size_t iter = 0; iter < maxIter; ++iter) {
		int firstPart = getFirstPart(tss);
		int firstNode = firstNode2Mv(iter, tss, firstPart);
		if (localOpt1 < smlOptCnt) {
			if (tss.bestObj < tss.curObj) { ++localOpt1; }
			mvOneNode(iter, tss, firstNode, firstPart);
		}
		else {
			if (tss.bestObj < tss.curObj) { ++localOpt2; }
			mvTwoNode(iter, tss, firstNode, firstPart);
			if (localOpt2 == smlOptCnt) {
				localOpt1 = localOpt2 = 0;
			}
		}

		int bq = getBalanceQuality(tss);
		if (tss.bestObj >= tss.curObj && bq <= bestBalance) {
			tss.bestObj = tss.curObj;
			partition.assign(tss.curSol.begin(), tss.curSol.end());
			noImpv = 0;
		}
		if (bq < bestBalance && tss.curObj == tss.bestObj) {
			bestBalance = bq;
		}

		noImpv++;
		if (iter % bigOptCnt == 0) {
			perturbe(iter, tss);
		}
	}
	cout << "obj=" << tss.bestObj << endl;
}

size_t Solver::check(vector<int> &sol, const shared_ptr<Graph> &pg) {
	size_t obj = 0;
	for (int v = 0; v < sol.size(); ++v) {
		for (auto &adj : pg->nodes[v].adjList) {
			if (sol[v] != sol[adj.vid]) {
				obj += adj.wgt;
			}
		}
	}
	return obj / 2;
}

double Solver::getImbalance() {
	vector<int> partWgts(ins.partNum, 0);
	for (int i = 0; i < ins.nodeNum; ++i) {
		partWgts[partition[i]] += pGraphs[0]->nodes[i].wgt;
	}
	int maxPartWgt = *std::max_element(partWgts.begin(), partWgts.end());
	int wgtSum = std::accumulate(partWgts.begin(), partWgts.end(), 0);
	int optWgt = 1 + wgtSum / ins.partNum;
	return 1.0*maxPartWgt / optWgt;
}

double Solver::getImbalance(const vector<int> &partWgts) {
	int maxPartWgt = *std::max_element(partWgts.begin(), partWgts.end());
	int wgtSum = std::accumulate(partWgts.begin(), partWgts.end(), 0);
	int optWgt = 1 + wgtSum / ins.partNum;
	return 1.0*maxPartWgt / optWgt;
}

void Solver::output(int runId) {
	double duration = (clock() - startTime) / static_cast<double>(CLOCKS_PER_SEC);
	string outfile = ins.solPath + ins.insName + "_p" + to_string(ins.partNum) + "." + to_string(runId) + ".sol";

	vector<int> partWgts(ins.partNum, 0);
	for (int i = 0; i < ins.nodeNum; ++i) {
		partWgts[partition[i]] += pGraphs[0]->nodes[i].wgt;
	}

	ofstream out(outfile);
	out << "Instance=" << ins.insName << " P=" << ins.partNum << " Obj=" << check(partition, pGraphs[0])
		<< " Imbalance=" << getImbalance(partWgts) << " Duration=" << duration << " Seed=" << seed << endl << endl;

	out << "wgt of parts:" << endl;
	for (int wgt : partWgts) { out << wgt << " "; }
	out << endl << endl;

	out << "solution:" << endl;
	for (int k : partition) { out << k << " "; }
	out << endl;

	out.close();
}
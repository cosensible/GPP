#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <set>
#include <map>

using namespace std;

struct AdjNode {
	int vid;
	int wgt;
	AdjNode(int id, int wgt) :vid(id), wgt(wgt) {}
};

struct Comp {
	bool operator()(const AdjNode &lhs, const AdjNode &rhs) {
		if (lhs.wgt != rhs.wgt) return lhs.wgt > rhs.wgt;
		return lhs.vid > rhs.vid;
	}
};

struct Node {
	int wgt;
	set<AdjNode, Comp> adjList;
	Node() = default;
	Node(int wgt) :wgt(wgt) {}
};

struct Graph {
	vector<Node> nodes;

	Graph(const vector<int> &ns, const map<pair<int, int>, int> &es) {
		nodes.resize(ns.size());
		for (size_t i = 0; i < ns.size(); ++i) {
			nodes[i].wgt = ns[i];
		}
		for (auto e = es.begin(); e != es.end(); ++e) {
			int src = e->first.first, dst = e->first.second, wgt = e->second;
			nodes[src].adjList.insert({ dst,wgt });
		}
	}

	void printGraph() {
		for (const Node &node : nodes) {
			cout << node.wgt << " ";
			for (const AdjNode &adj : node.adjList) {
				cout << adj.vid << " " << adj.wgt << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
};

#endif // !GRAPH_H

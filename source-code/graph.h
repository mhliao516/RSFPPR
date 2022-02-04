#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <queue>
using namespace std;

bool maxcmp(const pair<int, int>& a, const pair<int, int>& b){
    return a.first > b.first;
}

bool maxcmp1(const pair<double, pair<int, int>>& a, const pair<double, pair<int, int>>& b){
    return a.first > b.first;
}

class graph{
public:
    int n;
    long m;
    int** AdjList;
    double** edgeWeight;
    int* degree;
    double* weightdegree;

    int** alias_table;
    double** alias_probability;

    int* degree_order;

    graph()
    {
        n = m = 0;
    }
    ~graph()
    {
        for(int i=0; i<n; i++) {
            delete[] AdjList[i];
            delete[] edgeWeight[i];
        }
        delete[] AdjList;
        delete[] edgeWeight;
        delete[] degree;
        delete[] weightdegree;
    }

    void read_graph_weighted(string fn) {
        string filename = "data/" + fn + ".txt";
        ifstream infile(filename.c_str());
        infile >> n;
        cout << "n: " << n << endl;

        degree = new int[n];
        weightdegree = new double[n];
        double avg_d = 0.;
        for(int i=0; i<n; i++) {
            degree[i]=0;
            weightdegree[i] = 0.;
        }
        int u,v;
        double w;
        while(infile >> u >> v >> w) {
            assert((u<n)&&(v<n));
            degree[u]++;
            degree[v]++;
            weightdegree[u] += w;
            weightdegree[v] += w;
        }

        // for(int i=0; i<n; i++) {
        //     if(degree[i]==0) cout << "degree 0: " << i << endl;
        // }

        AdjList = new int*[n];
        edgeWeight = new double*[n];
        for(int i=0; i<n; i++) {
            AdjList[i] = new int[degree[i]];
            edgeWeight[i] = new double[degree[i]];
        }
        int* ptr = new int[n];
        for(int i=0; i<n; i++) ptr[i]=0;

        infile.clear();
        infile.seekg(0);

        infile >> n;
        while(infile >> u >> v >> w) {
            AdjList[u][ptr[u]] = v;
            AdjList[v][ptr[v]] = u;
            m++;
            edgeWeight[u][ptr[u]++] = w;
            edgeWeight[v][ptr[v]++] = w;
        }

        for(int i=0; i<n; i++) {
            avg_d += weightdegree[i]/n;
        }

        cout << "m: " << m << endl;
        cout << "avg degree: " << (double)2*m/n << endl;
        cout << "avg d: " << avg_d << endl;

        infile.close();
        delete[] ptr;
    }

    void read_graph(string fn) {
        string filename = "data/" + fn + ".txt";
        ifstream infile(filename.c_str());
        infile >> n;
        cout << "n: " << n << endl;

        degree = new int[n];
        weightdegree = new double[n];
        for(int i=0; i<n; i++) {
            degree[i]=0;
            weightdegree[i] = 0.;
        }
        int u,v;
        while(infile >> u >> v) {
            assert((u<n)&&(v<n));
            degree[u]++;
            degree[v]++;
            weightdegree[u] += 1.0;
            weightdegree[v] += 1.0;
        }

        AdjList = new int*[n];
        edgeWeight = new double*[n];
        for(int i=0; i<n; i++) {
            AdjList[i] = new int[degree[i]];
            edgeWeight[i] = new double[degree[i]];
        }
        int* ptr = new int[n];
        for(int i=0; i<n; i++) ptr[i]=0;

        infile.clear();
        infile.seekg(0);

        infile >> n;
        while(infile >> u >> v) {
            AdjList[u][ptr[u]] = v;
            AdjList[v][ptr[v]] = u;
            m++;
            edgeWeight[u][ptr[u]++] = 1.0;
            edgeWeight[v][ptr[v]++] = 1.0;
        }
        cout << "m: " << m << endl;

        infile.close();
        delete[] ptr;
    }

    void alias_preprocess() {
        alias_table = new int*[n];
        alias_probability = new double*[n];
        for(int i=0; i<n; i++) {
            alias_table[i] = new int[degree[i]];
            alias_probability[i] = new double[degree[i]];
        }
        queue<int> small;
        queue<int> large;

        for(int i=0; i<n; i++) {
            double* src = new double[degree[i]];
            for(int j=0; j<degree[i]; j++) {
                src[j] = edgeWeight[i][j];
            }
            double weight_avg = weightdegree[i]/(double)degree[i];
            for(int j=0; j<degree[i]; j++) {
                src[j] /= weight_avg;
                if(src[j]<1) small.push(j);
                else large.push(j);
            }
            while((!small.empty()) && (!large.empty())) {
                int small_idx = small.front();
                int large_idx = large.front();
                alias_table[i][small_idx] = large_idx;
                alias_probability[i][small_idx] = src[small_idx];
                small.pop();
                src[large_idx] = src[large_idx] + src[small_idx] - 1;
                if(src[large_idx]<1) {
                    large.pop();
                    small.push(large_idx);
                }
            }
            while(!small.empty()) {
                int small_idx = small.front();
                small.pop();
                alias_table[i][small_idx] = small_idx;
                alias_probability[i][small_idx] = 1.0;
            }
            while(!large.empty()) {
                int large_idx = large.front();
                large.pop();
                alias_table[i][large_idx] = large_idx;
                alias_probability[i][large_idx] = 1.0;
            }
        }

        cout << "alias preprocess finish" << endl;
    }

    void compute_degree_order() {
        degree_order = new int[n];
        vector<pair<int, int>> node_deg;
        for(int i=0; i<n; i++) {
            node_deg.push_back(make_pair(weightdegree[i],i));
        }
        sort(node_deg.begin(), node_deg.end(), maxcmp);
        for(int i=0; i<n; i++) {
            degree_order[i] = node_deg[i].second;
        }
    }

    int* get_degree_order() {
        return degree_order;
    }

    void sortAdjList() {
        for(int i=0; i<n; i++) {
            pair<double, pair<int, int>>* templist;
            templist = new pair<double, pair<int, int>>[degree[i]];
            for(int j=0; j<degree[i]; j++) {
                int tempNbr = j;
                int tempNbrID = AdjList[i][j];
                templist[j] = make_pair(weightdegree[tempNbr], make_pair(tempNbr, tempNbrID));
            }
            sort(templist, (templist+degree[i]), maxcmp1);
            for(int j=0; j<degree[i]; j++) {
                AdjList[i][j] = templist[j].second.second;
                templist[j].first = edgeWeight[i][templist[j].second.first];
            }
            for(int j=0; j<degree[i]; j++) {
                edgeWeight[i][j] = templist[j].first;
            }
        }
    }
};
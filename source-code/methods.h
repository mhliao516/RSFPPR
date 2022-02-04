#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <fstream>
#include <future>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <numeric>
#include <thread>
#include <sys/time.h> 
#include <time.h>
#include <cstring>
#include<iomanip>

#include "graph.h"
#include "random.h"

double get_current_time_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

void generateQueryNode(int num_query_nodes, int n, ofstream& outf) {
    for(int i=0; i<num_query_nodes; i++) {
        int node = drand() % n;
        outf << node << endl;
    }
}

void compare_results(graph& g, int s, string filename) {
    int n = g.n;
    double* gt = new double[n];
    double* power = new double[n];
    double* powerpush = new double[n];
    ifstream gt_f("result/" + filename + "/gt/" + to_string(s) + ".txt");
    ifstream power_f("result/" + filename + "/power/" + to_string(s) + ".txt");
    ifstream powerpush_f("result/" + filename + "/powerpush/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
        gt_f >> gt[i];
        power_f >> power[i];
        powerpush_f >> powerpush[i];
    }
    gt_f.close();
    power_f.close();
    powerpush_f.close();

    double l1_error = 0.;
    double max_error = 0.;

    for(int i=0; i<n; i++) {
        double err = abs(power[i] - gt[i]);
        if(max_error < err) max_error = err;
        l1_error += err;
    }
    cout << "max_error: " << max_error << endl;
    cout << "l1_error: " << l1_error << endl;

    max_error = 0.;
    for(int i=0; i<n; i++) {
        double err = abs(powerpush[i] - gt[i]);
        if(max_error < err) max_error = err;
        l1_error += err;
    }
    cout << "max_error: " << max_error << endl;
    cout << "l1_error: " << l1_error << endl;

    delete[] gt;
    delete[] power;
    delete[] powerpush;
}

void compare_results_prefix(graph& g, int s, double alpha, double eps, string filename, string method) {
    cout << method << endl;
    int n = g.n;
    double method_time;
    double* gt = new double[n];
    double* method_result = new double[n];
    ifstream gt_f("result/" + filename + "/gt/" + to_string(s) + to_string(alpha) + ".txt");
    ifstream method_f("result/" + filename + "/" + method + "/" + to_string(s) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        gt_f >> gt[i];
        method_f >> method_result[i];
    }
    method_f >> method_time;
    gt_f.close();
    method_f.close();
    double l1_error = 0.;
    double max_error = 0.;
    for(int i=0; i<n; i++) {
        double err = abs(method_result[i] - gt[i]);
        if(max_error < err) max_error = err;
        l1_error += err;
    }
    cout << "max_error: " << max_error << endl;
    cout << "l1_error: " << l1_error << endl;
    cout << "time: " << method_time << endl;
    delete[] gt;
    delete[] method_result;
}

void compare_results_prefix_back(graph& g, int s, string filename, string method, double alpha, double eps) {
    cout << method << endl;
    int n = g.n;
    double* gt = new double[n];
    double* method_result = new double[n];
    ifstream gt_f("result/" + filename + "/gtback/" + to_string(s) + to_string(alpha) + ".txt");
    ifstream method_f("result/" + filename + "/" + method + "/" + to_string(s) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        gt_f >> gt[i];
        method_f >> method_result[i];
    }
    gt_f.close();
    method_f.close();
    double l1_error = 0.;
    double max_error = 0.;
    for(int i=0; i<n; i++) {
        double err = abs(method_result[i] - gt[i]);
        if(max_error < err) max_error = err;
        l1_error += err;
    }
    cout << "max_error: " << max_error << endl;
    cout << "l1_error: " << l1_error << endl;
    delete[] gt;
    delete[] method_result;
}

void ground_truth_PageRank(graph& g, int s, double alpha, string filename) {
    cout << s << " degree " << g.degree[s] << " " << g.weightdegree[s] << endl; 
    double t_start = get_current_time_sec();
    int n = g.n;
    double l1_error = 1e-20;
    double* pi = new double[n];
    double* residual = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[s] = 1;
    int number_of_pushes = 0;
    double r_sum = 1;
    const int num_epoch = 8;
    for(int epoch = 0, iter = 0; epoch<num_epoch && r_sum>l1_error; ++epoch) {
        double l1_error_this_epoch = pow(l1_error, (1.0+epoch) / num_epoch);
        const double threshold = l1_error_this_epoch / g.m / 2;
        int num_iter = (int)log(2./l1_error_this_epoch)/log(1./(1.-alpha));
        for(int k=0; k<num_iter && r_sum>l1_error_this_epoch; ++iter, ++k) {
            for(int i=0; i<n; i++) {
                if(residual[i] > threshold*g.weightdegree[i]) {
                    pi[i] += alpha*residual[i];
                    r_sum -= alpha*residual[i];
                    double increment = (1-alpha)*residual[i] / g.weightdegree[i];
                    residual[i] = 0;
                    number_of_pushes += g.degree[i];
                    for(int j=0; j<g.degree[i]; j++) {
                        int neighbor = g.AdjList[i][j];
                        residual[neighbor] += g.edgeWeight[i][j]*increment;
                    }
                }
            }
        }
    }
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/gt/" + to_string(s) + to_string(alpha) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "ground truth time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] residual;
}

void power_push_PageRank(graph& g, int s, double l1_error, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    // double l1_error = 1e-8;
    double* pi = new double[n];
    double* residual = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[s] = 1;
    int number_of_pushes = 0;
    double r_sum = 1;
    const int num_epoch = 8;
    for(int epoch = 0, iter = 0; epoch<num_epoch && r_sum>l1_error; ++epoch) {
        double l1_error_this_epoch = pow(l1_error, (1.0+epoch) / num_epoch);
        const double threshold = l1_error_this_epoch / g.m / 2;
        int num_iter = (int)log(2./l1_error_this_epoch)/log(1./(1.-alpha));
        for(int k=0; k<num_iter && r_sum>l1_error_this_epoch; ++iter, ++k) {
            for(int i=0; i<n; i++) {
                if(residual[i] > threshold*g.degree[i]) {
                    pi[i] += alpha*residual[i];
                    r_sum -= alpha*residual[i];
                    double increment = (1-alpha)*residual[i] / g.weightdegree[i];
                    residual[i] = 0;
                    number_of_pushes += g.degree[i];
                    for(int j=0; j<g.degree[i]; j++) {
                        int neighbor = g.AdjList[i][j];
                        residual[neighbor] += g.edgeWeight[i][j]*increment;
                    }
                }
            }
        }
    }
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/powerpush/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "power push time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] residual;
}

void power_method_PageRank(graph& g, int s, double l1_error, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    // double l1_error = 1e-8;
    double* pi = new double[n];
    double* residual = new double[n];
    double* new_residuals = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[s] = 1;
    int number_of_pushes = 0;
    double r_sum = 1;
    for(int iter = 0; r_sum>l1_error; ++iter) {
        for(int i=0; i<n; i++) {
            number_of_pushes += g.degree[i];
            pi[i] += alpha*residual[i];
            r_sum -= alpha*residual[i];
            double increment = (1-alpha)*residual[i] / g.weightdegree[i];
            residual[i] = 0;
            for(int j=0; j<g.degree[i]; j++) {
                int neighbor = g.AdjList[i][j];
                new_residuals[neighbor] += g.edgeWeight[i][j]*increment;
            }
        }
        for(int i=0; i<n; i++) {
            residual[i] = new_residuals[i];
            new_residuals[i] = 0;
        }
    }
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/power/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "power method time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] residual;
    delete[] new_residuals;
}

void ground_truth_back_PageRank(graph& g, int t, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    int L = 2000;
    double* pi = new double[n];
    double* residual = new double[n];
    double* new_residuals = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[t] = 1;
    int number_of_pushes = 0;
    for(int iter = 0; iter<L; ++iter) {
        for(int i=0; i<n; i++) {
            number_of_pushes += g.degree[i];
            pi[i] += alpha*residual[i];
            double increment = (1-alpha)*residual[i];
            residual[i] = 0;
            for(int j=0; j<g.degree[i]; j++) {
                int neighbor = g.AdjList[i][j];
                new_residuals[neighbor] += g.edgeWeight[i][j]*increment/g.weightdegree[neighbor];
            }
        }
        for(int i=0; i<n; i++) {
            residual[i] = new_residuals[i];
            new_residuals[i] = 0;
        }
    }
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/gtback/" + to_string(t) + to_string(alpha) + ".txt");
    for(int i=0; i<n; i++) {
        fout << setprecision(16) << pi[i] << "\n";
    }
    fout << setprecision(16) << t_end - t_start << "\n";
    cout << "single target ground truth time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] residual;
    delete[] new_residuals;
}

int random_walk(graph& g, int s) {
    double p = prand();
    int index = 0;
    double index_p = 0.;
    int d = g.degree[s];
    while(index_p<p && index<d) {
        index_p += g.edgeWeight[s][index]/g.weightdegree[s];
        index++;
    }
    return g.AdjList[s][index-1];
}

int random_walk_alias(graph& g, int s) {
    int u = drand() % g.degree[s];
    if(prand() < g.alias_probability[s][u]) return g.AdjList[s][u];
    else return g.AdjList[s][g.alias_table[s][u]];
}

void random_walk_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    long num_walks = ceil(n*log(n)/eps/eps);
    cout << "walk num: " << num_walks << endl;
    double* pi = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
    }
    for(int i=0; i<num_walks; i++) {
        int u = s;
        while(prand()>alpha) {
            u = random_walk_alias(g, u);
        }
        pi[u] += 1./num_walks;
    }
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/randomwalk/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "random walk time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
}

int* loop_erased_walk(graph& g, double alpha) {
    int n = g.n;
    bool* intree = new bool[n];
    int* next = new int[n];
    int* root = new int[n];
    for(int i=0; i<n; i++) {
        intree[i] = false;
        next[i] = -1;
        root[i] = -1;
    }
    for(int i=0; i<n; i++) {
        if(g.degree[i] == 0) continue;
        int u = i;
        while(!intree[u]) {
            if(prand()<alpha) {
                intree[u] = true;
                root[u] = u;
            } else {
                // cout << "deg: " << g.degree[u] << endl;
                next[u] = random_walk_alias(g, u);
                u = next[u];
            }
        }
        int r = root[u];
        u = i;
        while(!intree[u]) {
            root[u] = r;
            intree[u] = true;
            u = next[u];
        }
    }

    return root;
}

void fora_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    long num_walks = ceil(n*log(n)/eps/eps);
    cout << "walk num: " << num_walks << endl;
    double* pi = new double[n];
    double* residual = new double[n];
    bool* is_in_queue = new bool[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
        is_in_queue[i] = false;
    }
    queue<int> q;
    residual[s] = 1;
    q.push(s);
    is_in_queue[s] = true;
    double push_threshold = eps/n;
    double r_sum = 1.0;

    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        r_sum -= alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[u];
            if(!is_in_queue[update] && (residual[update]/g.weightdegree[u] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    cout << "r_sum: " << r_sum << endl;
    long rw_count = 0;
    for(int id=0; id<n; id++) {
        if(residual[id]>0) {
            int num_id_rw = ceil(residual[id]*num_walks);
            rw_count += num_id_rw;
            double ppr_incre = residual[id]/(double)num_id_rw;
            for(int i=0; i<num_id_rw; i++) {
                int u = id;
                while(prand()>alpha) {
                    u = random_walk_alias(g, u);
                }
                pi[u] += ppr_incre;
            }
        }
    }
    cout << "rw_count: " << rw_count << endl;
 
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/fora/" + to_string(s) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "fora time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
}

void foral_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    long num_walks = ceil(n*log(n)/eps/eps);
    cout << "walk num: " << num_walks << endl;
    double* pi = new double[n];
    double* residual = new double[n];
    bool* is_in_queue = new bool[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
        is_in_queue[i] = false;
    }
    queue<int> q;
    residual[s] = 1;
    q.push(s);
    is_in_queue[s] = true;
    double push_threshold = eps/n;
    double r_sum = 1.0;

    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        r_sum -= alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[u];
            if(!is_in_queue[update] && (residual[update] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    cout << "r_sum: " << r_sum << endl;
    int rw_count = 0;
    int num_loop_erased_rw = push_threshold*num_walks;
    for(int i=0; i<num_loop_erased_rw; i++) {
        int* root = loop_erased_walk(g, alpha);
        rw_count++;
        for(int id=0; id<n; id++) {
            if(g.degree[id] == 0) continue;
            pi[root[id]] += residual[id]/(double)num_loop_erased_rw; 
        }
        delete[] root;
    }
    cout << "rw_count: " << rw_count << endl;
 
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/foral/" + to_string(s)+ to_string(alpha) + to_string(eps)  + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "foral time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
}

void foralv_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    long num_walks = ceil(n*log(n)/eps/eps);
    cout << "walk num: " << num_walks << endl;
    double* pi = new double[n];
    double* residual = new double[n];
    bool* is_in_queue = new bool[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
        is_in_queue[i] = false;
    }
    queue<int> q;
    residual[s] = 1;
    q.push(s);
    is_in_queue[s] = true;
    double push_threshold = eps/n;
    double r_sum = 1.0;

    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        r_sum -= alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[u];
            if(!is_in_queue[update] && (residual[update] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    cout << "r_sum: " << r_sum << endl;
    int rw_count = 0;
    int num_loop_erased_rw = push_threshold*num_walks;
    for(int i=0; i<num_loop_erased_rw; i++) {
        int* root = loop_erased_walk(g, alpha);
        rw_count++;
        double* partition_acc = new double[n]();
        double* residual_over_partition = new double[n]();
        for(int id=0; id<n; id++) {
            if(g.degree[id] == 0) continue;
            partition_acc[root[id]] += g.weightdegree[id]; 
            residual_over_partition[root[id]] += residual[id];
        }
        for(int id=0; id<n; id++) {
            if(g.degree[id] == 0) continue;
            pi[id] += g.weightdegree[id]*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_loop_erased_rw;
        }
        delete[] root;
    }
    cout << "rw_count: " << rw_count << endl;
 
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/foralv/" + to_string(s)+ to_string(alpha) + to_string(eps)  + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "foralv time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
}

void speedppr_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    double num_walks = n*log(n)/eps/eps;
    double push_threshold = 1./num_walks;
    double* pi = new double[n];
    double* residual = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[s] = 1;
    queue<int> q;
    bool* is_in_queue = new bool[n];
    q.push(s);
    is_in_queue[s] = true;

    double queue_threshold = (double)n*n/10./g.m;
    double initial_size = std::max(num_walks/1000./log(n), 1.0);
    double step_size = std::max(powf(initial_size, 1.0/3.0), 2.0f);
    cout << "queue threshold: " << queue_threshold << endl;
    cout << "initial size: " << initial_size << endl;
    cout << "step_size: " << step_size << endl;

    for(int scale_factor = initial_size; scale_factor >= 1 && q.size() < queue_threshold;) {
        while(!q.empty() && q.size() < queue_threshold) {
            int u = q.front();
            q.pop();
            is_in_queue[u] = false;
            if(residual[u] > push_threshold*g.weightdegree[u]*scale_factor) {
                pi[u] += alpha*residual[u];
                double increment = (1-alpha)*residual[u]/g.weightdegree[u];
                residual[u] = 0;
                for(int i=0; i<g.degree[u]; i++) {
                    int update = g.AdjList[u][i];
                    residual[update] += increment*g.edgeWeight[u][i];
                    if(!is_in_queue[update]) {
                        q.push(update);
                        is_in_queue[update] = true;
                    }
                }
            }
        }
        scale_factor /= step_size;
        cout << "scale factor = " << scale_factor << " q size: " << q.size() << endl;
        if(q.empty()) {
            for(int id=0; id<n; id++) {
                if(residual[id] > g.weightdegree[id]*push_threshold*scale_factor) {
                    q.push(id);
                    is_in_queue[id] = true;
                }
            }
        }
    }

    int queue_size = q.size();
    cout << "q_size: " << queue_size << endl;
    for(;queue_size > queue_threshold;) {
        queue_size = 0;
        for(int id=0; id<n; id++) {
            if(residual[id] > g.weightdegree[id]*push_threshold) {
                pi[id] += alpha*residual[id];
                double increment = (1-alpha)*residual[id] / g.weightdegree[id];
                residual[id] = 0;
                queue_size += g.degree[id];
                for(int j=0; j<g.degree[id]; j++) {
                    int neighbor = g.AdjList[id][j];
                    residual[neighbor] += g.edgeWeight[id][j]*increment;
                }
            }
        }
    }

    queue_size = 0;
    queue<int> empty;
    swap(q, empty);
    for(int i=0; i<n; i++) is_in_queue[i] = false;
    for(int i=0; i<n; i++) {
        if(residual[i]/g.weightdegree[i] > push_threshold) {
            q.push(i);
            is_in_queue[i] = true;
            queue_size++;
        }
    }
    cout << "queue size: " << queue_size << endl;
    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[u];
            if(!is_in_queue[update] && (residual[update]/g.weightdegree[u] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    long rw_count = 0;
    for(int id=0; id<n; id++) {
        if(residual[id]>0) {
            int num_id_rw = ceil(residual[id]*num_walks);
            rw_count += num_id_rw;
            double ppr_incre = residual[id]/(double)num_id_rw;
            for(int i=0; i<num_id_rw; i++) {
                int u = id;
                while(prand()>alpha) {
                    u = random_walk_alias(g, u);
                }
                pi[u] += ppr_incre;
            }
        }
    }
    cout << "rw_count: " << rw_count << endl;
 
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/speedppr/" + to_string(s)+ to_string(alpha) + to_string(eps)  + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "speedppr time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
}

void speedl_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    double num_walks = n*log(n)/eps/eps;
    double push_threshold = 1./num_walks*log(n);
    double* pi = new double[n];
    double* residual = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[s] = 1;
    queue<int> q;
    bool* is_in_queue = new bool[n];
    q.push(s);
    is_in_queue[s] = true;

    double queue_threshold = (double)n*n/10./g.m;
    double initial_size = std::max(num_walks/1000./log(n), 1.0);
    double step_size = std::max(powf(initial_size, 1.0/3.0), 2.0f);
    cout << "queue threshold: " << queue_threshold << endl;
    cout << "initial size: " << initial_size << endl;
    cout << "step_size: " << step_size << endl;

    for(int scale_factor = initial_size; scale_factor >= 1 && q.size() < queue_threshold;) {
        while(!q.empty() && q.size() < queue_threshold) {
            int u = q.front();
            q.pop();
            is_in_queue[u] = false;
            if(residual[u] > push_threshold*scale_factor) {
                pi[u] += alpha*residual[u];
                double increment = (1-alpha)*residual[u]/g.weightdegree[u];
                residual[u] = 0;
                for(int i=0; i<g.degree[u]; i++) {
                    int update = g.AdjList[u][i];
                    residual[update] += increment*g.edgeWeight[u][i];
                    if(!is_in_queue[update]) {
                        q.push(update);
                        is_in_queue[update] = true;
                    }
                }
            }
        }
        scale_factor /= step_size;
        cout << "scale factor = " << scale_factor << " q size: " << q.size() << endl;
        if(q.empty()) {
            for(int id=0; id<n; id++) {
                if(residual[id] > push_threshold*scale_factor) {
                    q.push(id);
                    is_in_queue[id] = true;
                }
            }
        }
    }

    int queue_size = q.size();
    cout << "q_size: " << queue_size << endl;
    for(;queue_size > queue_threshold;) {
        queue_size = 0;
        for(int id=0; id<n; id++) {
            if(residual[id] > push_threshold) {
                pi[id] += alpha*residual[id];
                double increment = (1-alpha)*residual[id] / g.weightdegree[id];
                residual[id] = 0;
                queue_size += g.degree[id];
                // queue_size++;
                for(int j=0; j<g.degree[id]; j++) {
                    int neighbor = g.AdjList[id][j];
                    residual[neighbor] += g.edgeWeight[id][j]*increment;
                }
            }
        }
    }

    queue_size = 0;
    queue<int> empty;
    swap(q, empty);
    for(int i=0; i<n; i++) is_in_queue[i] = false;
    for(int i=0; i<n; i++) {
        if(residual[i] > push_threshold) {
            q.push(i);
            is_in_queue[i] = true;
            queue_size++;
            // cout << i << " " << residual[i] << " " << g.degree[i] << " " << g.weightdegree[i] << endl;
        }
    }
    cout << "queue size: " << queue_size << endl;
    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[u];
            if(!is_in_queue[update] && (residual[update] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    int rw_count = 0;
    int num_loop_erased_rw = push_threshold*num_walks;
    cout << num_loop_erased_rw << endl;
    for(int i=0; i<num_loop_erased_rw; i++) {
        // out << i << endl;
        int* root = loop_erased_walk(g, alpha);
        // cout << i << endl;
        rw_count++;
        for(int id=0; id<n; id++) {
            if(g.degree[id]==0) continue; 
            pi[root[id]] += residual[id]/(double)num_loop_erased_rw; 
        }
        delete[] root;
    }
    cout << "rw_count: " << rw_count << endl;

    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/speedl/" + to_string(s) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "speedl time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] residual;
}

void speedlv_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    // cout << "node: " << s << " degree: " << g.degree[s] << endl;
    double t_start = get_current_time_sec();
    int n = g.n;
    double num_walks = n*log(n)/eps/eps;
    double push_threshold = 1./num_walks*log(n);
    double* pi = new double[n];
    double* residual = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[s] = 1;
    queue<int> q;
    bool* is_in_queue = new bool[n];
    q.push(s);
    is_in_queue[s] = true;

    double queue_threshold = (double)n*n/10./g.m;
    double initial_size = std::max(num_walks/1000./log(n), 1.0);
    double step_size = std::max(powf(initial_size, 1.0/3.0), 2.0f);
    cout << "queue threshold: " << queue_threshold << endl;
    cout << "initial size: " << initial_size << endl;
    cout << "step_size: " << step_size << endl;

    for(int scale_factor = initial_size; scale_factor >= 1 && q.size() < queue_threshold;) {
        while(!q.empty() && q.size() < queue_threshold) {
            int u = q.front();
            q.pop();
            is_in_queue[u] = false;
            if(residual[u] > push_threshold*scale_factor) {
                pi[u] += alpha*residual[u];
                double increment = (1-alpha)*residual[u]/g.weightdegree[u];
                residual[u] = 0;
                for(int i=0; i<g.degree[u]; i++) {
                    int update = g.AdjList[u][i];
                    residual[update] += increment*g.edgeWeight[u][i];
                    if(!is_in_queue[update]) {
                        q.push(update);
                        is_in_queue[update] = true;
                    }
                }
            }
        }
        scale_factor /= step_size;
        cout << "scale factor = " << scale_factor << " q size: " << q.size() << endl;
        if(q.empty()) {
            for(int id=0; id<n; id++) {
                if(residual[id] > push_threshold*scale_factor) {
                    q.push(id);
                    is_in_queue[id] = true;
                }
            }
        }
    }

    int queue_size = q.size();
    cout << "q_size: " << queue_size << endl;
    for(;queue_size > queue_threshold;) {
        queue_size = 0;
        for(int id=0; id<n; id++) {
            if(residual[id] > push_threshold) {
                pi[id] += alpha*residual[id];
                double increment = (1-alpha)*residual[id] / g.weightdegree[id];
                residual[id] = 0;
                queue_size += g.degree[id];
                for(int j=0; j<g.degree[id]; j++) {
                    int neighbor = g.AdjList[id][j];
                    residual[neighbor] += g.edgeWeight[id][j]*increment;
                }
            }
        }
    }

    queue_size = 0;
    queue<int> empty;
    swap(q, empty);
    for(int i=0; i<n; i++) is_in_queue[i] = false;
    for(int i=0; i<n; i++) {
        if(residual[i] > push_threshold) {
            q.push(i);
            is_in_queue[i] = true;
            queue_size++;
        }
    }
    cout << "queue size: " << queue_size << endl;
    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[u];
            if(!is_in_queue[update] && (residual[update] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    int rw_count = 0;
    int num_loop_erased_rw = push_threshold*num_walks;
    for(int i=0; i<num_loop_erased_rw; i++) {
        int* root = loop_erased_walk(g, alpha);
        rw_count++;
        double* partition_acc = new double[n]();
        double* residual_over_partition = new double[n]();
        for(int id=0; id<n; id++) {
            if(g.degree[id] == 0) continue;
            partition_acc[root[id]] += g.weightdegree[id]; 
            residual_over_partition[root[id]] += residual[id];
        }
        for(int id=0; id<n; id++) {
            if(g.degree[id] == 0) continue;
            pi[id] += g.weightdegree[id]*residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_loop_erased_rw;
        }
        delete[] root;
    }
    cout << "rw_count: " << rw_count << endl;

    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/speedlv/" + to_string(s) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "speedlv time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] residual;
}

void powerpush_simple_PageRank(graph& g, int s, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    double num_walks = n*log(n)/eps/eps;
    double push_threshold = 1./num_walks;
    eps = 1e-3;
    double* pi = new double[n];
    double* residual = new double[n];
    double* new_residuals = new double[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
    }
    residual[s] = 1;
    int number_of_pushes = 0;
    double r_sum = 1;
    int L = (int)ceil(log(eps)/log(1-alpha))+1;
    L = 10;
    for(int iter = 0; iter<L; ++iter) {
        for(int i=0; i<n; i++) {
            if(residual[i]>0) {
                number_of_pushes += g.degree[i];
                pi[i] += alpha*residual[i];
                r_sum -= alpha*residual[i];
                double increment = (1-alpha)*residual[i] / g.weightdegree[i];
                residual[i] = 0;
                for(int j=0; j<g.degree[i]; j++) {
                    int neighbor = g.AdjList[i][j];
                    new_residuals[neighbor] += g.edgeWeight[i][j]*increment;
                }
            }
        }
        for(int i=0; i<n; i++) {
            residual[i] = new_residuals[i];
            new_residuals[i] = 0;
        }
    }

    queue<int> q;
    bool* is_in_queue = new bool[n];
    int queue_size = 0;
    for(int i=0; i<n; i++) is_in_queue[i] = false;
    for(int i=0; i<n; i++) {
        if(residual[i]/g.weightdegree[i] > push_threshold) {
            q.push(i);
            is_in_queue[i] = true;
            queue_size++;
        }
    }
    cout << "queue size: " << queue_size << endl;
    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        r_sum -= alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[u];
            if(!is_in_queue[update] && (residual[update]/g.weightdegree[u] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/power/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "power method time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] residual;
    delete[] new_residuals;
} 

void back_PageRank(graph& g, int t, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    // int num_walks = ceil(n*log(n)/eps/eps);
    // cout << "walk num: " << num_walks << endl;
    double* pi = new double[n];
    double* residual = new double[n];
    bool* is_in_queue = new bool[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
        is_in_queue[i] = false;
    }
    queue<int> q;
    residual[t] = 1;
    q.push(t);
    is_in_queue[t] = true;
    double push_threshold = eps/n;

    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[update];
            if(!is_in_queue[update] && (residual[update] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/back/" + to_string(t) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << setprecision(16) << t_end - t_start << "\n";
    cout << "back time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
}

void rback_PageRank(graph& g, int t, double eps, double alpha, string filename) {
    int type = 1;
    double t_start = get_current_time_sec();
    int n = g.n;
    double* pi = new double[n];
    int pi_count = 0;
    int* H[2];
    int* U[2];
	int *candidate_set[2];
	int candidate_count[2];
	double *residue[2];
    candidate_count[0]=0;
	candidate_count[1]=0;
	H[0] = new int[n];
	H[1] = new int[n];
	U[0] = new int[n];
	U[1] = new int[n];
	candidate_set[0] = new int[n];
	candidate_set[1] = new int[n];
	residue[0] = new double[n];
	residue[1] = new double[n];
    for(int i=0; i<n; i++) {
        residue[0][i] = 0;
        residue[1][i] = 0;
        pi[i] = 0;
        H[0][i] = 0;
		H[1][i] = 0;
		U[0][i] = 0;
		U[1][i] = 0;
		candidate_set[0][i]=0;
		candidate_set[1][i]=0;
    }

    int templevel = 0;
    residue[0][t] = 1;
    candidate_set[0][0] = t;
    candidate_count[0] = 1;
    int L = (int)ceil(log(eps)/log(1-alpha))+1;

    while(templevel <= L) {
        int templevelID = templevel%2;
        int newlevelID = (templevel+1)%2;
        int candidateCnt = candidate_count[templevelID];
        if(candidateCnt == 0) break;
        candidate_count[templevelID] = 0;
        for(int j=0; j<candidateCnt; j++) {
            int tempNode = candidate_set[templevelID][j];
            double tempR = residue[templevelID][tempNode];
            U[templevelID][tempNode] = 0;
            residue[templevelID][tempNode] = 0;
            if(H[1][tempNode] == 0) {
                H[0][pi_count++] = tempNode;
                H[1][tempNode] = 1;
            }
            pi[tempNode] += alpha*tempR;

            if(templevel == L) continue;

            double ran = prand();
            for(int k=0; k<g.degree[tempNode]; k++) {
                int update = g.AdjList[tempNode][k];
                double update_deg = g.weightdegree[update];
                double incre = tempR*(1-alpha)*g.edgeWeight[tempNode][k]/update_deg;

                if(type==1) { // additive error
                    if(sqrt(update_deg)*eps <= (1-alpha)*tempR) {
                        residue[newlevelID][update] += incre;
                    } else {
                        if(sqrt(update_deg)*ran*eps <= (1-alpha)*tempR) {
                            double ran_incre = eps/sqrt(update_deg);
                            residue[newlevelID][update] += ran_incre;
                        } else {
                            break;
                        }
                    }
                } else { // relative error
                    if(update_deg*eps <= (1-alpha)*tempR) {
                        residue[newlevelID][update] += incre;
                    } else {
                        if(update_deg*eps*ran <= (1-alpha)*tempR) {
                            double ran_incre = eps;
                            residue[newlevelID][update] += ran_incre;
                        } else {
                            break;
                        }
                    }
                }

                if((U[newlevelID][update] == 0) && (residue[newlevelID][update]>eps)) {
                    U[newlevelID][update] = 1;
                    candidate_set[newlevelID][candidate_count[newlevelID]++] = update;
                }
            }
        }
        templevel++;
    }

    double t_end = get_current_time_sec();

    ofstream fout("result/" + filename + "/rback/" + to_string(t) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << setprecision(16) << t_end - t_start << "\n";
    cout << "rback time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
    delete[] H[0];
	delete[] H[1];
	delete[] U[0];
	delete[] U[1];
	delete[] residue[0];
	delete[] residue[1];
	delete[] candidate_set[0];
	delete[] candidate_set[1];
}

void backlv_PageRank(graph& g, int t, double eps, double alpha, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    double num_walks = n*log(n)/eps/eps;
    double* pi = new double[n];
    double* residual = new double[n];
    bool* is_in_queue = new bool[n];
    for(int i=0; i<n; i++) {
        pi[i] = 0;
        residual[i] = 0;
        is_in_queue[i] = false;
    }
    queue<int> q;
    residual[t] = 1;
    q.push(t);
    is_in_queue[t] = true;
    double scale_factor = 100.;
    double push_threshold = scale_factor*eps/n;

    while(!q.empty()) {
        int u = q.front();
        q.pop();
        is_in_queue[u] = false;
        pi[u] += alpha*residual[u];
        double increment = (1-alpha)*residual[u];
        residual[u] = 0;
        for(int w=0; w<g.degree[u]; w++) {
            int update = g.AdjList[u][w];
            residual[update] += increment*g.edgeWeight[u][w]/g.weightdegree[update];
            if(!is_in_queue[update] && (residual[update] > push_threshold)) {
                is_in_queue[update] = true;
                q.push(update);
            }
        }
    }

    int num_loop_erased_rw = ceil(push_threshold*num_walks*alpha);
    cout << "num of spanning forests: " << num_loop_erased_rw << endl;

    for(int i=0; i<num_loop_erased_rw; i++) {
        // cout << i << endl;
        int* root = loop_erased_walk(g, alpha);
        double* partition_acc = new double[n]();
        double* residual_over_partition = new double[n]();
        for(int id=0; id<n; id++) {
            if(g.degree[id] == 0) continue;
            partition_acc[root[id]] += g.weightdegree[id];
            residual_over_partition[root[id]] += residual[id]*g.weightdegree[id];
        }
        for(int id=0; id<n; id++) {
            if(g.degree[id] == 0) continue;
            pi[id] += residual_over_partition[root[id]]/partition_acc[root[id]]/(double)num_loop_erased_rw;
        }
        delete[] root;
        delete[] residual_over_partition;
        delete[] partition_acc;
    }

    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/backlv/" + to_string(t) + to_string(alpha) + to_string(eps) + ".txt");
    for(int i=0; i<n; i++) {
        fout << setprecision(16) << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "backlv time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
}

void compute_uniform_distribution(graph& g, int s, string filename) {
    double t_start = get_current_time_sec();
    int n = g.n;
    double* pi = new double[n]();
    for(int i=0; i<n; i++) {
        pi[i] = g.weightdegree[i]/2/g.m;
    }
 
    double t_end = get_current_time_sec();
    ofstream fout("result/" + filename + "/uniform/" + to_string(s) + ".txt");
    for(int i=0; i<n; i++) {
        fout << pi[i] << "\n";
    }
    fout << t_end - t_start << "\n";
    cout << "uniform time: " << t_end - t_start << "\n";
    fout.close();
    delete[] pi;
} 

void compare_results_prefix_plot(graph& g, int s, double alpha, string filename, string method, double time_result[][7], double error_result[][7], int &method_pointer) {
    cout << method << endl;
    int n = g.n;
    double method_time;
    double* gt = new double[n];
    double* method_result = new double[n];
    for(int i=0; i<5; i++) {
        double eps = time_result[i][0];
        ifstream gt_f("result/" + filename + "/gt/" + to_string(s) + to_string(alpha) + ".txt");
        ifstream method_f("result/" + filename + "/" + method + "/" + to_string(s) + to_string(alpha) + to_string(eps) + ".txt");
        for(int i=0; i<n; i++) {
            gt_f >> gt[i];
            method_f >> method_result[i];
        }
        method_f >> method_time;
        gt_f.close();
        method_f.close();
        double l1_error = 0.;
        double max_error = 0.;
        for(int i=0; i<n; i++) {
            double err = abs(method_result[i] - gt[i]);
            if(max_error < err) max_error = err;
            l1_error += err;
        }
        cout << "max_error: " << max_error << endl;
        cout << "l1_error: " << l1_error << endl;
        cout << "time: " << method_time << endl;
        time_result[i][method_pointer] += method_time;
        error_result[i][method_pointer] += l1_error;
    }
    
    delete[] gt;
    delete[] method_result;

    method_pointer++;
}

void plot_results(graph& g, int num_query_nodes, double alpha, string filename) {
    double time_result[5][7];
    double error_result[5][7];
    for(int i=0; i<5; i++) {
        for(int j=0; j<7; j++) {
            time_result[i][j] = 0;
            error_result[i][j] = 0;
        }
    }
    time_result[0][0] = 0.5;
    time_result[1][0] = 0.4;
    time_result[2][0] = 0.3;
    time_result[3][0] = 0.2;
    time_result[4][0] = 0.1;
    error_result[0][0] = 0.5;
    error_result[1][0] = 0.4;
    error_result[2][0] = 0.3;
    error_result[3][0] = 0.2;
    error_result[4][0] = 0.1;

    int method_pointer;

    ifstream infile("data/" + filename + ".query");
    for(int i=0; i<num_query_nodes; i++) {
        int s;
        infile >> s;
        method_pointer = 1;
        compare_results_prefix_plot(g, s, alpha, filename, "fora", time_result, error_result, method_pointer); // 2
        compare_results_prefix_plot(g, s, alpha, filename, "foral", time_result, error_result, method_pointer); // 3
        compare_results_prefix_plot(g, s, alpha, filename, "foralv", time_result, error_result, method_pointer); // 4
        compare_results_prefix_plot(g, s, alpha, filename, "speedppr", time_result, error_result, method_pointer); // 5
        compare_results_prefix_plot(g, s, alpha, filename, "speedl", time_result, error_result, method_pointer); // 6
        compare_results_prefix_plot(g, s, alpha, filename, "speedlv", time_result, error_result, method_pointer); // 7
    }
    infile.close();
    for(int i=0; i<5; i++) {
        for(int j=1; j<7; j++) {
            time_result[i][j] /= (double)num_query_nodes;
            error_result[i][j] /= (double)num_query_nodes;
        }
    }

    ofstream fout1("plot_result/new" + filename + to_string(alpha) + "time.txt");
    for(int i=0; i<5; i++) {
        for(int j=0; j<7; j++) {
            fout1 << time_result[i][j] << " ";
        }
        fout1 << endl;
    }
    fout1.close();

    ofstream fout2("plot_result/new" + filename + to_string(alpha) + "error.txt");
    for(int i=0; i<5; i++) {
        for(int j=0; j<7; j++) {
            fout2 << error_result[i][j] << " ";
        }
        fout2 << endl;
    }
    fout2.close();
}
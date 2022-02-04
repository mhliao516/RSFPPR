#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>    
#include <cstdlib>
#include <cstring>
#include <fstream>

#include "methods.h"

int main(int argc, char *argv[]) {

    SFMTinitialize();

    string filename;
    string algo = "RSFPPR";

    int num_query_nodes = 20;
    double alpha = 0.2;
    double l1_error = 1e-12;
    double eps = 0.1;
    double reps = 1e-12;

    int i=1;
    char *endptr;
    while(i<argc) {
        if(!strcmp(argv[i], "-d")) {
            filename = argv[++i];
        } else if(!strcmp(argv[i], "-algo")) {
            algo = argv[++i];
        } else if(!strcmp(argv[i], "-n")) {
            num_query_nodes = strtod(argv[++i], &endptr);
        } else if(!strcmp(argv[i], "-alpha")) {
            alpha = strtod(argv[++i], &endptr);
        } else if(!strcmp(argv[i], "-e")) {
            eps = strtod(argv[++i], &endptr);
        } else if(!strcmp(argv[i], "-re")) {
            reps = strtod(argv[++i], &endptr);
        }
        i++;
    }

    graph g;
    // g.read_graph_weighted(filename);
    g.read_graph(filename);

    if(algo == "GEN_QUERY") {
        ifstream infile("data/" + filename + ".query");
        if(infile.good()) {
            cerr << "query file already generated!" << endl;
            return 0;
        }
        ofstream outf("data/" + filename + ".query");
        generateQueryNode(num_query_nodes, g.n, outf);
        outf.close();
    } else if(algo == "GROUND_TRUTH") {
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s;
            infile >> s;
            ground_truth_PageRank(g, s, alpha, filename);
        }
        infile.close();
    } else if(algo == "GROUND_TRUTH_BACK") {
        g.compute_degree_order();
        int* order = g.get_degree_order();
        for(int i=0; i<num_query_nodes; i++) {
            int t = order[i];
            // cout << t << " " << g.degree[t] << endl;
            ground_truth_back_PageRank(g, t, alpha, filename);
        }
    } else if(algo == "ITERATIVE_METHOD") {
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s;
            infile >> s;
            power_method_PageRank(g, s, l1_error, alpha, filename);
            power_push_PageRank(g, s, l1_error, alpha, filename);
        }
        infile.close();
    } else if(algo == "RANDOM_WALK") {
        g.alias_preprocess();
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s;
            infile >> s;
            random_walk_PageRank(g, s, eps, alpha, filename);
        }
        infile.close();
    } else if(algo == "FORA") {
        g.alias_preprocess();
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s;
            infile >> s;
            fora_PageRank(g, s, eps, alpha, filename);
            foral_PageRank(g, s, eps, alpha, filename);
            foralv_PageRank(g, s, eps, alpha, filename);
        }
        infile.close();
    } else if(algo == "SPEEDPPR") {
        g.alias_preprocess();
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s;
            infile >> s;
            speedppr_PageRank(g, s, eps, alpha, filename);
            speedl_PageRank(g, s, eps, alpha, filename);
            speedlv_PageRank(g, s, eps, alpha, filename);
        }
        infile.close();
    } else if(algo == "BACK") {
        g.compute_degree_order();
        int* order = g.get_degree_order();
        g.alias_preprocess();
        for(int i=0; i<num_query_nodes; i++) {
            int t = order[i];
            back_PageRank(g, t, eps, alpha, filename);
            backlv_PageRank(g, t, eps, alpha, filename);
        }
    } else if(algo == "RBACK") {
        g.compute_degree_order();
        g.sortAdjList();
        g.alias_preprocess();
        int* order = g.get_degree_order();
        for(int i=0; i<num_query_nodes; i++) {
            int t = order[i];
            cout << "t: " << t <<  " " << g.degree[t] << endl;
            rback_PageRank(g, t, reps, alpha, filename);
        }
    } else if(algo == "COMPARE_RESULTS_BACK") {
        g.compute_degree_order();
        int* order = g.get_degree_order();
        for(int i=0; i<num_query_nodes; i++) {
            int t = order[i];
            compare_results_prefix_back(g, t, filename, "back", alpha, eps);
            compare_results_prefix_back(g, t, filename, "rback", alpha, reps);
            compare_results_prefix_back(g, t, filename, "backlv", alpha, eps);
        }
    } else if(algo == "COMPARE_RESULTS") {
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s;
            infile >> s;
            compare_results_prefix(g, s, alpha, eps, filename, "fora");
            compare_results_prefix(g, s, alpha, eps, filename, "foral");
            compare_results_prefix(g, s, alpha, eps, filename, "foralv");
            compare_results_prefix(g, s, alpha, eps, filename, "speedppr");
            compare_results_prefix(g, s, alpha, eps, filename, "speedl");
            compare_results_prefix(g, s, alpha, eps, filename, "speedlv");
        }
        infile.close();
    } else if(algo == "VERY_SMALL_ALPHA") {
        g.alias_preprocess();
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s;
            infile >> s;
            speedlv_PageRank(g, s, eps, alpha, filename);
            compute_uniform_distribution(g, s, filename);
        }
        infile.close();
    } else if(algo == "PLOT_RESULTS") {
        plot_results(g, num_query_nodes, alpha, filename);
    }

    return 0;
}

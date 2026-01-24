#ifndef _KDEFECTIVE_CLIQUE_BB_BITSET_
#define _KDEFECTIVE_CLIQUE_BB_BITSET_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"
#include "Dinic_dijkstra.h"
#include <bits/extc++.h>

class kDefectiveClique_BB {
private:
	ui n;
	ui CNT;
	ui *pstart;
	ui *pstart_R;
	ui *pend;
	ui *edges;
    ui *removed_level;

    ui *degree;
    ui *degree_in_S;

    ui K, UB, max_n;
    ui *best_solution;
    ui best_solution_size;

    ui *SR; // union of S and R, where S is at the front
    ui *SR_rid; // reverse ID for SR
    std::queue<ui> Qv;
	std::vector<ui> zero_degree_vertices;

    ui *neighbors;
    ui *nonneighbors;
    ui *buf;
    ui *buf1;
    ui *buf2;
	ui *buf3;
	ui *buf4;
	ui *buf5;
	ui *buf6;
	ui *buf7;
	ui *buf8;
	ui *buf9;
	ui *buf10;
	ui *buf11;
	ui history_max_n;
	ui two_stage;
	ui flow_prune;
	ui look_ahead_prune;
	ui color_prune;
	ui branch_tot;

	ListLinearHeap *calc_heap;

	DinicSolver *solver;

	ui *nonneighbor_1;
	ui *nonneighbor_2;
    char *vis;

	bool *matrix;

	std::set<ui> my_vec;

	const double eps = 1e-6;
	__gnu_pbds::gp_hash_table <long long,ui> color_table;
	__gnu_pbds::gp_hash_table <long long,ui> color_cs;

public:
    kDefectiveClique_BB() {
    	history_max_n = n = 0;
		CNT = 0;
		two_stage = 0;
		flow_prune = 0;
		look_ahead_prune = 0;
    	pstart = pstart_R = NULL;
		branch_tot = 0;
    	pend = NULL;
    	edges = NULL;
        degree = degree_in_S = NULL;
        
        best_solution = NULL;
        K = best_solution_size = UB = 0;

        SR = SR_rid = NULL;
        removed_level = NULL;

        neighbors = nonneighbors = NULL;
        buf = buf1 = buf2 = buf3 = buf4 = buf5 = buf6 = buf7 = buf8 = buf9 = buf10 = buf11 = NULL;
		nonneighbor_1 = nonneighbor_2 = NULL;
        vis = NULL;
		calc_heap = NULL;
		solver = NULL;
		matrix = NULL;
    }

    ~kDefectiveClique_BB() {
    	if(pstart != NULL) {
    		delete[] pstart;
    		pstart = NULL;
    	}
    	if(pstart_R != NULL) {
    		delete[] pstart_R;
    		pstart_R = NULL;
    	}
    	if(pend != NULL) {
    		delete[] pend;
    		pend = NULL;
    	}
    	if(edges != NULL) {
    		delete[] edges;
    		edges = NULL;
    	}
        if(degree != NULL) {
            delete[] degree;
            degree = NULL;
        }
        if(degree_in_S != NULL) {
            delete[] degree_in_S;
            degree_in_S = NULL;
        }
        if(best_solution != NULL) {
            delete[] best_solution;
            best_solution = NULL;
        }
        if(SR != NULL) {
            delete[] SR;
            SR = NULL;
        }
        if(SR_rid != NULL) {
            delete[] SR_rid;
            SR_rid = NULL;
        }
        if(removed_level != NULL) {
        	delete[] removed_level;
        	removed_level = NULL;
        }
        if(neighbors != NULL) {
        	delete[] neighbors;
        	neighbors = NULL;
        }
        if(nonneighbors != NULL) {
        	delete[] nonneighbors;
        	nonneighbors = NULL;
        }
        if(buf != NULL) {
        	delete[] buf;
        	buf = NULL;
        }
        if(buf1 != NULL) {
        	delete[] buf1;
        	buf1 = NULL;
        }
        if(buf2 != NULL) {
        	delete[] buf2;
        	buf2 = NULL;
        }
        if(buf3 != NULL) {
        	delete[] buf3;
        	buf3 = NULL;
        }
        if(buf4 != NULL) {
        	delete[] buf4;
        	buf4 = NULL;
        }
        if(buf5 != NULL) {
        	delete[] buf5;
        	buf5 = NULL;
        }
        if(buf6 != NULL) {
        	delete[] buf6;
        	buf6 = NULL;
        }
        if(buf7 != NULL) {
        	delete[] buf7;
        	buf7 = NULL;
        }
        if(buf8 != NULL) {
        	delete[] buf8;
        	buf8 = NULL;
        }
        if(buf9 != NULL) {
        	delete[] buf9;
        	buf9 = NULL;
        }
        if(buf10 != NULL) {
        	delete[] buf10;
        	buf10 = NULL;
        }
        if(buf11 != NULL) {
        	delete[] buf11;
        	buf11 = NULL;
        }
        if(vis != NULL) {
        	delete[] vis;
        	vis = NULL;
        }
		if(nonneighbor_1 != NULL){
			delete[] nonneighbor_1;
			nonneighbor_1 = NULL;
		}
		if(nonneighbor_2 != NULL){
			delete[] nonneighbor_2;
			nonneighbor_2 = NULL;
		}
		if(calc_heap != NULL){
			delete calc_heap;
			calc_heap = NULL;
		}
		if(solver != NULL){
			delete solver;
			solver = NULL;
		}
		if(matrix != NULL){
			delete[] matrix;
			matrix = NULL;
		}
    }

	void allocate_memory(ui _max_n,ui _K){
		max_n = _max_n;K = _K;
		
		degree = new ui[max_n+1];
		degree_in_S = new ui[max_n+1];
		best_solution = new ui[max_n+1];
		SR = new ui[max_n+1];
		SR_rid = new ui[max_n+1];
		removed_level = new ui[max_n+1];
		pstart_R = new ui[max_n+1];
		pend = new ui[max_n+1];
		if(max_n > K+1) {
			neighbors = new ui[max_n+1];
			nonneighbors = new ui[max_n+1];
			nonneighbor_1 = new ui[max_n+1];
			nonneighbor_2 = new ui[max_n+1];
			buf = new ui[max_n+1];
			buf1 = new ui[max_n+1];
			buf2 = new ui[max_n+1];
			buf3 = new ui[max_n+1];
			buf4 = new ui[max_n+1];
			buf5 = new ui[max_n+1];
			buf6 = new ui[max_n+1];
			buf7 = new ui[max_n+1];
			buf8 = new ui[max_n+1];
			buf9 = new ui[max_n+1];
			buf10 = new ui[max_n+1];
			buf11 = new ui[max_n+1];
			calc_heap = new ListLinearHeap(max_n*4,max_n*4-1);
		}
		else {
			neighbors = new ui[K+2];
			nonneighbors = new ui[K+2];
			nonneighbor_1 = new ui[K+2];
			nonneighbor_2 = new ui[K+2];
			buf = new ui[K+2];
			buf1 = new ui[K+2];
			buf2 = new ui[K+2];
			buf3 = new ui[K+2];
			buf4 = new ui[K+2];
			buf5 = new ui[K+2];
			buf6 = new ui[K+2];
			buf7 = new ui[K+2];
			buf8 = new ui[K+2];
			buf9 = new ui[K+2];
			buf10 = new ui[K+2];
			buf11 = new ui[K+2];
			calc_heap = new ListLinearHeap(K*4,K*4-1);
		}
		vis = new char[max_n+1];
		return;
	}

	ui part_cnt = 0;

	void calc_heap_update(int id,int x){
		int x_key = calc_heap->get_key(x),y_key = 0;
		int new_key = (x_key/4+1)*4;
		if((x_key&3) > 1) new_key += 1;
		// std::cerr << " calc update " << x << " " << x_key << " " << new_key << std::endl;
		calc_heap->increment(x,new_key-x_key);
		if(nonneighbor_1[x] == id) std::swap(nonneighbor_1[x],nonneighbor_2[x]);
		// std::cerr << " calc update " << x << " " << x_key << std::endl;
		int y = nonneighbor_1[x];
		nonneighbor_2[x] = max_n + 5;
		// std::cerr << " calc update " << x << " " << y << std::endl;
		if(y < max_n && degree_in_S[y] == degree_in_S[x]){
			y_key = calc_heap->get_key(y);
			if(nonneighbor_1[y] < max_n && nonneighbor_2[y] < max_n && (y_key&3) == 3)
			calc_heap->decrement(y,1);
		}
		return;
	}

	void calc_ans_on_D_set(ui S_end, ui D_end, ui R_end, ui level){
		ui *tag = buf,*rid = buf1;
		for(ui i = S_end;i < D_end;++i){
			if(nonneighbor_2[SR[i]] != max_n + 5 && SR_rid[nonneighbor_2[SR[i]]] >= D_end) nonneighbor_2[SR[i]] = max_n + 5;
			if(nonneighbor_1[SR[i]] != max_n + 5 && SR_rid[nonneighbor_1[SR[i]]] >= D_end) nonneighbor_1[SR[i]] = max_n + 5;
		}
		for(ui i = S_end;i < D_end;++i){
			tag[SR[i]] = 0;
			if(nonneighbor_2[SR[i]] != max_n + 5){
				tag[SR[i]] = 2;
				if(nonneighbor_1[SR[i]] != max_n + 5 && degree_in_S[SR[i]] < degree_in_S[nonneighbor_1[SR[i]]]) {
					tag[SR[i]] -= 1;
				}
				else if(nonneighbor_2[SR[i]] != max_n + 5 && degree_in_S[SR[i]] < degree_in_S[nonneighbor_2[SR[i]]]) {
					tag[SR[i]] -= 1;
				}
			}
		}
		calc_heap->init_by_minus_three(D_end - S_end,(S_end + 3)*4,SR + S_end,degree_in_S,tag,S_end);
		ui id,key,now_edge = compute_missing_edges_in_S(S_end);
		for(ui i = S_end;i < D_end;++i) tag[SR[i]] = 0;
		for(ui i = S_end;i < D_end;++i){
			calc_heap->pop_min(id,key);
			swap_pos(i,SR_rid[id]);
			now_edge += key/4;
			tag[id] = 1;
			if(nonneighbor_1[id] < max_n && SR_rid[nonneighbor_1[id]] < D_end && tag[nonneighbor_1[id]] == 0) calc_heap_update(id,nonneighbor_1[id]);
			if(nonneighbor_2[id] < max_n && SR_rid[nonneighbor_2[id]] < D_end && tag[nonneighbor_2[id]] == 0) calc_heap_update(id,nonneighbor_2[id]);
			if(now_edge <= K){
				if(two_stage){
					if(i >= best_solution_size){
						store_a_kDefectiveClique(i + 1);
						// std::cerr << "R_end.." << std::endl;
					}
				}
				else{
					ui extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor((sqrt((2*i+1)*(2*i+1)+8*(K-now_edge))-(2*i + 1)) / 2 + eps));
					if(i + extra_num >= best_solution_size){
						store_a_kDefectiveClique(i + 1,extra_num);
					}
				}
			}
			else break;
		}
		return;
	}

	void partition_D_set(ui S_end, ui &D_end, ui &X_end, ui R_end, ui level){
		++part_cnt;
		D_end = X_end = S_end;
		ui *tag = buf,part_flag = 0;
		for(ui i = S_end;i < R_end;++i){
			nonneighbor_1[SR[i]] = max_n + 5;
			tag[SR[i]] = 0;
			if(degree_in_S[SR[i]] == S_end) swap_pos(i,D_end),++D_end;
		}
		bool flag = 0;
		if(S_end != D_end){
			for(ui i = S_end;i < D_end;++i){
				flag = 0;
				++part_flag;
				for(ui j = pstart_R[SR[i]];j < pend[SR[i]] && removed_level[edges[j]] >= level;++j){
					tag[edges[j]] = part_flag;
				}
				for(ui j = i+1;j < D_end;++j){
					if(tag[SR[j]] != part_flag){
						if(nonneighbor_1[SR[j]] != max_n + 5) nonneighbor_1[SR[j]] = max_n + 6;
						else nonneighbor_1[SR[j]] = SR[i];
						if(nonneighbor_1[SR[i]] != max_n + 5) flag = 1;
						else nonneighbor_1[SR[i]] = SR[j];
					}
				}
				if(flag == 1) continue;
				swap_pos(i,X_end);++X_end;
			}
			D_end = X_end;
		}
		else{
			for(ui i = S_end;i < R_end;++i){
				flag = 0;
				++part_flag;
				for(ui j = pstart_R[SR[i]];j < pend[SR[i]] && removed_level[edges[j]] >= level;++j){
					tag[edges[j]] = part_flag;
				}
				for(ui j = i+1;j < R_end;++j){
					if(tag[SR[j]] != part_flag){
						if(nonneighbor_1[SR[j]] != max_n + 5) nonneighbor_1[SR[j]] = max_n + 6;
						else nonneighbor_1[SR[j]] = SR[i];
						if(nonneighbor_1[SR[i]] != max_n + 5) flag = 1;
						else nonneighbor_1[SR[i]] = SR[j];
					}
				}
				if(flag == 1) continue;
				swap_pos(i,X_end);++X_end;
			}
		}
		return;
	}

	ui choose_vertex_by_neighbours(ui &S_end,ui D_end, ui &R_end, ui level){
		if(S_end == 0) return n;
		ui min_C_u = SR[S_end],max_C_u = SR[S_end],max_S_u = SR[0];
		for(ui i = 1; i < S_end; ++ i) if(degree_in_S[SR[i]] < degree_in_S[max_S_u]) max_S_u = SR[i];
		for(ui i = D_end; i < R_end; ++ i){
			if(degree[SR[i]] > degree[min_C_u]) {min_C_u = SR[i];break;}
		}
		if(R_end - degree[min_C_u] - 1 <= S_end - degree_in_S[max_S_u] - 1) return min_C_u;
		return n;
	}

    void load_graph(ui _n, ui *_pstart, ui *_pend, ui *_edges, ui _two_stage,std::vector<ui> zero_point_array) {
		two_stage = _two_stage;
    	n = _n;
        ui m = 0;
        for(ui i = 0;i < n;i ++) m += _pend[i] - _pstart[i];

        assert(pstart == NULL);
		if(pstart != NULL) delete[] pstart;pstart = NULL;
		if(edges != NULL) delete[] edges;edges = NULL;

        pstart = new ui[n+1];
        edges = new ui[m];

		flow_prune = (n <= 1000 ? 1 : 0);
		look_ahead_prune = 0;
		#ifdef COLOR_PRUNE
		color_prune = 1;
		#else
		color_prune = 0;
		#endif
		if(solver != NULL) delete solver;solver = NULL;
		solver = new DinicSolver(_n,look_ahead_prune);

        m = 0;
        for(ui i = 0;i < n;i ++) {
        	pstart[i] = m;
        	for(ui j = _pstart[i];j < _pend[i];j ++){
				edges[m ++] = _edges[j];
			}
        }
        pstart[n] = m;

        // printf("load graph of size n=%u, m=%u (undirected), density=%.5lf\n", n, m/2, double(m)/n/(n-1));
    }

	void load_graph(ui _n,std::vector<std::pair<ui,ui>> vp){
		n = _n;
		if(matrix != NULL) delete[] matrix;
		if(n > history_max_n) history_max_n = n,matrix = new bool[n*n];
		for(ui i = 0,lim = vp.size();i < lim;++i){
			matrix[n*vp[i].first + vp[i].second] = 1;
			matrix[n*vp[i].second + vp[i].first] = 1;
		}
		return;
	}

    void kDefectiveClique(ui _K, ui _UB, std::vector<ui> &kDC) {
        K = _K;
        UB = _UB;
        if(K == 0) {
        	printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
        	return ;
        }
        best_solution_size = kDC.size();

		// std::cerr << n << " " << max_n << std::endl;

		memset(vis, 0, sizeof(char)*(n+1));
		memset(degree_in_S, 0, sizeof(ui)*(n+1));
		for(ui i = 0;i < n;i ++) removed_level[i] = n;
		for(ui i = 0;i < n;i ++) SR[i] = SR_rid[i] = i;
		for(ui i = 0;i < n;i ++) {
			pstart_R[i] = pstart[i];
			pend[i] = pstart[i+1];
		}
		while(!Qv.empty()) Qv.pop();

		zero_degree_vertices.clear();

		for(ui i = 0;i < n;i ++) {
			degree[i] = pstart[i+1] - pstart[i];
			// printf("degree %u: %u\n", i, degree[i]);
			if(degree[i] + K < best_solution_size) {
				removed_level[i] = 0;
				Qv.push(i);
			}
			if(!two_stage && degree[i] == 0){
				removed_level[i] = 0;
				Qv.push(i);
				zero_degree_vertices.push_back(i);
			}
		}

		ui R_end = n;
		remove_vertices_and_prune(0, R_end, 0);

		if(R_end != n) printf("Initially pruned %u vertices!\n", n - R_end);

		if(two_stage){
			for(ui i = pstart[0];i < pend[0];++i) ++degree_in_S[edges[i]];
			if(SR_rid[0] < R_end)
			BB_search(1,R_end,1);
		}

		else{
			BB_search(0,R_end,1);
		}

		// std::cerr << "branch_tot:" << branch_tot << std::endl;

		std::sort(kDC.begin(),kDC.end());
        if(best_solution_size > kDC.size()) {
            kDC.clear();
            for(int i = 0;i < best_solution_size;i ++) kDC.push_back(best_solution[i]);
        }
		// std::cerr << "branch_tot:" << branch_tot << std::endl;
		return;
    }

    int main(int argc, char *argv[]) {
	    if(argc < 3) {
		    printf("Usage: [1]exe [2]k [3]dir [4 option] lb_of_max_kDefectiveClique_size\n");
		    return 0;
	    }
        K = atoi(argv[1]);
        UB = 100000000;
        if(K == 0) {
        	printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
        	return 0;
        }
        std::vector<ui> kDC;
        if(argc >= 4) {
        	best_solution_size = atoi(argv[3]);
        	printf("initial lb: %u\n", best_solution_size);
        	for(ui i = 0;i < best_solution_size;i ++) kDC.pb(0);
        }
        readGraph_binary(argv[2]);
        printf("Finish reading graph\n");
        Timer t;
		two_stage = 0;
        kDefectiveClique(K, UB, kDC);
        printf("Maximum %u-DefectiveClique-modified size: %u, time excluding reading: %s (micro seconds)\n", K, best_solution_size, Utility::integer_to_string(t.elapsed()).c_str());
        return 0;
    }

private:
    void readGraph_binary(char* dir) {
	    FILE *f = Utility::open_file( (std::string(dir) + std::string("/b_degree.bin")).c_str(), "rb");

	    ui tt;
	    fread(&tt, sizeof(ui), 1, f);
	    if(tt != sizeof(ui)) {
		    printf("sizeof int is different: edge.bin(%u), machine(%lu)\n", tt, sizeof(ui));
		    return ;
	    }
	    ui m;
	    fread(&n, sizeof(ui), 1, f);
	    fread(&m, sizeof(ui), 1, f);
	    printf("n = %u, m = %u\n", n, m/2);

	    ui *degree = new ui[n];
	    fread(degree, sizeof(ui), n, f);
	    fclose(f);

	    f = Utility::open_file( (std::string(dir) + std::string("/b_adj.bin")).c_str(), "rb");
	    pstart = new ui[n+1];
	    edges = new ui[m];
	    pstart[0] = 0;
	    for(ui i = 0;i < n;i ++) {
	    	pstart[i+1] = pstart[i] + degree[i];
	    	if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
		}
	    fclose(f);

	    delete[] degree;
    }

    void initialization(ui S_end, ui R_end, ui level) { // only reorganizes neighbors of vertices of R
    	assert(level > 0);
    	while(!Qv.empty()) Qv.pop();
    	ui start_idx = S_end;
    	if(start_idx > 0) -- start_idx;
    	//print_array("SR", SR, 0, R_end, 0);
    	for(ui i = start_idx;i < R_end;i ++) {
    		ui u = SR[i];
    		ui non_neighbors_n = 0, end = pstart_R[u];
    		//if(u == 0) print_neighbors(u, pstart, pend, edges);
    		for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level-1;j ++) {
    			assert(removed_level[edges[j]] == level-1||removed_level[edges[j]] == n);
    			if(removed_level[edges[j]] >= level) {
    				edges[end ++] = edges[j];
    				if(SR_rid[edges[j]] < S_end) {
    					std::swap(edges[pstart_R[u]], edges[end-1]);
    					++ pstart_R[u];
    				}
    			}
    			else{
					nonneighbors[non_neighbors_n ++] = edges[j];
				}
    		}
#ifndef NDEBUG
    		if(degree[u] != end - pstart[u]) {
    			printf("u: %u, pstart[u+1]-pstart[u]: %u, degree[u]: %u, end - pstart[u]: %u\n", u, pstart[u+1]-pstart[u], degree[u], end - pstart[u]);
    		}
    		assert(degree[u] == end-pstart[u]);
#endif
    		for(ui j = 0;j < non_neighbors_n;j ++){
				edges[end ++] = nonneighbors[j];
			}
    		assert((end < pend[u]&&removed_level[edges[end]] < level-1)||end == pend[u]);
#ifndef NDEBUG
    		for(ui j = end;j < pend[u];j ++) {
    			if(removed_level[edges[j]] >= level) printf("removed_level[edges[j]]: %u, level: %u\n", removed_level[edges[j]], level);
    			assert(removed_level[edges[j]] < level);
    		}
#endif
    	}
    }

	ui debug_node_degree(ui u,ui R_end,ui level){
		ui ans = 0;
		for(ui i = pstart[u];i < pend[u];++i){
			if(SR_rid[edges[i]] < R_end){
				++ans;
			}
		}
		return ans;
	}

    void compute_a_heuristic_solution_and_prune(ui &R_end, ui level) {
    	// the following computes the degeneracy ordering and a heuristic solution
#ifndef NDEBUG
    	for(ui i = 0;i < R_end;i ++) {
    		ui u = SR[i];
    		assert(degree[u] + K >= best_solution_size);
    		assert(pstart[u] == pstart_R[u]);
    		ui end = pstart[u];
    		while(end < pend[u]&&removed_level[edges[end]] >= level) ++ end;
    		for(ui j = end;j < pend[u];j ++) {
    			if(removed_level[edges[j]] >= level) printf("removed_level[edges[j]]: %u, level: %u\n", removed_level[edges[j]], level);
    			assert(removed_level[edges[j]] < level);
    		}
    		if(degree[u] != end - pstart[u]) printf("degree[u]: %u, %u\n", degree[u], end-pstart[u]);
    		assert(degree[u] == end - pstart[u]);
    	}
#endif
		ui *core = neighbors;
		ui *rid = nonneighbors;
		ui *id = buf;
		ui *t_degree = buf1;
		ui total_edges = 0;
		for(ui i = 0;i < R_end;i ++) {
			id[i] = 0;
			t_degree[SR[i]] = degree[SR[i]];
			assert(t_degree[SR[i]] < R_end);
			total_edges += degree[SR[i]];
		}
		for(ui i = 0;i < R_end;i ++) ++ id[t_degree[SR[i]]];
		for(ui i = 1;i < R_end;i ++) id[i] += id[i-1];

		for(ui i = 0;i < R_end;i ++) rid[SR[i]] = -- id[t_degree[SR[i]]];
		for(ui i = 0;i < R_end;i ++) id[rid[SR[i]]] = SR[i];

		ui *degree_start = buf2;
		for(ui i = 0, j = 0;i <= R_end;i ++) {
			while(j < R_end&&t_degree[id[j]] < i) ++ j;
			degree_start[i] = j;
		}

		ui max_core = 0;
		for(ui i = 0;i < R_end;i ++) {
			ui u = id[i];
			assert(degree_start[t_degree[u]] == i);
			if(t_degree[u] > max_core) max_core = t_degree[u];
			core[u] = max_core;

			long long t_n = R_end - i;

			if(t_n*(t_n-1)/2 <= total_edges/2 + K&&R_end - i > best_solution_size) {
				best_solution_size = R_end - i;
				for(ui j = i;j < R_end;j ++) best_solution[j-i] = id[j];
				printf("Degen find a solution of size %u\n", best_solution_size);
			}

			++ degree_start[t_degree[u]];
			if(t_degree[u] == 0) continue;

			degree_start[t_degree[u]-1] = degree_start[t_degree[u]];
			for(ui j = pstart[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(rid[edges[j]] > i) {
				ui v = edges[j];
				ui pos1 = degree_start[t_degree[v]], pos2 = rid[v];
				std::swap(id[pos1], id[pos2]);
				rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
				++ degree_start[t_degree[v]];
				-- t_degree[v];
				total_edges -= 2;
			}
		}

		assert(Qv.empty());
		for(ui i = 0;i < R_end;i ++) if(core[SR[i]] + K < best_solution_size) {
			assert(removed_level[SR[i]] > level);
			removed_level[SR[i]] = level;
			Qv.push(SR[i]);
		}
    }

    ui degeneracy_ordering_and_coloring_adj(ui S_end, ui R_end, ui level, ui *color) {
		ui *rid = buf;
		ui *id = buf1;
		ui *t_degree = color;
		ui max_degree = 0;
		for(ui i = S_end;i < R_end;i ++) {
			ui &d = t_degree[SR[i]] = 0;
			for(ui j = pstart_R[SR[i]];j < pend[SR[i]]&&removed_level[edges[j]] >= level;j ++) {
				assert(SR_rid[SR[i]] >= S_end);
				if(SR_rid[edges[j]] < R_end) ++ d;
			}
			assert(t_degree[SR[i]] == degree[SR[i]] - degree_in_S[SR[i]]);
			if(d > max_degree) max_degree = d;
		}
		memset(id, 0, sizeof(ui)*(max_degree+1));
		for(ui i = S_end;i < R_end;i ++) ++ id[t_degree[SR[i]]];
		for(ui i = 1;i <= max_degree;i ++) id[i] += id[i-1];

		for(ui i = S_end;i < R_end;i ++) rid[SR[i]] = -- id[t_degree[SR[i]]];
		for(ui i = S_end;i < R_end;i ++) id[rid[SR[i]]] = SR[i];

		ui *degree_start = buf2;
		for(ui i = 0, j = 0;i <= max_degree;i ++) {
			while(j < R_end&&t_degree[id[j]] < i) ++ j;
			degree_start[i] = j;
		}

		ui max_core = 0;
		for(ui i = 0;i < R_end-S_end;i ++) {
			ui u = id[i];
			assert(degree_start[t_degree[u]] == i);
			if(t_degree[u] > max_core) max_core = t_degree[u];

			++ degree_start[t_degree[u]];
			if(t_degree[u] == 0) continue;

			degree_start[t_degree[u]-1] = degree_start[t_degree[u]];
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end&&rid[edges[j]] > i) {
				ui v = edges[j];
				ui pos1 = degree_start[t_degree[v]], pos2 = rid[v];
				std::swap(id[pos1], id[pos2]);
				rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
				++ degree_start[t_degree[v]];
				-- t_degree[v];
			}
		}

    	ui max_color = 0;
    	for(ui i = R_end-S_end;i > 0;i --) {
    		ui u = id[i-1];
    		for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) {
    			if(SR_rid[edges[j]] < R_end&&rid[edges[j]] >= i) vis[color[edges[j]]] = 1;
    		}
    		for(ui j = 0;;j ++) if(!vis[j]) {
    			color[u] = j;
    			if(j > max_color) max_color = j;
    			break;
    		}
    		for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) {
    			if(SR_rid[edges[j]] < R_end&&rid[edges[j]] >= i) vis[color[edges[j]]] = 0;
    		}
    	}

    	return max_color + 1;
    }

    bool is_kDefectiveClique(ui R_end) {
    	ui total_edges = 0;
    	long long all_edges = R_end;
    	for(ui i = 0;i < R_end;i ++){
			total_edges += degree[SR[i]];
		}
    	return all_edges*(R_end-1)/2 <= total_edges/2 + K;
    }

	bool is_kDefectiveClique(ui R_end,ui &zero_end){
    	ui total_edges = 0,max_end = zero_end;
		zero_end = 0;
    	for(ui i = 0;i < R_end;i ++) total_edges += R_end - degree[SR[i]] - 1;
		if(total_edges > 2*K) return false;
		total_edges /= 2;
		for(ui i = 1;i <= max_end;i ++){
			total_edges += R_end + i - 1;
			if(total_edges > K) break;
			zero_end = i;
		}
		total_edges /= 2;
		return true;
	}

    void store_a_kDefectiveClique(ui S_end) {
		// std::cerr << best_solution_size << " " << S_end << std::endl;
    	assert(S_end > best_solution_size);
		best_solution_size = S_end;
		for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
		ui* tag = buf,lost_edges=0,lost_deg = 0;
		memset(tag,0,sizeof(ui)*(S_end+1));
		// printf("kDefectiveClique begin\n");
		for(ui i=0;i<S_end;++i){
			ui u = SR[i],cnt=0;
			for(ui j=0;j<S_end;++j) tag[SR[j]] = 0;
			tag[SR[i]] = 1;
			for(ui j=pstart[u];j<pend[u];++j){
				if(SR_rid[edges[j]]<S_end) ++cnt;
				tag[edges[j]] = 1;
			}
			for(ui j=0;j<S_end;++j){
				if(tag[SR[j]] == 0 && SR[i] < SR[j]){
					++lost_deg;
				}
			}
			// printf("%u ",u);
		}
		// printf("\n");
		// printf("new_size(bitset):%u,lost_edges: %u\n",S_end, lost_edges/2);
		// std::cerr << compute_missing_edges_in_S(S_end) << std::endl;
#ifndef NDEBUG
		printf("Find a kDefectiveClique of size: %u\n", best_solution_size);
#endif
    }

    void store_a_kDefectiveClique(ui S_end,ui zero_num) {
    	assert(S_end > best_solution_size);
		best_solution_size = S_end + zero_num;
		for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
		for(ui i = 0;i < zero_num;++i) best_solution[S_end+i] = zero_degree_vertices[i];
		ui* tag = buf,lost_edges=0,lost_deg = 0;
		memset(tag,0,sizeof(ui)*(S_end+zero_num+1));
		for(ui i=0;i<S_end+zero_num;++i){
			ui u = SR[i],cnt=0;
			for(ui j=0;j<S_end;++j) tag[SR[j]] = 0;
			tag[SR[i]] = 1;
			for(ui j=pstart[u];j<pend[u];++j){
				if(SR_rid[edges[j]]<S_end) ++cnt;
				tag[edges[j]] = 1;
			}
			lost_edges += S_end-cnt-1;
		}
#ifndef NDEBUG
		printf("Find a kDefectiveClique of size: %u\n", best_solution_size);
#endif
    }

	int move_noneighbours_of_u_from_R_to_S(ui u, ui &S_end, ui &R_end, ui level) {
		ui* tag = buf;
		int tot = 0;
		memset(tag,0,sizeof(ui)*(R_end+1));
		for(ui i=pstart_R[u];i<pend[u]&&removed_level[edges[i]]>=level;++i){
			tag[SR_rid[edges[i]]] = 1;
		}
        ui missing_edges = compute_missing_edges_in_S(S_end);
		for(ui i=S_end;i<R_end;++i){
			if(!tag[i]){
				if(S_end - degree_in_S[SR[i]] + missing_edges > K) return -tot;
				++tot;
				missing_edges += S_end - degree_in_S[SR[i]];
				move_u_from_R_to_S(SR[i],S_end,R_end,level);
				if(remove_vertices_and_prune(S_end,R_end,level)) return -tot;
			}
		}
		return tot;
	}

	int move_noneighbours_of_u_from_R_to_S_except_v(ui u, ui v,std::vector<ui> vec, ui &S_end, ui &R_end, ui level){
		int tot = 0;
        ui missing_edges = compute_missing_edges_in_S(S_end);
		for(ui i = 0,lim = vec.size();i<lim;++i){
			if(vec[i] == v) continue;
			if(SR_rid[vec[i]] >= R_end) return -tot;
			if(S_end - degree_in_S[vec[i]] + missing_edges > K) return -tot;
			++tot;
			missing_edges += S_end - degree_in_S[vec[i]];
			move_u_from_R_to_S(vec[i],S_end,R_end,level);
			if(remove_vertices_and_prune(S_end,R_end,level)) return -tot;
		}
		return tot;
	}

	void extract_non_neighbours_of_u_from_R(ui u, ui &S_end, ui &R_end, ui level, std::vector<ui> &vec) {
		ui* tag = buf;
		memset(tag,0,sizeof(ui)*(R_end+1));
		tag[SR_rid[u]] = 1;
		for(ui i=pstart_R[u];i<pend[u]&&removed_level[edges[i]]>=level;++i){
			tag[SR_rid[edges[i]]] = 1;
		}
		for(ui i=S_end;i<R_end;++i){
			if(!tag[i]){
				vec.push_back(SR[i]);
			}
		}
		return;
	}

	
	void extract_non_neighbours_of_u(ui u, ui &S_end, ui &R_end, ui level, std::vector<ui> &vec) {
		ui* tag = buf;
		memset(tag,0,sizeof(ui)*(R_end+1));
		tag[SR_rid[u]] = 1;
		for(ui i=pstart_R[u];i<pend[u]&&removed_level[edges[i]]>=level;++i){
			tag[SR_rid[edges[i]]] = 1;
		}
		for(ui i=0;i<R_end;++i){
			if(!tag[i]){
				vec.push_back(SR[i]);
			}
		}
		return;
	}

	
	
	void extract_neighbours_of_u_from_vec(ui u, ui &S_end, ui &R_end, ui level,std::vector<ui> &src_vec, std::vector<ui> &vec) {
		ui* tag = buf;
		memset(tag,0,sizeof(ui)*(R_end+1));
		tag[SR_rid[u]] = 1;
		for(ui i=pstart_R[u];i<pend[u]&&removed_level[edges[i]]>=level;++i){
			tag[SR_rid[edges[i]]] = 1;
		}
		for(ui i=0,lim=src_vec.size();i<lim;++i){
			if(!tag[src_vec[i]]){
				vec.push_back(src_vec[i]);
			}
		}
		return;
	}

	
	int choose_branch_vertex_based_on_MDC_2_non(ui S_end,ui D_end, ui R_end,ui level){
		for(ui i=S_end;i<D_end;++i){
			if(degree[SR[i]]+3>=R_end && degree_in_S[SR[i]] + 1 >= S_end){
				// std::cerr << SR[i] << " " << degree[SR[i]] << " " << R_end << std::endl;
				return SR[i];
			}
		}
		return n;
	}
	
	void partition_C1_set(ui S_end,ui &D_end,ui R_end,ui level){
		D_end = S_end;
		for(ui i = S_end;i < R_end;++i){
			if(degree_in_S[SR[i]] + 1 <= S_end) swap_pos(i,D_end),++D_end;
		}
		return;
	}

	ui check_tot = 0;

	ui in_check = 0;
	ui out_check = 0;

    void BB_search(ui S_end, ui R_end, ui level) {
		++branch_tot;
		// if(S_end == 0)
		// std::cerr << R_end << std::endl;
		// std::cerr << S_end << " " << R_end << " " << level << " " << zero_degree_vertices.size() << std::endl;
		// for(ui i = 0;i < S_end;++i){
		// 	std::cerr << SR[i] << " ";
		// }
		// std::cerr << std::endl;
		bool t_flag = false;
		
		// for(ui i = S_end;i < R_end;++i){
		// 	std::cerr << S_end - degree_in_S[SR[i]] << " ";
		// }
		// std::cerr << std::endl;
		// if(t_flag == true)
		// std::cerr << " (S_end) " << S_end << " (R_begin) ";

		// t_flag = false;

		// for(ui i = 0;i < R_end;++i){
		// 	if(my_vec.find(SR[i]) != my_vec.end())
		// 	std::cerr << SR[i] << " ",t_flag = true;
		// }
		// if(t_flag == true)
		// std::cerr << std::endl;
		++branch_tot;
		ui cnt = 0;
#ifndef NDEBUG
    	assert(compute_missing_edges_in_S(S_end) <= K);
    	for(ui i = 0;i < R_end;i ++) assert(degree[SR[i]] + K >= best_solution_size);
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif
		if(two_stage){
			if(best_solution_size >= UB||R_end <= best_solution_size) return ;
			fflush(stdout);
			if(S_end > best_solution_size){
				store_a_kDefectiveClique(S_end);
						// std::cerr << "R_end.." << std::endl;
			}
			if(R_end > best_solution_size&&is_kDefectiveClique(R_end)){
				store_a_kDefectiveClique(R_end);
			}
			if(R_end <= best_solution_size+1) return ;
		}
		else{
			ui S_deg = compute_missing_edges_in_S(S_end);
			if(best_solution_size >= UB) return;
			ui extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(S_end?(sqrt((2*S_end-1)*(2*S_end-1)+8*(K-S_deg))-(2*S_end - 1)) / 2:0 + eps));
			if(S_end + extra_num > best_solution_size){
				store_a_kDefectiveClique(S_end,extra_num);
			}
			extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(R_end?(sqrt((2*R_end-1)*(2*R_end-1)+8*(K-S_deg))-(2*R_end - 1)) / 2:0 + eps));
			if(R_end + extra_num > best_solution_size && is_kDefectiveClique(R_end,extra_num) && R_end + extra_num > best_solution_size ){
				store_a_kDefectiveClique(R_end,extra_num);
			}
			if(R_end + extra_num <= best_solution_size + 1) return ;
		}
		// std::cerr << "checkmate 1" << std::endl;
        initialization(S_end, R_end, level);
#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			for(ui j = pstart[u];j < pstart_R[u];j ++) assert(SR_rid[edges[j]] < S_end);
			assert(degree_in_S[u] == pstart_R[u] - pstart[u]);
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end) {
				assert(SR_rid[edges[j]] >= S_end);
				++ cnt;
			}
			assert(degree[u] == cnt + degree_in_S[u]);
		}
#endif

		// std::cerr << "checkmate 0.5" << std::endl;
		ui flag = 0,old_R_end = R_end,loss_edge = compute_missing_edges_in_S(S_end);
		for(ui i = S_end;i < R_end;i ++){
			ui u = SR[i];
			if(degree[u] + 2 >= R_end && S_end - degree_in_S[u] + loss_edge <= K){
				// std::cerr << u << " " << degree[u] + 2 << " -deg- " << R_end << std::endl;
				loss_edge += S_end - degree_in_S[u];
				++flag;
				move_u_from_R_to_S(u,S_end,R_end,level);
				if(S_end > best_solution_size)  store_a_kDefectiveClique(S_end);
				if(remove_vertices_and_prune(S_end, R_end, level)){
					while(flag--) move_u_from_S_to_R(S_end,R_end,level);
					restore_R(S_end,R_end,old_R_end,level);
					return;
				}
			}
		}

		// std::cerr << "checkmate 1" << std::endl;
		if(flag){
			BB_search(S_end,R_end,level+1);
			while(flag--) move_u_from_S_to_R(S_end,R_end,level);
			restore_R(S_end,R_end,old_R_end,level);
			return;
		}
		// std::cerr << "checkmate 2" << std::endl;
        old_R_end = R_end;
        assert(Qv.empty());
        bool terminate = collect_removable_vertices_based_on_degree_in_S_opt(S_end, R_end, level);
        if(terminate||remove_vertices_and_prune(S_end, R_end, level)) {
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }
        // if(S_end == 0) {
        // 	assert(Qv.empty());
        // 	compute_a_heuristic_solution_and_prune(R_end, level);
        // 	if(level == 1) printf("First level kdefectiveclique size: %lu\n", best_solution_size);
        // 	if(remove_vertices_and_prune(S_end, R_end, level)) {
        // 		restore_R(S_end, R_end, old_R_end, level);
        // 		return ;
        // 	}
        // }
		// std::cerr << "checkmate 3" << std::endl;
        if(S_end > 0) {
        	assert(Qv.empty());
        	collect_removable_vertices_based_on_vertex_pair(S_end, R_end, level); // remove vertices from R
			if(remove_vertices_and_prune(S_end, R_end, level)) {
        		restore_R(S_end, R_end, old_R_end, level);
				return ;
			}
        }
		if(degree_based_prune(S_end, R_end, level)) {
			restore_R(S_end, R_end, old_R_end, level);
			return ;
		}
		// std::cerr << "checkmate 4" << std::endl;
		if(two_stage){
			if(R_end > best_solution_size&&is_kDefectiveClique(R_end)){
				store_a_kDefectiveClique(R_end);
			}
			if(R_end <= best_solution_size+1){
				restore_R(S_end, R_end, old_R_end, level);
				return;
			}
		}
		else{
			ui extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(R_end?(sqrt((2*R_end-1)*(2*R_end-1)+8*K)-(2*R_end - 1)) / 2:0 + eps));
			if(R_end + extra_num > best_solution_size&&is_kDefectiveClique(R_end,extra_num) && R_end + extra_num > best_solution_size){
				store_a_kDefectiveClique(R_end,extra_num);
			}
			if(R_end + extra_num <= best_solution_size+1){
				restore_R(S_end, R_end, old_R_end, level);
				return;
			}
		}
		// std::cerr << "checkmate 5" << std::endl;
		if(R_end == S_end) return;

		// if(S_end == 0)
		// std::cerr << S_end << " -S_end- " << R_end << std::endl;

		flag = 0;
		if(color_prune)
		while(flag = coloring_based_prune(S_end, R_end, level)){
			// if(S_end == 0)
			// std::cerr << S_end << " -col- " << R_end << std::endl;
			// std::cerr << "checkmate 5.5 " << flag << std::endl;
			if(flag == 2){
				// std::cerr << " checkmate 6.5 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
				restore_R(S_end, R_end, old_R_end, level);
				return;
			}
			else{
				// std::cerr << "checkmate 6" << std::endl;
				if(remove_vertices_and_prune(S_end, R_end, level)) {
					// std::cerr << " checkmate 5.7 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
					restore_R(S_end, R_end, old_R_end, level);
					return ;
				}
				if(two_stage){
					// std::cerr << " checkmate 5.6 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
					if(R_end > best_solution_size&&is_kDefectiveClique(R_end)){
						store_a_kDefectiveClique(R_end);
						// std::cerr << "R_end.." << std::endl;
					}
					if(R_end <= best_solution_size+1) {
						// std::cerr << " checkmate 5.8 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
						restore_R(S_end, R_end, old_R_end, level);
						return ;
					}
				}
				else{
					// std::cerr << " checkmate 5.6 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
					ui extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(R_end?(sqrt((2*R_end-1)*(2*R_end-1)+8*K)-(2*R_end - 1)) / 2:0 + eps));
					if(R_end + extra_num > best_solution_size&&is_kDefectiveClique(R_end,extra_num) && R_end + extra_num > best_solution_size){
						store_a_kDefectiveClique(R_end);
						// std::cerr << "R_end.." << std::endl;
					}
					if(R_end + extra_num <= best_solution_size+1) {
						// std::cerr << " checkmate 6.5 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
						restore_R(S_end, R_end, old_R_end, level);
						// std::cerr << " checkmate 5.8 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
						return ;
					}
				}
			}
		}

		
#ifdef MDC_BRANCH //重写MDC的分支
		ui D_end = S_end,X_end = S_end,u=n,t_old_R_end,m_edges = compute_missing_edges_in_S(S_end),temp_m_edges = m_edges,t_old_R_end_2;
		partition_C1_set(S_end,D_end,R_end,level);
		u = choose_branch_vertex_based_on_MDC_2_non(S_end,D_end, R_end, level);
		// std::cerr << "checkmate 6" << std::endl;
		// std::cerr << "choose u:" << u << " " << n << " " << S_end << " " << R_end << degree[u] << " " << R_end << std::endl;
		if(u != n){
			// std::cerr << "2-non branch" << std::endl;
			if(degree[u] + 2 >= R_end){
				// std::cerr << "plan 1:" << u << std::endl;
				if(m_edges + S_end - degree_in_S[u] <= K){
					move_u_from_R_to_S(u,S_end,R_end,level);
					if(remove_vertices_and_prune(S_end, R_end, level)){
						move_u_from_S_to_R(S_end,R_end,level);
						restore_R(S_end,R_end,old_R_end,level);
						// std::cerr << "return.. plan 1 n finish" << std::endl;
						return;
					}
					// std::cerr << "1 branch" << std::endl;
					BB_search(S_end,R_end,level+1);
					move_u_from_S_to_R(S_end,R_end,level);
					restore_R(S_end,R_end,old_R_end,level);
				}
				// std::cerr << "return.. plan 1 y finish" << std::endl;
				return;
			}else{
				// std::cerr << "plan 2:" << u << " " << degree[u] << " " << R_end << std::endl;
				if(m_edges + S_end - degree_in_S[u] <= K){
					t_old_R_end = R_end;
					move_u_from_R_to_S(u,S_end,R_end,level);
					// std::cerr << "2 branch" << std::endl;
					if(!remove_vertices_and_prune(S_end, R_end, level)){
						BB_search(S_end,R_end,level+1);
					}
					move_u_from_S_to_R(S_end,R_end,level);
					restore_R(S_end,R_end,t_old_R_end,level);
				}
				std::vector<ui> non_neigh_vec,non_neigh_v,non_neigh_w;
				extract_non_neighbours_of_u(u,S_end,R_end,level,non_neigh_vec);
				if(non_neigh_vec.size() == 2){
					// std::cerr << "plan 2 A:" << std::endl;
					int v = non_neigh_vec[0],w = non_neigh_vec[1],flag = 0;
					extract_non_neighbours_of_u_from_R(v,S_end,R_end,level,non_neigh_v);
					extract_non_neighbours_of_u_from_R(w,S_end,R_end,level,non_neigh_w);
					for(auto p:non_neigh_v) if(p==w) {flag=1;break;}
					if(2*S_end - degree_in_S[v] - degree_in_S[w] == 0 && flag == 0){
						// std::cerr << "into" << std::endl;
						removed_level[u] = level;
						Qv.push(u);
						for(auto p:non_neigh_v) if(removed_level[p] > level) {
							removed_level[p] = level;
							Qv.push(p);
						}
						for(auto p:non_neigh_w) if(removed_level[p] > level) {
							removed_level[p] = level;
							Qv.push(p);
						}

						temp_m_edges = m_edges;
						if(remove_vertices_and_prune(S_end,R_end,level) || removed_level[w] >= level){
							restore_R(S_end,R_end,old_R_end,level);
							return;
						}
						if(temp_m_edges + S_end - degree_in_S[w] > K){
							restore_R(S_end,R_end,old_R_end,level);
							return;
						}
						temp_m_edges += S_end - degree_in_S[w];
						move_u_from_R_to_S(w,S_end,R_end,level);

						if(remove_vertices_and_prune(S_end,R_end,level) || removed_level[v] >= level){
							restore_R(S_end,R_end,old_R_end,level);
							return;
						}
						if(temp_m_edges + S_end - degree_in_S[v] > K){
							restore_R(S_end,R_end,old_R_end,level);
							return;
						}
						temp_m_edges += S_end - degree_in_S[v];
						move_u_from_R_to_S(v,S_end,R_end,level);

						if(!remove_vertices_and_prune(S_end,R_end,level))
						BB_search(S_end,R_end,level+1);

						move_u_from_S_to_R(S_end,R_end,level);
						move_u_from_S_to_R(S_end,R_end,level);
					}
					// std::cerr << "return.. A" << std::endl;
				}
				else if(non_neigh_vec.size() == 1 && S_end - degree_in_S[u] == 1){
					// std::cerr << "plan 2 B:" << std::endl;
					int v = non_neigh_vec[0];
					if(S_end == degree_in_S[v]){
						removed_level[u] = level;
						Qv.push(u);
						extract_non_neighbours_of_u_from_R(v,S_end,R_end,level,non_neigh_v);
						for(auto p:non_neigh_v) if(removed_level[p] > level){
							removed_level[p] = level;
							Qv.push(p);
						}
						if(remove_vertices_and_prune(S_end,R_end,level) || removed_level[v] >= level || temp_m_edges + S_end - degree_in_S[v] > K){
							restore_R(S_end,R_end,old_R_end,level);
							return;
						}
						move_u_from_R_to_S(v,S_end,R_end,level);
						// std::cerr << "B branch" << std::endl;
						if(!remove_vertices_and_prune(S_end,R_end,level))
						BB_search(S_end,R_end,level+1);
						move_u_from_S_to_R(S_end,R_end,level);
					}
					// std::cerr << "return.. B" << std::endl;
				}
				else{
					printf("non_neigh_vec size: %lu %lu %lu\n",u,non_neigh_vec.size(),S_end - degree_in_S[u]);
					puts("Wrong!!!! please try to debug.\n");
				}
				restore_R(S_end,R_end,old_R_end,level);
				return;
			}
		}else if(D_end != R_end){
			u = SR[D_end];
			// std::cerr << "D_end branch " << u << std::endl;
			for(int i = D_end;i < R_end;++i)if(removed_level[SR[i]] > level && degree_in_S[SR[i]] < degree_in_S[u]){
				u = SR[i];
			}
			t_old_R_end = R_end;
			move_u_from_R_to_S(u,S_end,R_end,level);
			// std::cerr << "D branch 1" << std::endl;
			if(!remove_vertices_and_prune(S_end,R_end,level))
			BB_search(S_end,R_end,level+1);
			move_u_from_S_to_R(S_end,R_end,level);
			restore_R(S_end,R_end,t_old_R_end,level);
			if(removed_level[u] > level){
				Qv.push(u);
				removed_level[u] = level;
			}
			if(!remove_vertices_and_prune(S_end,R_end,level))
			BB_search(S_end,R_end,level+1);
			restore_R(S_end,R_end,old_R_end,level);
			return;
		}else{
			u = SR[S_end];
			for(int i=S_end+1;i<R_end;++i)if(degree[SR[i]]-degree_in_S[SR[i]] > degree[u]-degree_in_S[u]){
				u = SR[i];
			}
			std::vector<ui> P1 = {u},P2,non_neigh_u;
			bool flagNnbSu = (S_end - degree_in_S[u] == 1);

			extract_non_neighbours_of_u_from_R(u,S_end,R_end,level,non_neigh_u);
			for(auto p:non_neigh_u){
				if(flagNnbSu && S_end == degree_in_S[p])
					P1.push_back(p);
				else
					P2.push_back(p);
			}

			for(auto p:P1)if(SR_rid[p] >= S_end && SR_rid[p] < R_end){
				if(m_edges + S_end - degree_in_S[p] <= K){
					t_old_R_end = R_end;
					move_u_from_R_to_S(p,S_end,R_end,level);
					if(!remove_vertices_and_prune(S_end,R_end,level))
					BB_search(S_end,R_end,level+1);
					restore_R(S_end,R_end,t_old_R_end,level);
					remove_u_from_S_with_prune(S_end,R_end,level);
				}
			}

			int temp_m_edges = 0;
			std::vector<ui> v_neigh;
			for(auto v:P2)if(SR_rid[v] >= S_end && SR_rid[v] < R_end){
				t_old_R_end = R_end;
				move_u_from_R_to_S(v,S_end,R_end,level);
				temp_m_edges = m_edges + S_end - degree_in_S[v];

				extract_neighbours_of_u_from_vec(v,S_end,R_end,level,P2,v_neigh);
				if(!remove_vertices_and_prune(S_end,R_end,level) && temp_m_edges <= K)
				for(auto w:v_neigh)if(SR_rid[w] >= S_end && SR_rid[w] < R_end){
					if(temp_m_edges + S_end - degree_in_S[w] > K) continue;
					t_old_R_end_2 = R_end;
					move_u_from_R_to_S(w,S_end,R_end,level);
					if(!remove_vertices_and_prune(S_end,R_end,level))
					BB_search(S_end,R_end,level+1);
					restore_R(S_end,R_end,t_old_R_end_2,level);
					remove_u_from_S_with_prune(S_end,R_end,level);
				}
				restore_R(S_end,R_end,t_old_R_end,level);
				remove_u_from_S_with_prune(S_end,R_end,level);
			}
			
			restore_R(S_end,R_end,old_R_end,level);
			return;
		}
#else // 基于kDC的分支


		// std::cerr << "checkmate 6" << std::endl;
		ui D_end = S_end,X_end = S_end;
		#ifdef NEW_BRANCH
		// partition_D_set(S_end,D_end,X_end,R_end,level);
		// calc_ans_on_D_set(S_end,X_end,R_end,level);
		#endif
		if(D_end == R_end){
			return restore_R(S_end,R_end,old_R_end,level);
		}
		// std::cerr << "checkmate 7" << std::endl;
		if(two_stage){
			if(R_end <= best_solution_size + 1) return restore_R(S_end, R_end, old_R_end, level);
		}
		else{
			ui extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(R_end?(sqrt((2*R_end-1)*(2*R_end-1)+8*K)-(2*R_end - 1)) / 2:0 + eps));
			if(R_end + extra_num <= best_solution_size + 1){
				restore_R(S_end, R_end, old_R_end, level);
				return;
			}
		}
        ui u = n; // u is the branching vertex
        ui must_include = 0;
        for(ui i = D_end;i < R_end;i ++){
			if(degree[SR[i]] + 2 >= R_end) {
				if(u == n || must_include !=1)
				u = SR[i];
       			must_include = 1;
       			break;
       		}
		}
		// std::cerr << "checkmate 8" << std::endl;
		if(u == n) u = choose_branch_vertex_based_on_degree(S_end, D_end, R_end, level);
		assert(u != n&&SR_rid[u] >= S_end);
		assert(SR_rid[u] < R_end);
		assert(SR[SR_rid[u]] == u);
        assert(degree[u] + K >= best_solution_size);
	#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			for(ui j = pstart[u];j < pstart_R[u];j ++) assert(SR_rid[edges[j]] < S_end);
			assert(degree_in_S[u] == pstart_R[u] - pstart[u]);
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end) {
				assert(SR_rid[edges[j]] >= S_end);
				++ cnt;
			}
			assert(degree[u] == cnt + degree_in_S[u]);
		}
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
	#endif

        //printf("here3\n");

        // the first branch includes u into S
        ui pre_best_solution_size = best_solution_size, t_old_R_end = R_end;
        move_u_from_R_to_S(u, S_end, R_end, level);
        if(!remove_vertices_and_prune(S_end, R_end, level)){
			BB_search(S_end, R_end, level+1);
		}
        if(must_include == 1){
        	move_u_from_S_to_R(S_end, R_end, level);
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }

	#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			for(ui j = pstart[u];j < pstart_R[u];j ++) assert(SR_rid[edges[j]] < S_end);
			//assert(degree_in_S[u] == pstart_R[u] - pstart[u]);
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end) {
				if(SR_rid[edges[j]] >= S_end) ++ cnt;
			}
			assert(degree[u] == cnt + degree_in_S[u]);
		}
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
	#endif
        // the second branch excludes u from S
        assert(Qv.empty());
		restore_R(S_end, R_end, t_old_R_end, level);
        bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
		// std::cerr << "checkmate 9" << std::endl;

		if(pruned){
			restore_R(S_end, R_end, old_R_end, level);
			return;
		}
        if(best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices(S_end, R_end, level);
        if(!remove_vertices_and_prune(S_end, R_end, level)){
			BB_search(S_end, R_end, level+1);
		}
        restore_R(S_end, R_end, old_R_end, level);

#endif
		return;
    }

	bool FLLAG = 0;

    bool collect_removable_vertices_based_on_degree_in_S(ui S_end, ui R_end, ui level) {
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++){
			++ cnt[S_end - degree_in_S[SR[i]]];
		}
    	assert(R_end > best_solution_size);

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui remaining_vertices_n = best_solution_size - S_end, idx = S_end;
    	assert(missing_edges <= K&&S_end <= best_solution_size);
    	assert(Qv.empty());
    	for(ui i = 0;i <= S_end;i ++) {
    		if(cnt[i] < remaining_vertices_n) {
    			remaining_vertices_n -= cnt[i];
    			missing_edges += i*cnt[i];
    		}
    		else {
    			missing_edges += i*remaining_vertices_n;
    			idx = i;

    			ui next_value = i;
    			if(cnt[i] == remaining_vertices_n) {
    				for(i ++;i <= S_end;i ++) if(cnt[i]) {
    					next_value = i;
    					break;
    				}
    			}
    			if(missing_edges + next_value > K) return true;

    			break;
    		}
    	}
    	for(ui i = S_end;i < R_end;i ++) if(S_end-degree_in_S[SR[i]] > idx&&missing_edges+S_end-degree_in_S[SR[i]] > K) {
    		removed_level[SR[i]] = level;
    		Qv.push(SR[i]);
		}
    	return false;
    }
	
    bool collect_removable_vertices_based_on_degree_in_S_opt(ui S_end, ui R_end, ui level) {
    	if(S_end >= R_end) return false;
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++){
			++ cnt[S_end - degree_in_S[SR[i]]];
		}
    	assert(R_end > best_solution_size);

    	ui *candidates = buf1;
    	for(ui i = 1;i <= S_end;i ++) cnt[i] += cnt[i-1];
    	for(ui i = S_end;i < R_end;i ++) candidates[-- cnt[S_end-degree_in_S[SR[i]]]] = SR[i];
#ifndef NDEBUG
    	for(ui i = 1;i < R_end-S_end;i ++) assert(degree_in_S[candidates[i]] <= degree_in_S[candidates[i-1]]);
#endif

    	ui *ids = buf;
    	ui ids_n = 0;

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui remaining_vertices_n = best_solution_size - S_end;
    	assert(missing_edges <= K&&S_end <= best_solution_size);
    	assert(Qv.empty());
    	if(remaining_vertices_n >= R_end-S_end) return true;
    	ui t_missing_edges = missing_edges;
    	for(ui i = 0;i <= remaining_vertices_n;i ++) t_missing_edges += S_end - degree_in_S[candidates[i]];
    	if(t_missing_edges > K) return true;

    	ui *suffix_sum = buf2;
    	suffix_sum[R_end-S_end] = 0;
    	for(ui i = 0;i < R_end-S_end;i ++) {
    		ui j = R_end-S_end-1-i;
    		suffix_sum[j] = suffix_sum[j+1] + S_end - degree_in_S[candidates[j]];
    	}

    	ui *ids_rid = buf3;
    	for(ui i = S_end;i < R_end;i ++) ids_rid[SR[i]] = R_end;

    	ui n_same_lower = R_end;
    	for(ui i = 0;i < R_end-S_end&&remaining_vertices_n;i ++) {
    		ui u = candidates[i];
    		if(missing_edges + S_end - degree_in_S[u] > K||ids_n + R_end-S_end-i <= remaining_vertices_n) {
    			removed_level[u] = level;
    			Qv.push(u);
    			continue;
    		}

    		t_missing_edges = missing_edges;
    		ui end_idx = i + 1;
    		if(ids_n < remaining_vertices_n) end_idx = i + 1 + remaining_vertices_n - ids_n;
    		t_missing_edges += suffix_sum[i] - suffix_sum[end_idx];
    		if(t_missing_edges > K) {
    			removed_level[u] = level;
    			Qv.push(u);
    			continue;
    		}

    		if(t_missing_edges+R_end-S_end-1+degree_in_S[u]-degree[u] > K) {
    			ui n_neighbors_more = 0, n_neighbors_same_lower = 0, n_neighbors_same_higher = 0;
    			ui t_n_same_lower = 0;
    			if(ids_n >= remaining_vertices_n) {
    				ui cut_off = degree_in_S[ids[remaining_vertices_n-1]];
    				if(n_same_lower >= R_end) {
    					n_same_lower = 0;
    					for(ui j = 0;j < remaining_vertices_n;j ++) if(degree_in_S[ids[j]] == cut_off) ++ n_same_lower;
    				}
    				t_n_same_lower = n_same_lower;
    				for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(removed_level[edges[j]] > level) {
    					if(degree_in_S[edges[j]] > cut_off) ++ n_neighbors_more;
    					if(degree_in_S[edges[j]] == cut_off) {
    						if(ids_rid[edges[j]] < remaining_vertices_n) ++ n_neighbors_same_lower;
    						else ++ n_neighbors_same_higher;
    					}
    				}
    			}
    			else {
    				ui cut_off = degree_in_S[candidates[i+remaining_vertices_n-ids_n]];
    				for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(removed_level[edges[j]] > level) {
    					if(degree_in_S[edges[j]] > cut_off) ++ n_neighbors_more;
    					if(degree_in_S[edges[j]] == cut_off) {
    						if(SR_rid[edges[j]] <= i+remaining_vertices_n-ids_n) ++ n_neighbors_same_lower;
    						else ++ n_neighbors_same_higher;
    					}
    				}
    				for(ui j = i+remaining_vertices_n-ids_n;j > i&&degree_in_S[candidates[j]] == cut_off;j --) ++ t_n_same_lower;
    				for(ui j = ids_n;j > 0&&degree_in_S[ids[j-1]] == cut_off;j --) ++ t_n_same_lower;
    			}

    			t_missing_edges += remaining_vertices_n-t_n_same_lower-n_neighbors_more;
				// if(remaining_vertices_n < t_n_same_lower + n_neighbors_more){
				// 	std::cerr << "n_neighbors_same_lower: " << remaining_vertices_n << " t_n_same_lower: " << t_n_same_lower + n_neighbors_more << std::endl;
				// }
				if(n_neighbors_same_lower < t_n_same_lower){
					ui non_neighbors_n = t_n_same_lower - n_neighbors_same_lower;
					if(non_neighbors_n > n_neighbors_same_higher) t_missing_edges += non_neighbors_n - n_neighbors_same_higher;
				}

    			if(t_missing_edges > K) {
    				removed_level[u] = level;
    				Qv.push(u);
    				continue;
    			}
    		}

    		if(ids_n < remaining_vertices_n) missing_edges += S_end - degree_in_S[u];
    		ids[ids_n++] = u;
    		ids_rid[u] = ids_n-1;
    	}

    	return false;
    }

    void collect_removable_vertices_based_on_vertex_pair(ui S_end, ui R_end, ui level) {
       	assert(S_end >= 1&&Qv.empty());
       	ui u = SR[S_end-1], neighbors_n = 0, non_neighbors_n = 0;
       	get_neighbors_and_non_neighbors_in_R(u, S_end, R_end, level, neighbors_n, non_neighbors_n);
       	for(ui i = 0;i < neighbors_n;i ++) vis[neighbors[i]] = 1;
       	for(ui i = 0;i < non_neighbors_n;i ++) vis[nonneighbors[i]] = 2;

       	ui missing_edges_in_S = compute_missing_edges_in_S(S_end);
       	assert(missing_edges_in_S <= K);

       	for(ui i = 0;i < neighbors_n;i ++) {
       		ui v = neighbors[i];
			if(removed_level[v] < level) continue;
       		ui common_neighbors = 0, exclusive_non_neighbor_u = 0;
       		for(ui j = pstart_R[v];j < pend[v]&&removed_level[edges[j]] >= level;j ++) {
       			if(vis[edges[j]] == 1) ++ common_neighbors;
       			else if(vis[edges[j]] == 2) ++ exclusive_non_neighbor_u;
       		}
       		ui exclusive_non_neighbor_v = neighbors_n - 1 - common_neighbors;
       		ui exclusive_non_neighbors = exclusive_non_neighbor_u + exclusive_non_neighbor_v;
       		ui common_non_neighbors = non_neighbors_n - exclusive_non_neighbor_u;

       		ui UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			UB += mmin(tmp, common_non_neighbors);
       		}
       		if(UB <= best_solution_size) {
       			removed_level[v] = level;
       			vis[v] = 0;
       			Qv.push(v);
       			neighbors[i] = neighbors[-- neighbors_n];
       			-- i;
       		}
       	}

       	++ missing_edges_in_S;
       	assert(missing_edges_in_S <= K||non_neighbors_n == 0);
       	for(ui i = 0;i < non_neighbors_n;i ++) {
   			ui v = nonneighbors[i];
			if(removed_level[v] < level) continue;
   			ui common_neighbors = 0, exclusive_non_neighbor_u = 0;
   			for(ui j = pstart[v];j < pend[v]&&removed_level[edges[j]] >= level;j ++) {
				if(vis[edges[j]] == 1) ++ common_neighbors;
				else if(vis[edges[j]] == 2) ++ exclusive_non_neighbor_u;
			}
			ui exclusive_non_neighbor_v = neighbors_n - common_neighbors;
       		ui exclusive_non_neighbors = exclusive_non_neighbor_u + exclusive_non_neighbor_v;
			ui common_non_neighbors = non_neighbors_n - 1 - exclusive_non_neighbor_u;

       		ui UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			UB += mmin(tmp, common_non_neighbors);
       		}
   			if(UB <= best_solution_size) {
   				removed_level[v] = level;
   				vis[v] = 0;
				// std::cerr << "???" << std::endl;
   				Qv.push(v);
   				nonneighbors[i] = nonneighbors[-- non_neighbors_n];
   				-- i;
   			}
   		}
       	for(ui i = 0;i < neighbors_n;i ++) vis[neighbors[i]] = 0;
       	for(ui i = 0;i < non_neighbors_n;i ++) vis[nonneighbors[i]] = 0;
    }

    bool collect_removable_vertices(ui S_end, ui R_end, ui level) {
		for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K < best_solution_size) return true;

		for(ui i = S_end;i < R_end;i ++) if(removed_level[SR[i]] > level){
			ui v = SR[i];
			if(degree[v] + K < best_solution_size) {
				removed_level[v] = level;
				Qv.push(v);
			}
		}

		return false;
	}

    bool degree_based_prune(ui S_end, ui R_end, ui level) {
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++){
			++ cnt[S_end - degree_in_S[SR[i]]];
		}

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui t_UB = S_end;
    	for(ui i = 0;i <= S_end;i ++) for(ui j = 0;j < cnt[i];j ++) {
    		if(missing_edges + i > K) break;

    		missing_edges += i;
    		++ t_UB;
    	}

		// if(t_UB <= best_solution_size)
		// std::cerr << "deg??" << std::endl;

    	return t_UB <= best_solution_size;
    }

	ui netflow_double_coloring_prune(ui S_end, ui R_end, ui level, ui *color, ui *color_2, ui *color_buf_1, ui *color_buf_2,  ui color_n,ui *tag, ui &reduce_val, ui &max_color, ui *net_loss){
		// std::cerr << "color1" << std::endl;
		ui calc_time = 0;
		ui *color_array = buf2;
		ui *lost_edge = buf3;
		ui *deg = buf5;
		ui *id = buf6;
		static int cs = 0,inner_cs = 0;
		++ cs;
		inner_cs = cs;
		for(ui i = 0;i <= S_end;++i) deg[i] = 0;
		for(ui i = S_end;i < R_end;++i){
			color_2[SR[i]] = color_n - 1;
			color_array[i - S_end] = 0;
			lost_edge[SR[i]] = S_end - degree_in_S[SR[i]];
			tag[SR[i]] = 0;
			++deg[lost_edge[SR[i]]];
			color_buf_1[i - S_end] = 0;
			color_buf_2[i - S_end] = 0;
		}
		for(ui i = 1;i <= S_end;++i) deg[i] += deg[i-1];
		for(ui i = S_end;i < R_end;++i) id[S_end+(--deg[lost_edge[SR[i]]])] = SR[i];
		color_array[R_end - S_end] = 0;
		color_table.clear();
		color_cs.clear();
		for(ui i = 0;i < color_n;++i) color_table[i] = 1;
		for(ui i = S_end,j = S_end;i < R_end;++i){
			ui u = SR[i];
			ui color_p = 0,las_p = 0;
			++calc_time;
			color_array[0] = calc_time;
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;++j){
				if(SR_rid[edges[j]] < R_end){
  					color_array[color_2[edges[j]]-color_n + 1] = calc_time;
				}
			}
			while(color_array[color_p] == calc_time){
				las_p = color_p;
				color_p = color_table[1ll*color_p * color_n + color[u]];
				if(color_table.find(1ll*color_p * color_n + color[u]) == color_table.end())
				color_table[1ll*color_p * color_n + color[u]] = color_p + 1;
				// std::cerr << las_p << " -col- " << color_p << " " << 1ll*color_p * color_n + color[u] << std::endl;
			}
			if(color_table.find(1ll*color_p * color_n + color[u]) == color_table.end())
			color_table[1ll*color_p * color_n + color[u]] = color_p + 1;
			color_table[1ll*las_p * color_n + color[u]] = color_table[1ll*color_p * color_n + color[u]];
			color_2[u] = color_p + color_n - 1;
			max_color = std::max(max_color, color_2[u]+1);
		}
		// std::cerr << "begin..." << std::endl;
		// std::cerr << inner_cs << " " << cs << " " << max_color << " " << color_n << std::endl;
		cs = inner_cs;
		// std::cerr << "color2" << std::endl;
		solver->init(color_n,max_color - color_n,R_end-S_end,color,color_2,lost_edge,SR+S_end,K,best_solution_size - S_end,S_end * level);

		// std::cerr << inner_cs << " " << cs << " " << max_color << " " << color_n << std::endl;
		ui res = solver->run(tag,net_loss,color_buf_1,color_buf_2, reduce_val,K);
		// std::cerr << "color3" << std::endl;
		return res;
	}

    ui coloring_based_prune(ui S_end, ui R_end, ui level) {
		#ifndef COLOR_BASED_PRUNE
		return 0;
		#endif
#ifndef NDEBUG
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif
    	ui *color = neighbors;
		// std::cerr << "color 1" << std::endl;
    	ui color_n = degeneracy_ordering_and_coloring_adj(S_end, R_end, level, color),max_color = color_n;
		ui flag = 0;
		ui *color_tag_old = buf;
		ui *net_loss = buf1;
		ui *tag = buf4;
		ui *color_tag = buf6;
		ui *color_tag_old_2 = buf7;
		ui *color_2 = buf8;
		ui *color_tag_2 = buf9;
		ui *color_buf_1 = buf10;
		ui *color_buf_2 = buf11;
    	// printf("color_n: %u\n", color_n);

    	ui missing_edges_n = compute_missing_edges_in_S(S_end);
    	assert(missing_edges_n <= K);
		// std::cerr << "color 2" << std::endl;
    	if(S_end + color_n + K - missing_edges_n <= best_solution_size) return 2;
		ui net_missing_edges = 0, reduce_val = 0;


		// std::cerr << "color 3" << " " << flow_prune << " " << two_stage << std::endl;
		if(flow_prune){
			net_missing_edges = netflow_double_coloring_prune(S_end, R_end, level, color, color_2, color_buf_1, color_buf_2, color_n, tag, reduce_val, max_color, net_loss);
			// std::cerr << "color pruning.." << std::endl;
			net_missing_edges += reduce_val;
			if(net_missing_edges + missing_edges_n > K) return 2;
			// std::cerr << "finish color prune" << std::endl;

			if(look_ahead_prune){
				ui tag_1 = 0,tag_2 = 0;
				for(ui i = 0;i < color_n; ++i) color_tag_old[i] = color_buf_1[i];
				for(ui i = 0;i < max_color - color_n; ++i)color_tag_old_2[i] = color_buf_2[i];
				for(ui i = S_end;i < R_end;++i){
					if(tag[SR[i]] == 1) ++tag_1,++color_tag_old[color[SR[i]]],++color_tag_old_2[color_2[SR[i]] - color_n];
					if(tag[SR[i]] == 3) ++tag_2;
				}
				for(ui i = S_end;i < R_end;++i)if(removed_level[SR[i]] > level){
					if(tag[SR[i]] != 1 && tag[SR[i]] != 3 && net_missing_edges + missing_edges_n + net_loss[SR[i]] > K + reduce_val){
						removed_level[SR[i]] = level;
						Qv.push(SR[i]);
						flag = 1;
						continue;
					}

					ui tag_deg = 0, tag_deg_1=0,tag_deg_2=0, tag_3 = 0;
					memcpy(color_tag,color_tag_old,sizeof(ui)*color_n);
					memcpy(color_tag_2,color_tag_old_2,sizeof(ui)*(max_color-color_n));
					for(ui j = pstart_R[SR[i]];j < pend[SR[i]]&&removed_level[edges[j]] >= level;++j)if(SR_rid[edges[j]] < R_end){
						if(SR_rid[edges[j]] < R_end && tag[edges[j]] > 0 && tag[edges[j]] < 3){
							if(color_tag[color[edges[j]]])	++tag_deg_1,--color_tag[color[edges[j]]];
							if(color_tag_2[color_2[edges[j]] - color_n])	++tag_deg_2,--color_tag_2[color_2[edges[j]] - color_n];
						}
						else if(SR_rid[edges[j]] < R_end && tag[edges[j]] >= 3) ++ tag_3;
					}
					tag_3 = tag_2 - tag_3;
					tag_deg = std::min(tag_deg_1,tag_deg_2);
					if(net_missing_edges + missing_edges_n + S_end - degree_in_S[SR[i]] + std::max((int)tag_1 - (int)tag_deg,0) + tag_3 > K + reduce_val && removed_level[SR[i]] > level){
						removed_level[SR[i]] = level;
						Qv.push(SR[i]);
						flag = 1;
					}
				}
			}
			else
			{
				for(ui i = S_end;i < R_end;++i)if(tag[SR[i]] != 1){
					ui tag_deg = 0;
					if((net_missing_edges + missing_edges_n + S_end - degree_in_S[SR[i]] > K + reduce_val ||
						net_missing_edges + missing_edges_n + best_solution_size > degree[SR[i]] + K + reduce_val) 
					&& removed_level[SR[i]] > level){
						removed_level[SR[i]] = level;
						Qv.push(SR[i]);
						flag = 1;
					}
				}
			}
			// std::cerr << "finish color prune" << std::endl;
		}

#ifndef NDEBUG
    	// for(ui i = 0;i < n;i ++) assert(!vis[i]);
    	for(ui i = S_end;i < R_end;i ++) vis[color[SR[i]]] ;
    	for(ui i = 0;i < color_n;i ++) assert(vis[i]);
    	for(ui i = S_end;i < R_end;i ++) vis[color[SR[i]]] = 0;
    	//for(ui i = S_end;i < R_end;i ++) if(color[SR[i]] == 16) printf("color of %u is 16\n", SR[i]);
#endif

    	ui *head = nonneighbors;
    	ui *next = buf;
    	for(ui i = 0;i <= S_end;i ++) head[i] = n;
    	for(ui i = S_end;i < R_end;i ++) {
    		ui u = SR[i], non_neighbors_n = S_end - degree_in_S[SR[i]];
    		next[u] = head[non_neighbors_n];
    		head[non_neighbors_n] = u;
    	}

    	ui *color_head = buf1;
    	ui *color_next = buf2;
		ui *checker = buf3;
		ui *color_deg = buf4;
		for(ui i = 0;i <= R_end;i ++) checker[i] = 0;
    	for(ui i = 0;i < color_n;i ++) color_head[i] = n, color_deg[i] = 0;
    	for(ui i = 0;i <= S_end;i ++) for(ui u = head[S_end-i];u != n;u = next[u]) {
    		ui c = color[u];
    		assert(c < color_n);
    		color_next[u] = color_head[c];
    		color_head[c] = u;
    	}

    	for(ui i = 0;i <= K;i ++) head[i] = n;
    	for(ui i = 0;i < color_n;i ++) {
    		ui u = color_head[i];
    		assert(u != n&&color[u] == i);
    		color_head[i] = color_next[u];

    		ui non_neighbors_n = S_end - degree_in_S[u];
    		assert(non_neighbors_n <= K);
    		next[u] = head[non_neighbors_n];
    		head[non_neighbors_n] = u;
			color_tag[i] = 0;
    	}

    	ui min_key = 0, LB = best_solution_size - S_end,M_end = S_end,B_end = S_end;
		std::vector<ui> same_degree_point_vec;
    	while(LB) {
    		while(min_key <= K+1&&head[min_key] == n) ++ min_key,same_degree_point_vec.clear();
    		if(min_key > K+1||missing_edges_n + min_key > K) break;
    		missing_edges_n += min_key;
    		-- LB;
    		ui u = head[min_key];
			checker[u] = 1;
    		assert(u != n);
    		head[min_key] = next[u];
			++ color_deg[color[u]];
    		ui new_u = color_head[color[u]];
			swap_pos(M_end, SR_rid[u]);++M_end;
			same_degree_point_vec.push_back(u);
    		if(new_u == n) continue;

    		color_head[color[u]] = color_next[new_u];
    		ui non_neighbors_n = S_end - degree_in_S[u];
    		ui new_key = S_end - degree_in_S[new_u] + 1 + min_key - non_neighbors_n;
			if(new_key > K) new_key = K+1;
    		next[new_u] = head[new_key];
    		head[new_key] = new_u;
    	}
		// std::cerr << "color 4" << " " << flow_prune << " " << two_stage << std::endl;
		B_end = M_end;
		ui prune_key = min_key;
		if(min_key >= K+1 && LB){
			return 2;
		}
    	while(min_key <= K &&head[min_key] == n) ++ min_key;
		assert(min_key <= K+1);
		if(missing_edges_n + min_key > K){
			return 2;
		}

		#ifdef COLOR_AHEAD_PRUNE

		while(head[prune_key] != n){
			ui u = head[prune_key];
			swap_pos(B_end, SR_rid[u]),++B_end;
			head[prune_key] = next[u];
			same_degree_point_vec.push_back(u);
		}
		for(auto u:same_degree_point_vec) color_tag[color[u]] = 1;
		if(prune_key < K){
			while(head[prune_key+1] != n){
				ui u = head[prune_key+1];
				if(color_tag[color[u]])
				swap_pos(B_end, SR_rid[u]),++B_end;
				head[prune_key+1] = next[u];
			}
		}
		for(auto u:same_degree_point_vec){
			ui next_p = color_head[color[u]];
			color_head[color[u]] = color_next[u];
			while(next_p != n && degree_in_S[next_p] == degree_in_S[u]){
				swap_pos(B_end, SR_rid[next_p]),++B_end;
				next_p = color_next[next_p];
			}
		}
		// std::cerr << "color 6" << " " << flow_prune << " " << two_stage << std::endl;

		for(ui i = M_end; i < R_end; i++) {
			ui u = SR[i],non = M_end - S_end,nei = 0,other = 0;
			for(ui j = pstart_R[u]; j < pend[u] && removed_level[edges[j]] >= level; j++) {
				
				if(SR_rid[edges[j]] >= S_end && SR_rid[edges[j]] < M_end) {
					--non;
				}
				else if(SR_rid[edges[j]] >= M_end && SR_rid[edges[j]] < B_end && edges[j] != u){
					++nei;
				}
			}
			if(missing_edges_n + S_end - degree_in_S[SR[i]] + std::max(nei,non) - nei > K && removed_level[SR[i]] > level){
 				removed_level[SR[i]] = level;
				Qv.push(SR[i]);
				flag = 1;
			}
		}

		// if(S_end == 0)
		// std::cerr << S_end << " " << M_end << " " << B_end << " " << R_end << " " << flag << " " << missing_edges_n << std::endl;

		#endif

    	return flag;
    }

    void get_neighbors_and_non_neighbors_in_R(ui u, ui S_end, ui R_end, ui level, ui &neighbors_n, ui &non_neighbors_n) {
    	neighbors_n = non_neighbors_n = 0;
    	for(ui i = pstart_R[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) vis[edges[i]] = 1;
    	for(ui i = S_end;i < R_end;i ++) if(SR[i] != u) {
    		if(vis[SR[i]]) neighbors[neighbors_n++] = SR[i];
    		else nonneighbors[non_neighbors_n++] = SR[i];
    	}
    	for(ui i = pstart_R[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) vis[edges[i]] = 0;
    }

    void move_u_from_R_to_S(ui u, ui &S_end, ui R_end, ui level) {
    	assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end&&SR[SR_rid[u]] == u);
    	swap_pos(S_end, SR_rid[u]);
    	++ S_end;

    	for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
    		++ degree_in_S[edges[i]];
    	}

        ui missing_edges = compute_missing_edges_in_S(S_end);
        assert(missing_edges <= K);
        assert(Qv.empty());
        for(ui i = S_end;i < R_end;i ++) if(removed_level[SR[i]] > level && S_end - degree_in_S[SR[i]] + missing_edges > K) {
        	removed_level[SR[i]] = level;
        	Qv.push(SR[i]);
        }
    }

    void move_u_from_S_to_R(ui &S_end, ui R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[-- S_end];
		for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++){
			if(SR_rid[edges[i]] < R_end) {
				ui v = edges[i];
				-- degree_in_S[v];
				if(SR_rid[v] >= S_end&&pstart_R[v] > pstart[v]){
					for(int j = pstart_R[v]-1;j >= (int)pstart[v]; --j){
						if(edges[j]==u){
							--pstart_R[v];
							std::swap(edges[j],edges[pstart_R[v]]);
							break;
						}
					}
				}
			}
		}
	}

	//from addtional add of u I should consider to pass tot to the function

    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		removed_level[u] = level;
		
		bool ret = false;
		for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++){
			if(SR_rid[edges[i]] < R_end) {
				ui v = edges[i];
				-- degree_in_S[v];
				if(SR_rid[v] >= S_end&&pstart_R[v] > pstart[v]){
					for(int j = pstart_R[v]-1;j >= (int)pstart[v]; --j){
						if(edges[j]==u){
							--pstart_R[v];
							std::swap(edges[j],edges[pstart_R[v]]);
							break;
						}
					}
				}
				-- degree[v];
				if(degree[v] + K < best_solution_size) {
					if(SR_rid[v] < S_end) ret = true;
					else {
						assert(removed_level[v] > level);
						removed_level[v] = level;
						Qv.push(v);
					}
				}
			}
		}
		return ret;
	}

    bool remove_vertices_and_prune(ui S_end, ui &R_end, ui level) {
		while(!Qv.empty()) {
			ui u = Qv.front(); Qv.pop(); // remove u
			assert(SR[SR_rid[u]] == u);
			assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
			-- R_end;
			swap_pos(SR_rid[u], R_end);

			bool terminate = false;
			for(ui i = pstart[u];i < pend[u] && removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
				ui w = edges[i];
				assert(degree[w] > 0);
				-- degree[w];
				if(degree[w] + K < best_solution_size) {
					if(SR_rid[w] < S_end) terminate = true; // UB
					else if(removed_level[w] > level) { // RR
						removed_level[w] = level;
						Qv.push(w);
					}
				}
			}
			if(terminate) return true;
		}

        return false;
    }

    void restore_R(ui S_end, ui &R_end, ui old_R_end, ui level) {
        while(!Qv.empty()) {
            ui u = Qv.front();
			Qv.pop();
            assert(removed_level[u] == level&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
            removed_level[u] = n;
        }
        while(R_end < old_R_end) { // insert u back into R
            ui u = SR[R_end];
            assert(removed_level[u] == level&&SR_rid[u] == R_end);
            removed_level[u] = n;
			degree[u] = degree_in_S[u] = 0;
            for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) {
            	assert(SR[SR_rid[edges[i]]] == edges[i]);
            	if(SR_rid[edges[i]] < R_end) ++ degree[edges[i]],++degree[u];
				if(SR_rid[edges[i]] < S_end) ++ degree_in_S[u];
            }
            ++ R_end;
        }
		return;
    }
    
    ui compute_missing_edges_in_S(ui S_end,bool flag = 0) {
    	ui res = 0,res2 = 0;
    	for(ui i = 0;i < S_end;i ++) {
    		assert(degree_in_S[SR[i]] < S_end);
    		res += S_end - 1 - degree_in_S[SR[i]];
    	}
		// if(flag){
		// 	for(ui i = 0;i < S_end;i ++) {
		// 		std::cerr << degree_in_S[SR[i]] << " ? " << SR[i] << std::endl;
		// 	}
		// 	std::cerr << std::endl;
		// }
    	return res/2;
    }

    void swap_pos(ui i, ui j) {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex_based_on_degree(ui S_end, ui D_end, ui R_end, ui level) {
    	//printf("start choosing branching vertices\n");
    	ui u = SR[D_end];

		ui neigh_tot = 0,neigh_history = R_end;
		
		for(ui i = D_end+1;i < R_end;i ++) {
			ui v = SR[i];
			#ifdef NEW_BRANCH
			neigh_tot = 0;
			if(degree_in_S[v] == S_end && S_end != 0){
				for(ui j = pstart_R[v];j < pend[v] && removed_level[edges[j]] >= level;++j){
					if(SR_rid[edges[j]] < R_end && degree_in_S[edges[j]] == S_end)
					++neigh_tot;
				}
				if(degree_in_S[u] != S_end || neigh_tot < neigh_history) u=v,neigh_history = neigh_tot;
			}
			else{
				if(degree_in_S[v] < degree_in_S[u] || (degree_in_S[v] == degree_in_S[u] && degree[v] < degree[u])) u=v,neigh_history = neigh_tot;
			}
			#else
			if(degree_in_S[u] == S_end) {
				if(degree_in_S[v] != S_end||degree[v] > degree[u]) u = v;
			}
			else if(degree_in_S[v] != S_end&&(degree[v] > degree[u]||(degree[v] == degree[u]&&degree_in_S[v] < degree_in_S[u]))) u = v;	
			#endif
		}

#ifndef NDEBUG
		ui missing_edges_n = compute_missing_edges_in_S(S_end);
		assert(u != n&&missing_edges_n + S_end - degree_in_S[u] <= K);
		assert(degree[u] + 2 < R_end);
#endif

		return u;
    }

    void print_neighbors(ui u, const ui *pstart, const ui *pend, const ui *edges) {
    	std::vector<ui> neighbors;
    	for(ui i = pstart[u];i < pend[u];i ++) neighbors.push_back(edges[i]);
    	//sort(neighbors.begin(), neighbors.end());
    	printf("neighbors of %u:", u);
    	for(ui i = 0;i < neighbors.size();i ++) printf(" %u(%u)", neighbors[i], removed_level[neighbors[i]]);
    	printf("\n");
    }

    void print_array(const char *str, const ui *array, ui idx_start, ui idx_end, ui l) {
    	for(ui i = 0;i < l;i ++) printf(" ");
    	printf("%s:", str);
    	for(ui i = idx_start;i < idx_end;i ++) printf(" %u", array[i]);
    	printf("\n");
    }
};

#endif

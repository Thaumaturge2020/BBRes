#ifndef _KDEFECTIVE_CLIQUE_BB_MATRIX_
#define _KDEFECTIVE_CLIQUE_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"
#include "boost/dynamic_bitset.hpp"
#include "Dinic_dijkstra.h"
#include <bits/extc++.h>

class kDefectiveClique_BB_matrix {
private:
	ui n;
	ui CNT;
	char *matrix;
	long long matrix_size;
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
	ui flow_prune,look_ahead_prune,color_prune;
	std::vector<ui> correct_S_array;
	std::set<ui> correct_S_set;
public:
	ui branch_tot;
	ui flow_init_time;
	ui flow_run_time;
	ui flow_color_time;
	ui flow_jump_time;
	ui flow_ahead_time;
	ui degen_time;
	ui flow_inside_ahead_time;
	double total_ti;
private:
	DinicSolver *solver;

	ListLinearHeap *calc_heap;

	ui *nonneighbor_1,*nonneighbor_2;
    char *vis;

	__gnu_pbds::gp_hash_table <long long,ui> color_table;
	__gnu_pbds::gp_hash_table <long long,ui> color_cs;

	std::set<ui> my_vec;

	const double eps = 1e-6;

public:
    kDefectiveClique_BB_matrix() {
		total_ti = 0;
    	history_max_n = n = 0;
		CNT = 0;
		two_stage = 0;
		flow_prune = 0;
		look_ahead_prune = 0;
		color_prune = 0;
		branch_tot = 0;
		degen_time = 0;
		flow_inside_ahead_time = 0;
		flow_init_time = 0;
		flow_run_time = 0;
		flow_color_time = 0;
		flow_jump_time = 0;
		flow_ahead_time = 0;
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

		color_table.clear();
		color_cs.clear();
    }

    ~kDefectiveClique_BB_matrix() {
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
		correct_S_array = {0,38,3,4,6,7,10,8,33,17,12,83,78,25,65,87,49,72,26,13,32,45,84,55,102,66,68,56,108,42,23,16,88,71};
		correct_S_set.clear();
		for(auto element:correct_S_array)
		correct_S_set.insert(element);
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

	bool debug_check_correct(ui S_end,ui R_end,ui level){
		ui ans_size = correct_S_set.size();
		for(ui i=0;i<S_end;++i){
			if(correct_S_set.find(SR[i]) == correct_S_set.end()) return 1;
		}
		for(ui i=S_end;i<R_end;++i)if(removed_level[SR[i]] > level){
			if(correct_S_set.find(SR[i]) != correct_S_set.end()) --ans_size;
		}
		if(ans_size != S_end) {std::cerr << S_end << " " << ans_size << " " << std::endl;return 0;}
		return 1;
	}

	void partition_D_set(ui S_end, ui &D_end, ui &X_end, ui R_end, ui level){
		++part_cnt;
		D_end = X_end = S_end;
		bool flag = 0;
		for(ui i = S_end;i < R_end;++i){
			nonneighbor_1[SR[i]] = nonneighbor_2[SR[i]] = max_n + 5;
			if(degree_in_S[SR[i]] == S_end) swap_pos(i,D_end),++D_end;
		}
		int T_end = (D_end == S_end ? R_end : D_end);
		for(ui i = S_end;i < T_end;++i){
			flag = 0;
			char *t_matrix = matrix + SR[i]*n;
			if(nonneighbor_2[SR[i]] == max_n + 6) continue;
			for(ui j = i+1;j < T_end;++j){
				if(!t_matrix[SR[j]]){	
					if(nonneighbor_1[SR[j]] == max_n + 5)
					nonneighbor_1[SR[j]] = SR[i];
					else if(nonneighbor_2[SR[j]] == max_n + 5)
					nonneighbor_2[SR[j]] = SR[i];
					else nonneighbor_2[SR[j]] = max_n + 6;

					if(nonneighbor_1[SR[i]] == max_n+5)
					nonneighbor_1[SR[i]] = SR[j];
					else if(nonneighbor_2[SR[i]] == max_n+5)
					nonneighbor_2[SR[i]] = SR[j];
					else flag = 1;
				}
			}
			if(flag == 1) continue;
			swap_pos(i,X_end);++X_end;
		}
		if(D_end != S_end)
		D_end = X_end;
		return;
	}
	
	void partition_C1_set(ui S_end,ui &D_end,ui R_end,ui level){
		D_end = S_end;
		for(ui i = S_end;i < R_end;++i){
			if(degree_in_S[SR[i]] + 1 <= S_end) swap_pos(i,D_end),++D_end;
		}
		return;
	}

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
		// std::cerr << " check D 0 " << std::endl;
		for(ui i = S_end;i < D_end;++i){
			if(nonneighbor_2[SR[i]] != max_n + 5 && SR_rid[nonneighbor_2[SR[i]]] >= D_end) nonneighbor_2[SR[i]] = max_n + 5;
			if(nonneighbor_1[SR[i]] != max_n + 5 && SR_rid[nonneighbor_1[SR[i]]] >= D_end) nonneighbor_1[SR[i]] = max_n + 5;
		}
		// std::cerr << " check D 0.5 " << std::endl;
		for(ui i = S_end;i < D_end;++i){
			tag[SR[i]] = 0;
			if(nonneighbor_2[SR[i]] != max_n + 5){
				tag[SR[i]] = 3;
				if(nonneighbor_1[SR[i]] != max_n + 5 && degree_in_S[SR[i]] < degree_in_S[nonneighbor_1[SR[i]]]) {
					tag[SR[i]] -= 1;
				}
				if(nonneighbor_2[SR[i]] != max_n + 5 && degree_in_S[SR[i]] < degree_in_S[nonneighbor_2[SR[i]]]) {
					tag[SR[i]] -= 1;
				}
			}
		}
		// std::cerr << " check D 1 " << std::endl;
		calc_heap->init_by_minus_three(D_end - S_end,(S_end + 3)*4,SR + S_end,degree_in_S,tag,S_end);
		ui id,key,now_edge = compute_missing_edges_in_S(S_end);
		// std::cerr << " check D 2 " << std::endl;
		for(ui i = S_end;i < D_end;++i) tag[SR[i]] = 0;
		for(ui i = S_end;i < D_end;++i){
			// std::cerr << " check D 2.5 " << " " << i << std::endl;
			calc_heap->pop_min(id,key);
			swap_pos(i,SR_rid[id]);
			// if(S_end == 12 && level == 36)
			// std::cerr << " check D " << " " << i << " " << id << " " << key << " " << now_edge << " " << calc_heap->get_key(2) << " " << calc_heap->get_key(3) << " " << calc_heap->get_key(17) << " " << calc_heap->get_key(18) << std::endl;
			// std::cerr << " check D 2.66 " << " " << i << " " << nonneighbor_1[id] << " " << nonneighbor_2[id] << std::endl;
			now_edge += key/4;
			tag[id] = 1;
			if(nonneighbor_1[id] < max_n && SR_rid[nonneighbor_1[id]] < D_end && tag[nonneighbor_1[id]] == 0) calc_heap_update(id,nonneighbor_1[id]),nonneighbor_1[id] = max_n + 5;
			if(nonneighbor_2[id] < max_n && SR_rid[nonneighbor_2[id]] < D_end && tag[nonneighbor_2[id]] == 0) calc_heap_update(id,nonneighbor_2[id]),nonneighbor_2[id] = max_n + 5;
			// std::cerr << " check D 2.68 " << " " << i << std::endl;
			if(now_edge <= K){
				if(two_stage){
					if(i >= best_solution_size){
						store_a_kDefectiveClique(i + 1);
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
			// std::cerr << " check D 2.7 " << " " << i << std::endl;
		}
		return;
	}

	// ui choose_vertex_by_neighbours(ui &S_end,ui D_end, ui &R_end, ui level){
	// 	if(S_end == 0) return n;
	// 	ui min_C_u = SR[S_end],max_S_u = SR[0];
	// 	for(ui i = 1; i < S_end; ++ i) if(degree_in_S[SR[i]] < degree_in_S[max_S_u]) max_S_u = SR[i];
	// 	for(ui i = D_end; i < R_end; ++ i){
	// 		if(degree[SR[i]] > degree[min_C_u]) {min_C_u = SR[i];break;}
	// 	}
	// 	if(R_end - degree[min_C_u] - 1 <= S_end - degree_in_S[max_S_u] - 1) return min_C_u;
	// 	return n;
	// }

	void load_graph(ui _n,  const std::vector<std::pair<ui,ui> > &vp,ui _two_stage){
		two_stage = _two_stage;
    	n = _n;
		
		if(matrix == NULL) {
			matrix_size = n*(long long)n;
			matrix = new char[matrix_size];
		}
		
		if(((long long)n)*n > matrix_size) {
			do {
				matrix_size *= 2;
			} while(((long long)n)*n > matrix_size);
			delete[] matrix; matrix = new char[matrix_size];
		}
		
		memset(matrix, 0, sizeof(char)*((long long)n)*n);
		for(ui i = 0;i < n;i ++) degree[i] = 0;
		for(ui i = 0;i < vp.size();i ++) {
			assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
			ui a = vp[i].first, b = vp[i].second;
			assert(!matrix[a*n+b]);
			degree[a] ++;
			degree[b] ++;
			matrix[a*n + b] = matrix[b*n + a] = 1;
		}

		#ifdef FLOW_PRUNE
		flow_prune = 1;
		#ifdef COLOR_AHEAD_PRUNE
		look_ahead_prune = 1;
		#else
		look_ahead_prune = 0;
		#endif
		#else
		flow_prune = 0;
		look_ahead_prune = 0;
		#endif
		#ifdef COLOR_PRUNE
		color_prune = 1;
		#else
		color_prune = 0;
		#endif
		if(solver != NULL) delete solver;solver = NULL;
		solver = new DinicSolver(_n,look_ahead_prune);

        // printf("load graph of size n=%u, m=%u (undirected), density=%.5lf\n", n, m/2, double(m)/n/(n-1));
    }

    void kDefectiveClique(ui _K, ui _UB, std::vector<ui> &kDC) {
		my_vec = {0,72,61,38,37,48,49,18,21,22,17,44,14,51,68,57,25,104,9};
        K = _K;
        UB = _UB;
        if(K == 0) {
        	printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
        	return ;
        }
        best_solution_size = kDC.size();

		memset(vis, 0, sizeof(char)*(n+1));
		memset(degree_in_S, 0, sizeof(ui)*(n+1));
		for(ui i = 0;i < n;i ++) removed_level[i] = n;
		for(ui i = 0;i < n;i ++) SR[i] = SR_rid[i] = i;
		while(!Qv.empty()) Qv.pop();

		zero_degree_vertices.clear();

		for(ui i = 0;i < n;i ++) {
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
			for(ui i = 0;i < n;++i) if(matrix[SR[i]]) ++degree_in_S[SR[i]];
			if(SR_rid[0] < R_end)
			BB_search(1,R_end,1);
		}

		else{
			BB_search(0,R_end,1);
		}

		std::sort(kDC.begin(),kDC.end());
        if(best_solution_size > kDC.size()) {
            kDC.clear();
            for(int i = 0;i < best_solution_size;i ++) kDC.push_back(best_solution[i]);
        }
		return;
    }

private:

	ui debug_node_degree(ui u,ui R_end){
		ui ans = 0;
		char* t_matrix = matrix + u*n;
		for(ui i = 0;i < R_end;++i){
			if(t_matrix[SR[i]]){
				++ans;
			}
		}
		return ans;
	}

    void compute_a_heuristic_solution_and_prune(ui R_end, ui level) {
    	// the following computes the degeneracy ordering and a heuristic solution
#ifndef NDEBUG
    	for(ui i = 0;i < R_end;i ++) assert(degree[SR[i]] + K >= best_solution_size);
    	for(ui i = 0;i < R_end;i ++) assert(!vis[SR[i]]);
#endif
		ui *core = neighbors;
		ui *t_degree = nonneighbors;
		ui total_edges = 0;
		for(ui i = 0;i < R_end;i ++) {
			t_degree[SR[i]] = degree[SR[i]];
			assert(t_degree[SR[i]] < R_end);
			total_edges += degree[SR[i]];
		}

		ui max_core = 0;
		for(ui i = 0;i < R_end;i ++) {
			ui u, min_degree = n;
			for(ui j = 0;j < R_end;j ++) if(!vis[SR[j]]&&t_degree[SR[j]] < min_degree) {
				u = SR[j];
				min_degree = t_degree[u];
			}
			assert(min_degree < n);
			if(t_degree[u] > max_core) max_core = t_degree[u];
			core[u] = max_core;

			long long t_n = R_end - i;

			if(t_n*(t_n-1)/2 <= total_edges/2 + K&&R_end - i > best_solution_size) {
				best_solution_size = R_end - i;
				ui cnt = 0;
				for(ui j = 0;j < R_end;j ++) if(!vis[SR[j]]) best_solution[cnt ++] = SR[j];
				assert(cnt == best_solution_size);
				printf("Degen find a solution of size %u\n", best_solution_size);
			}

			vis[u] = 1;
			if(t_degree[u] == 0) continue;
			char *t_matrix = matrix + u*n;
			for(ui j = 0;j < R_end;j ++) if(!vis[SR[j]]&&t_matrix[SR[j]]) {
				-- t_degree[SR[j]];
				total_edges -= 2;
			}
		}
		for(ui i = 0;i < R_end;i ++) vis[SR[i]] = 0;

		assert(Qv.empty());
		for(ui i = 0;i < R_end;i ++) if(core[SR[i]] + K < best_solution_size) {
			assert(removed_level[SR[i]] > level);
			removed_level[SR[i]] = level;
			Qv.push(SR[i]);
		}
    }

    ui degeneracy_ordering_and_coloring(ui S_end, ui R_end, ui *color) {
		ui *rid = buf;
		ui *peel_sequence = buf1;
		ui *t_degree = color;
		for(ui i = S_end;i < R_end;i ++) {
			t_degree[SR[i]] = degree[SR[i]] - degree_in_S[SR[i]];
#ifndef NDEBUG
			ui d = 0;
			char *t_matrix = matrix + SR[i]*n;
			assert(!vis[SR[i]]);
			for(ui j = S_end;j < R_end;j ++) if(t_matrix[SR[j]]) ++ d;
			assert(d == t_degree[SR[i]]);
#endif
		}

		for(ui i = S_end;i < R_end;i ++) {
			ui u, min_degree = n;
			for(ui j = S_end;j < R_end;j ++) if(!vis[SR[j]]&&t_degree[SR[j]] < min_degree) {
				u = SR[j];
				min_degree = t_degree[u];
			}
			assert(min_degree < n);
			peel_sequence[i-S_end] = u;
			rid[u] = i-S_end;

			vis[u] = 1;
			if(t_degree[u] == 0) continue;
			char *t_matrix = matrix + u*n;
			for(ui j = S_end;j < R_end;j ++) if(!vis[SR[j]]&&t_matrix[SR[j]]) -- t_degree[SR[j]];
		}
		for(ui i = S_end;i < R_end;i ++) vis[SR[i]] = 0;

    	ui max_color = 0;
    	for(ui i = R_end-S_end;i > 0;i --) {
    		ui u = peel_sequence[i-1];
			char *t_matrix = matrix + u*n;
			for(ui j = S_end;j < R_end;j ++) if(rid[SR[j]] >= i&&t_matrix[SR[j]]) vis[color[SR[j]]] = 1;
    		for(ui j = 0;;j ++) if(!vis[j]) {
    			color[u] = j;
    			if(j > max_color) max_color = j;
    			break;
    		}
			for(ui j = S_end;j < R_end;j ++) if(rid[SR[j]] >= i&&t_matrix[SR[j]]) vis[color[SR[j]]] = 0;
    	}

    	return max_color + 1;
    }

    /*ui degeneracy_ordering_and_coloring(ui S_end, ui R_end, ui *color) {
		ui *rid = buf;
		ui *peel_sequence = buf1;
		ui *color_array = buf2;
		ui *t_degree = color;
		for(ui i = 0;i < R_end - S_end;i ++) t_degree[i] = color_array[i] = 0;
		for(ui i = S_end;i < R_end;i ++) {
			++t_degree[degree[SR[i]] - degree_in_S[SR[i]]];
#ifndef NDEBUG
			ui d = 0;
			char *t_matrix = matrix + SR[i]*n;
			assert(!vis[SR[i]]);
			for(ui j = S_end;j < R_end;j ++) if(t_matrix[SR[j]]) ++ d;
			assert(d == t_degree[SR[i]]);
#endif
		}
		for(ui i = 1;i < R_end - S_end;i ++) t_degree[i] += t_degree[i-1];

		for(ui i = S_end;i < R_end;i ++){
			peel_sequence[--t_degree[degree[SR[i]] - degree_in_S[SR[i]]]] = SR[i];
			rid[SR[i]] = t_degree[degree[SR[i]] - degree_in_S[SR[i]]];
		}
// 		for(ui i = S_end;i < R_end;i ++) {
// 			t_degree[SR[i]] = degree[SR[i]] - degree_in_S[SR[i]];
// #ifndef NDEBUG
// 			ui d = 0;
// 			char *t_matrix = matrix + SR[i]*n;
// 			assert(!vis[SR[i]]);
// 			for(ui j = S_end;j < R_end;j ++) if(t_matrix[SR[j]]) ++ d;
// 			assert(d == t_degree[SR[i]]);
// #endif
// 		}

// 		for(ui i = S_end;i < R_end;i ++) {
// 			ui u, min_degree = n;
// 			for(ui j = S_end;j < R_end;j ++) if(!vis[SR[j]]&&t_degree[SR[j]] < min_degree) {
// 				u = SR[j];
// 				min_degree = t_degree[u];
// 			}
// 			assert(min_degree < n);
// 			peel_sequence[i-S_end] = u;
// 			rid[u] = i-S_end;

// 			vis[u] = 1;
// 			if(t_degree[u] == 0) continue;
// 			char *t_matrix = matrix + u*n;
// 			for(ui j = S_end;j < R_end;j ++) if(!vis[SR[j]]&&t_matrix[SR[j]]) -- t_degree[SR[j]];
// 		}
// 		for(ui i = S_end;i < R_end;i ++) vis[SR[i]] = 0;




    	ui max_color = 0,color_time = 0;
    	for(ui i = R_end-S_end;i > 0;i --) {
    		ui u = peel_sequence[i-1];
			char *t_matrix = matrix + u*n;
			++color_time;
			for(ui j = i;j < R_end-S_end;j ++) if(t_matrix[peel_sequence[j]]) color_array[color[peel_sequence[j]]] = color_time;
    		for(ui j = 0;;j ++) if(color_array[j] < color_time) {
    			color[u] = j;
    			if(j > max_color) max_color = j;
    			break;
    		}
    	}

    	return max_color + 1;
    }*/

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
    	assert(S_end > best_solution_size);
		best_solution_size = S_end;
		for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
		ui* tag = buf,lost_edges=0,lost_deg = 0;
		memset(tag,0,sizeof(ui)*(S_end+1));
		// printf("kDefectiveClique begin\n");
		for(ui i=0;i<S_end;++i){
			ui u = SR[i],cnt=0;
			char* t_matrix = matrix+u*n;
			for(ui j=0;j<S_end;++j){
				if(t_matrix[SR[j]]) ++cnt;
				else ++lost_edges;
			}
			--lost_edges;
			printf("%u ",u);
		}
		printf("\n");
		printf("new_size:%u,lost_edges: %u\n",S_end, lost_edges/2);
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
			char* t_matrix = matrix+u*n;
			for(ui j=0;j<S_end;++j){
				if(t_matrix[SR[j]]) ++cnt;
			}
			lost_edges += S_end-cnt-1;
		}
#ifndef NDEBUG
		printf("Find a kDefectiveClique of size: %u\n", best_solution_size);
#endif
    }

	int move_noneighbours_of_u_from_R_to_S(ui u, ui &S_end, ui &R_end, ui level) {
		int tot = 0;
		char *t_matrix = matrix + u*n;
        ui missing_edges = compute_missing_edges_in_S(S_end);
		for(ui i=S_end;i<R_end;++i)if(!t_matrix[SR[i]]){
			if(S_end - degree_in_S[SR[i]] + missing_edges > K) return -tot;
			++tot;
			missing_edges += S_end - degree_in_S[SR[i]];
			move_u_from_R_to_S(SR[i],S_end,R_end,level);
			if(remove_vertices_and_prune(S_end,R_end,level)) return -tot;
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
		char *t_matrix = matrix + u*n;
		for(ui i=S_end;i<R_end;++i){
			if(!t_matrix[SR[i]] && SR[i] != u && removed_level[SR[i]] > level){
				vec.push_back(SR[i]);
			}
		}
		return;
	}

	void extract_non_neighbours_of_u(ui u, ui &S_end, ui &R_end, ui level, std::vector<ui> &vec) {
		char *t_matrix = matrix + u*n;
		for(ui i=0;i<R_end;++i){
			// std::cerr << SR[i] << " " << t_matrix[SR[i]] << " " << removed_level[SR[i]] << " " << level << std::endl;
			if(!t_matrix[SR[i]] && SR[i] != u && removed_level[SR[i]] > level){
				vec.push_back(SR[i]);
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

	ui check_tot = 0;

	ui in_check = 0;
	ui out_check = 0;

	ui CORRECT_FLAG = 0;

    void BB_search(ui S_end, ui R_end, ui level) {
		// std::cerr << "check 1 " << S_end << " ... " << R_end << " " << level << std::endl;
		// for(int i=0;i<S_end;++i){
		// 	printf("%u ",SR[i]);
		// }
		// printf("|");
		// for(int i=S_end;i<R_end;++i){
		// 	printf("%u ",SR[i]);
		// }
		// puts("");
		bool t_flag = false;
		++branch_tot;
		ui cnt = 0;
		
		// for(ui i = 0;i < R_end;i ++){
		// 	char* t_matrix = matrix + SR[i]*n;
		// 	ui cnt = 0;
		// 	for(ui j = 0;j < S_end;j ++){
		// 		if(t_matrix[SR[j]]) ++cnt;
		// 	}
		// 	if(degree_in_S[SR[i]] != cnt){
		// 		std::cerr << "error in degree_in_S " << SR[i] << " " << degree_in_S[SR[i]] << " " << cnt << std::endl;
		// 	}
		// }

#ifndef NDEBUG
    	assert(compute_missing_edges_in_S(S_end) <= K);
    	for(ui i = 0;i < R_end;i ++) assert(degree[SR[i]] + K >= best_solution_size);
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif
		if(two_stage){
			if(best_solution_size >= UB||R_end <= best_solution_size) return ;
			fflush(stdout);
			if(S_end > best_solution_size){
				// std::cerr << "S_end" << std::endl;
				store_a_kDefectiveClique(S_end);
			}
			if(R_end > best_solution_size&&is_kDefectiveClique(R_end)){
				// std::cerr << "R_end" << std::endl;
				store_a_kDefectiveClique(R_end);
			}
			if(R_end <= best_solution_size+1) return ;
		}
		else{
			ui S_deg = compute_missing_edges_in_S(S_end);
			if(best_solution_size >= UB) return;
			ui extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(S_end?(sqrt((2*S_end-1)*(2*S_end-1)+8*(K-S_deg))-(2*S_end - 1)) / 2:0 + eps));
			if(S_end + extra_num > best_solution_size){
				// std::cerr << "S_end" << std::endl;
				store_a_kDefectiveClique(S_end,extra_num);
			}
			extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(R_end?(sqrt((2*R_end-1)*(2*R_end-1)+8*(K-S_deg))-(2*R_end - 1)) / 2:0 + eps));
			if(R_end + extra_num > best_solution_size && is_kDefectiveClique(R_end,extra_num)){
				// std::cerr << "R_end" << std::endl;
				store_a_kDefectiveClique(R_end,extra_num);
			}
			if(R_end + extra_num <= best_solution_size + 1) return ;
		}
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

		ui flag = 0,old_R_end = R_end,loss_edge = compute_missing_edges_in_S(S_end);
		for(ui i = S_end;i < R_end;i ++){
			ui u = SR[i];
			if(degree[u] + 2 >= R_end && S_end - degree_in_S[u] + loss_edge <= K){
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
		if(flag){
			BB_search(S_end,R_end,level+1);
			if(CORRECT_FLAG == 1){
				CORRECT_FLAG = 2;
				// std::cerr << S_end << " case(1) " << R_end << " " << level << std::endl;
			}
			while(flag--) move_u_from_S_to_R(S_end,R_end,level);
			restore_R(S_end,R_end,old_R_end,level);
			return;
		}
        old_R_end = R_end;
        if(S_end == 0) {
        	assert(Qv.empty());
        	compute_a_heuristic_solution_and_prune(R_end, level);
        	if(level == 1) printf("First level kdefectiveclique size: %lu\n", best_solution_size);
        	if(remove_vertices_and_prune(S_end, R_end, level)) {
        		restore_R(S_end, R_end, old_R_end, level);
        		return ;
        	}
        }
        assert(Qv.empty());
        bool terminate = collect_removable_vertices_based_on_degree_in_S_opt(S_end, R_end, level);
        if(terminate||remove_vertices_and_prune(S_end, R_end, level)) {
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }
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
		if(two_stage){
			if(R_end > best_solution_size && is_kDefectiveClique(R_end)){
				store_a_kDefectiveClique(R_end);
			}
			if(R_end <= best_solution_size+1){
				restore_R(S_end, R_end, old_R_end, level);
				return;
			}
		}
		else{
			ui extra_num = std::min((ui)zero_degree_vertices.size(),(ui)floor(R_end?(sqrt((2*R_end-1)*(2*R_end-1)+8*K)-(2*R_end - 1)) / 2:0 + eps));
			if(R_end > best_solution_size && is_kDefectiveClique(R_end,extra_num)){
				store_a_kDefectiveClique(R_end,extra_num);
			}
			if(R_end + extra_num <= best_solution_size+1){
				restore_R(S_end, R_end, old_R_end, level);
				return;
			}
		}
		if(R_end == S_end) return;
		flag = 0;
		while(flag = coloring_based_prune(S_end, R_end, level)){
			// std::cerr << "checkmate 5.5 " << flag << std::endl;
			if(flag == 2){
				restore_R(S_end, R_end, old_R_end, level);
				return;
			}
			else{
				if(remove_vertices_and_prune(S_end, R_end, level)) {
					// std::cerr << " checkmate 5.7 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
					restore_R(S_end, R_end, old_R_end, level);
					return ;
				}
				if(two_stage){
					// std::cerr << " checkmate 5.6 " << S_end << " " << R_end << " " << level << " " << branch_tot << std::endl;
					if(R_end > best_solution_size && is_kDefectiveClique(R_end)){
						store_a_kDefectiveClique(R_end);
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
					if(R_end + extra_num > best_solution_size && is_kDefectiveClique(R_end,extra_num)){
						store_a_kDefectiveClique(R_end);
					}
					if(R_end + extra_num <= best_solution_size+1) {
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
				for(ui i = 0;i < R_end;i ++){
					char* t_matrix = matrix + SR[i]*n;
					ui cnt = 0;
					for(ui j = 0;j < S_end;j ++){
						if(t_matrix[SR[j]]) ++cnt;
					}
					// if(degree_in_S[SR[i]] != cnt){
					// 	std::cerr << "error in plan 2 degree_in_S " << SR[i] << " " << degree_in_S[SR[i]] << " " << cnt << std::endl;
					// }
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
				// for(ui i = 0;i < R_end;i ++){
				// 	char* t_matrix = matrix + SR[i]*n;
				// 	ui cnt = 0;
				// 	for(ui j = 0;j < S_end;j ++){
				// 		if(t_matrix[SR[j]]) ++cnt;
				// 	}
				// 	if(degree_in_S[SR[i]] != cnt){
				// 		std::cerr << "error in plan 2 A degree_in_S " << SR[i] << " " << degree_in_S[SR[i]] << " " << cnt << std::endl;
				// 	}
				// }
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
			for(auto v:P2)if(SR_rid[v] >= S_end && SR_rid[v] < R_end){
				t_old_R_end = R_end;
				move_u_from_R_to_S(v,S_end,R_end,level);
				temp_m_edges = m_edges + S_end - degree_in_S[v];
				if(!remove_vertices_and_prune(S_end,R_end,level) && temp_m_edges <= K)
				for(auto w:P2)if(SR_rid[w] >= S_end && SR_rid[w] < R_end && matrix[v*n+w]){
					if(temp_m_edges + S_end - degree_in_S[w] > K) continue;
					t_old_R_end_2 = R_end;
					
					// for(ui i = 0;i < R_end;i ++){
					// 	char* t_matrix = matrix + SR[i]*n;
					// 	ui cnt = 0;
					// 	for(ui j = 0;j < S_end;j ++){
					// 		if(t_matrix[SR[j]]) ++cnt;
					// 	}
					// }
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
		ui D_end = S_end,X_end = S_end;
		#ifdef NEW_BRANCH
		partition_D_set(S_end,D_end,X_end,R_end,level);
		calc_ans_on_D_set(S_end,X_end,R_end,level);
		#endif
		if(D_end == R_end){
			return restore_R(S_end,R_end,old_R_end,level);
		}
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
	#endif

        //printf("here3\n");

        // the first branch includes u into S
        ui pre_best_solution_size = best_solution_size, t_old_R_end = R_end;
        move_u_from_R_to_S(u, S_end, R_end, level);
        if(!remove_vertices_and_prune(S_end, R_end, level)){
			BB_search(S_end, R_end, level+1);
			if(CORRECT_FLAG == 1){
				CORRECT_FLAG = 2;
			}
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
	#endif
        // the second branch excludes u from S
        assert(Qv.empty());
		restore_R(S_end, R_end, t_old_R_end, level);
        bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
		if(pruned){
			restore_R(S_end, R_end, old_R_end, level);
			return;
		}
		
		// std::cerr << "check 11 " << S_end << " ... " << R_end << " " << level << std::endl;

        if(best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices(S_end, R_end,level);

        if(!remove_vertices_and_prune(S_end, R_end, level)){
			BB_search(S_end, R_end, level+1);
			if(CORRECT_FLAG == 1){
				CORRECT_FLAG = 2;
			}
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
			// std::cerr << S_end << " .. " << degree_in_S[SR[i]] << std::endl;
			++ cnt[S_end - degree_in_S[SR[i]]];
		}
    	assert(R_end > best_solution_size);

    	ui *candidates = buf1;
    	for(ui i = 1;i <= S_end;i ++) cnt[i] += cnt[i-1];
    	for(ui i = S_end;i < R_end;i ++){
			candidates[-- cnt[S_end-degree_in_S[SR[i]]]] = SR[i];
		}
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
    			char *t_matrix = matrix + u*n;
    			if(ids_n >= remaining_vertices_n) {
    				ui cut_off = degree_in_S[ids[remaining_vertices_n-1]];
    				ui same_degree_nonneighbors = 0;
    				for(ui j = 0;j < remaining_vertices_n&&t_missing_edges+remaining_vertices_n-j > K;j ++) if(!t_matrix[ids[j]]) {
    					++ t_missing_edges;
    					if(degree_in_S[ids[j]] == cut_off) ++ same_degree_nonneighbors;
    				}
    				for(ui j = remaining_vertices_n;j < ids_n&&t_missing_edges > K;j ++) {
    					if(degree_in_S[ids[j]] != cut_off||same_degree_nonneighbors == 0||t_missing_edges-same_degree_nonneighbors > K) break;
    					if(t_matrix[ids[j]]) {
    						-- same_degree_nonneighbors;
    						-- t_missing_edges;
    					}
    				}
    				for(ui j = i+1;j < R_end-S_end&&t_missing_edges > K;j ++) {
    					if(degree_in_S[candidates[j]] != cut_off||same_degree_nonneighbors == 0||t_missing_edges-same_degree_nonneighbors > K) break;
    					if(t_matrix[candidates[j]]) {
    						-- same_degree_nonneighbors;
    						-- t_missing_edges;
    					}
    				}
    			}
    			else {
    				ui cut_off = degree_in_S[candidates[i+remaining_vertices_n-ids_n]];
    				ui same_degree_nonneighbors = 0;
    				for(ui j = 0;j < ids_n&&t_missing_edges+remaining_vertices_n-j > K;j ++) if(!t_matrix[ids[j]]) {
    					++ t_missing_edges;
    					if(degree_in_S[ids[j]] == cut_off) ++ same_degree_nonneighbors;
    				}
    				for(ui j = i+1;j <= i+remaining_vertices_n-ids_n&&t_missing_edges+remaining_vertices_n-ids_n-j+i+1>K;j ++) if(!t_matrix[candidates[j]]) {
    					++ t_missing_edges;
    					if(degree_in_S[candidates[j]] == cut_off) ++ same_degree_nonneighbors;
    				}
    				for(ui j = i+1+remaining_vertices_n-ids_n;j < R_end-S_end&&t_missing_edges > K;j ++) {
    					if(degree_in_S[candidates[j]] != cut_off||same_degree_nonneighbors == 0||t_missing_edges-same_degree_nonneighbors > K) break;
    					if(t_matrix[candidates[j]]) {
    						-- same_degree_nonneighbors;
    						-- t_missing_edges;
    					}
    				}
    			}

    			if(t_missing_edges > K) {
    				removed_level[u] = level;
    				Qv.push(u);
    				continue;
    			}
    		}

    		if(ids_n < remaining_vertices_n) missing_edges += S_end - degree_in_S[u];
    		ids[ids_n++] = u;
    	}

    	return false;
    }

    void collect_removable_vertices_based_on_vertex_pair(ui S_end, ui R_end, ui level) {
       	assert(S_end >= 1&&Qv.empty());
       	ui u = SR[S_end-1], neighbors_n = 0, non_neighbors_n = 0;
       	char *t_matrix = matrix + u*n;
       	for(ui i = S_end;i < R_end;i ++) {
       		if(t_matrix[SR[i]]) neighbors[neighbors_n ++] = SR[i];
       		else nonneighbors[non_neighbors_n ++] = SR[i];
       	}

       	ui missing_edges_in_S = compute_missing_edges_in_S(S_end);
       	assert(missing_edges_in_S <= K);

       	for(ui i = 0;i < neighbors_n;i ++) {
       		ui v = neighbors[i];
       		ui common_neighbors = 0, common_non_neighbors = 0;
       		t_matrix = matrix + v*n;
       		for(ui j = 0;j < neighbors_n;j ++) if(t_matrix[neighbors[j]]) ++ common_neighbors;
       		for(ui j = 0;j < non_neighbors_n;j ++) if(!t_matrix[nonneighbors[j]]) ++ common_non_neighbors;
       		ui exclusive_non_neighbors = neighbors_n + non_neighbors_n - 1 - common_neighbors - common_non_neighbors;

       		ui local_UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			local_UB += mmin(tmp, common_non_neighbors);
       		}
       		if(local_UB <= best_solution_size) {
       			removed_level[v] = level;
       			Qv.push(v);
       			neighbors[i] = neighbors[-- neighbors_n];
       			-- i;
       		}
       	}

       	++ missing_edges_in_S;
       	assert(missing_edges_in_S <= K||non_neighbors_n == 0);
       	for(ui i = 0;i < non_neighbors_n;i ++) {
   			ui v = nonneighbors[i];
   			ui common_neighbors = 0, common_non_neighbors = 0;
   			t_matrix = matrix + v*n;
       		for(ui j = 0;j < neighbors_n;j ++) if(t_matrix[neighbors[j]]) ++ common_neighbors;
       		for(ui j = 0;j < non_neighbors_n;j ++) if(j != i&&!t_matrix[nonneighbors[j]]) ++ common_non_neighbors;
       		ui exclusive_non_neighbors = neighbors_n + non_neighbors_n - 1 - common_neighbors - common_non_neighbors;

       		ui local_UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			local_UB += mmin(tmp, common_non_neighbors);
       		}
   			if(local_UB <= best_solution_size) {
   				removed_level[v] = level;
   				Qv.push(v);
   				nonneighbors[i] = nonneighbors[-- non_neighbors_n];
   				-- i;
   			}
   		}
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

    	return t_UB <= best_solution_size;
    }

	ui netflow_double_coloring_prune(ui S_end, ui R_end, ui level, ui *color,ui *color_2, ui *color_buf_1, ui *color_buf_2, ui color_n,ui *tag, ui &reduce_val,ui &max_color, ui *net_loss){
		ui *color_array = buf2;
		ui *lost_edge = buf3;
		ui *deg = buf5;
		ui *id = buf6,calc_time = 0;
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
			char* t_matrix = matrix + u*n;
			color_array[0] = calc_time;
			for(ui j = S_end,v;j < i;++j){
				v = SR[j];
				if(t_matrix[v]){
  					color_array[color_2[v]-color_n + 1] = calc_time;
				}
			}
			
			while(color_array[color_p] == calc_time){
				las_p = color_p;
				color_p = color_table[1ll*color_p * color_n + color[u]];
				if(color_table.find(1ll*color_p * color_n + color[u]) == color_table.end())
				color_table[1ll*color_p * color_n + color[u]] = color_p + 1;
			}
			if(color_table.find(1ll*color_p * color_n + color[u]) == color_table.end())
			color_table[1ll*color_p * color_n + color[u]] = color_p + 1;
			color_table[1ll*las_p * color_n + color[u]] = color_table[1ll*color_p * color_n + color[u]];
			color_2[u] = color_p + color_n - 1;
			max_color = std::max(max_color, color_2[u]+1);
		}
		cs = inner_cs;
		solver->init(color_n,max_color - color_n,R_end-S_end,color,color_2,lost_edge,SR+S_end,K,best_solution_size - S_end,S_end * level);
		ui res = solver->run(tag,net_loss,color_buf_1,color_buf_2,reduce_val,K);
		flow_inside_ahead_time += solver->ahead_time;
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
    	ui color_n = degeneracy_ordering_and_coloring(S_end, R_end, color),max_color = color_n;
		ui flag = 0;
		ui *color_tag_old = buf;
		ui *net_loss = buf1;
		ui *tag = buf4;
		ui *color_tag = buf6,color_deg_2 = 0;
		ui *color_tag_old_2 = buf7;
		ui *color_2 = buf8;
		ui *color_tag_2 = buf9;
		ui *color_buf_1 = buf10;
		ui *color_buf_2 = buf11;

    	ui missing_edges_n = compute_missing_edges_in_S(S_end);
    	assert(missing_edges_n <= K);
    	if(S_end + color_n + K - missing_edges_n <= best_solution_size) return 2;
		ui net_missing_edges = 0, reduce_val = 0;

		if(flow_prune){
			net_missing_edges = netflow_double_coloring_prune(S_end, R_end, level, color, color_2,color_buf_1,color_buf_2, color_n, tag, reduce_val, max_color, net_loss);
			net_missing_edges += reduce_val;
			if(net_missing_edges + missing_edges_n > K) return 2;

			if(look_ahead_prune){
				ui tag_1 = 0;
				for(ui i = 0;i < color_n; ++i) color_tag_old[i] = color_buf_1[i];
				for(ui i = 0;i < max_color - color_n; ++i)color_tag_old_2[i] = color_buf_2[i];
				for(ui i = S_end;i < R_end;++i)if(tag[SR[i]] == 1) ++tag_1,++color_tag_old[color[SR[i]]],++color_tag_old_2[color_2[SR[i]] - color_n];
				for(ui i = S_end;i < R_end;++i)if(removed_level[SR[i]] > level){
					ui tag_deg = 0, tag_deg_2 = 0,tag_deg_1 = 0, tag_3 = 0;
					if(tag[SR[i]] != 1 && tag[SR[i]] != 3 && net_missing_edges + missing_edges_n + net_loss[SR[i]] > K + reduce_val){
						removed_level[SR[i]] = level;
						Qv.push(SR[i]);
						flag = 1;
						continue;
					}
					char *t_matrix = matrix + SR[i]*n;
					memcpy(color_tag,color_tag_old,sizeof(ui)*color_n);
					memcpy(color_tag_2,color_tag_old_2,sizeof(ui)*(max_color-color_n));
					for(ui j = S_end;j < R_end;++j){
						if(t_matrix[SR[j]] && tag[SR[j]] > 0 && tag[SR[j]] < 3){
							if(color_tag[color[SR[j]]])	++tag_deg_1,--color_tag[color[SR[j]]];
							if(color_tag_2[color_2[SR[j]] - color_n])	++tag_deg_2,--color_tag_2[color_2[SR[j]] - color_n];
						}
						else if(!t_matrix[SR[j]] && tag[SR[j]] == 3){
							++tag_3;
						}
					}
					tag_deg = std::min(tag_deg_1,tag_deg_2);
					// tag_deg = tag_deg_1;
					if(net_missing_edges + missing_edges_n + S_end - degree_in_S[SR[i]] + tag_3 + std::max((int)tag_1 - (int)tag_deg,0) > K + reduce_val){
						removed_level[SR[i]] = level;
						Qv.push(SR[i]);
						flag = 1;
					}
				}
			}
		}

		#ifndef COLOR_PRUNE
		return flag;
		#endif

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


		for(ui i = M_end; i < R_end; i++) {
			ui u = SR[i],non = M_end - S_end,nei = 0,other = 0;
			char *t_matrix = matrix + u*n;
			for(ui i = 0;i<color_n;++i) color_tag[i] = 0;
			for(ui j = S_end; j < M_end; j++) {
				if(t_matrix[SR[j]]) {
					--non;
				}
			}
			for(ui j = M_end; j < B_end; j++){
				if(t_matrix[SR[j]] && color_tag[color[SR[j]]] == 0) {
					++nei;
					color_tag[color[SR[j]]] = 1;
				}
			}
			if(missing_edges_n + S_end - degree_in_S[SR[i]] + std::max(nei,non) - nei > K && removed_level[SR[i]] > level){
 				removed_level[SR[i]] = level;
				Qv.push(SR[i]);
				flag = 1;
			}
		}

		#endif

    	return flag;
    }

    bool move_u_from_R_to_S(ui u, ui &S_end, ui R_end, ui level) {
    	assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end&&SR[SR_rid[u]] == u);
    	swap_pos(S_end, SR_rid[u]);
    	++ S_end;

    	char *t_matrix = matrix + u*n;
    	for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) ++ degree_in_S[SR[i]];

        ui missing_edges = compute_missing_edges_in_S(S_end);
        assert(missing_edges <= K);

        assert(Qv.empty());
        for(ui i = S_end;i < R_end;i ++) if(S_end - degree_in_S[SR[i]] + missing_edges > K) {
        	removed_level[SR[i]] = level;
        	Qv.push(SR[i]);
        }

        return false;
    }

    void move_u_from_S_to_R(ui &S_end, ui R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[-- S_end];
		char *t_matrix = matrix + u*n;
		for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) -- degree_in_S[SR[i]];
	}

	//from addtional add of u I should consider to pass tot to the function

    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		removed_level[u] = level;

		bool ret = false;
		char *t_matrix = matrix + u*n;

		for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
			ui v = SR[i];
			-- degree_in_S[v];
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
		return ret;
	}
    bool remove_vertices_and_prune(ui S_end, ui &R_end, ui level) {
		while(!Qv.empty()) {
			ui u = Qv.front(); Qv.pop(); // remove u
			assert(SR[SR_rid[u]] == u);
			assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
			-- R_end;
			swap_pos(SR_rid[u], R_end);

			// printf("remove %u with degree %u\n", u, degree[u]);

			bool terminate = false;
			char *t_matrix = matrix + u*n;
			for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
				ui w = SR[i];
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
            ui u = Qv.front(); Qv.pop();
            assert(removed_level[u] == level&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
            removed_level[u] = n;
        }
        while(R_end < old_R_end) { // insert u back into R
            ui u = SR[R_end];
            assert(removed_level[u] == level&&SR_rid[u] == R_end);
            removed_level[u] = n;
            char *t_matrix = matrix + u*n;
			degree[u] = degree_in_S[u] = 0;
            for(ui i = S_end;i < R_end;i ++) if(t_matrix[SR[i]]) ++ degree[SR[i]],++ degree[u];
            for(ui i = 0;i < S_end;i ++) if(t_matrix[SR[i]]) ++ degree[SR[i]],++ degree[u],++ degree_in_S[u];

            ++ R_end;
            // printf("restore %u with degree %u\n", u, degree[u]);
        }
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
		Timer t_init;
    	ui u = SR[D_end];

		ui neigh_tot = 0,neigh_history = R_end;
		
		for(ui i = D_end+1;i < R_end;i ++) {
			ui v = SR[i];
			#ifdef NEW_BRANCH
			char *t_matrix = matrix+v*n;
			neigh_tot = 0;
			if(degree_in_S[v] == S_end){
				for(ui j = D_end;j < R_end;++j){
					if(j!=i && t_matrix[SR[j]] && degree_in_S[SR[j]] == S_end)
					++neigh_tot;
				}
				if(degree_in_S[u] != S_end || neigh_tot < neigh_history) u=v,neigh_history = neigh_tot;
			}
			else{
				if(degree_in_S[v] > degree_in_S[u] || (degree_in_S[v] == degree_in_S[u] && degree[v] < degree[u])) u=v,neigh_history = neigh_tot;
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
		total_ti += t_init.elapsed();

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

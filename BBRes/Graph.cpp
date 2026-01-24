#include "Graph.h"
#include "kDefectiveClique_BB_bitset.h"
#include "kDefectiveClique_BB_matrix.h"

using namespace std;

Graph::Graph(const char *_dir, const int _K) {
	dir = string(_dir);
	K = _K;

	n = m = 0;

	pstart = nullptr;
	edges = nullptr;

	kDefectiveClique.clear();

	pend = pend_buf = nullptr;
	tri_cnt = edges_pointer = Qe = nullptr;
	deleted = nullptr;
	vp.clear();
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(pend_buf != nullptr) {
		delete[] pend_buf;
		pend_buf = nullptr;
	}
	if(tri_cnt != nullptr) {
		delete[] tri_cnt;
		tri_cnt = nullptr;
	}
	if(edges_pointer != nullptr) {
		delete[] edges_pointer;
		edges_pointer = nullptr;
	}
	if(Qe != nullptr) {
		delete[] Qe;
		Qe = nullptr;
	}
	if(deleted != nullptr) {
		delete[] deleted;
		deleted = nullptr;
	}
}

void Graph::read_graph_binary() {
	printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
	FILE *f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	printf("\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("m not equal sum of degrees\n");
#endif

	fclose(f);

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart == nullptr) pstart = new ept[n+1];
	if(edges == nullptr) edges = new ui[m];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) {
			fread(edges+pstart[i], sizeof(int), degree[i], f);

			// remove self loops and parallel edges
			ui *buff = edges+pstart[i];
			sort(buff, buff+degree[i]);
			ui idx = 0;
			for(ui j = 0;j < degree[i];j ++) {
				if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
				if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;
				buff[idx ++] = buff[j];
			}
			degree[i] = idx;
		}

		pstart[i+1] = pstart[i] + degree[i];
	}

	fclose(f);

	delete[] degree;
}

void Graph::read_graph() {
	printf("# Start reading graph from an edgelist file, Require files \"edges.txt\"\n");
	printf("# Note that this function is not optimized. Reading from a binary file will be faster\n");
	FILE *f = Utility::open_file((dir + string("/edges.txt")).c_str(), "r");

	fscanf(f, "%u%u", &n, &m);
	m *= 2;
	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	for(ui i = 0;i < m/2;i ++) {
		ui a, b;
		fscanf(f, "%u%u", &a, &b);
		if(a >= n || b >= n) {
			printf("!!! Vertex IDs must be between 0 and n-1. Exit !!!\n");
			return ;
		}
		vp.pb(mp(a,b));
		vp.pb(mp(b,a));
	}
	sort(vp.begin(), vp.end());

	if(pstart != nullptr) delete[] pstart;
	pstart = new ept[n+1];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m];

	pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		pstart[i+1] = pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) edges[pstart[i+1] ++] = vp[idx ++].second;
	}

	fclose(f);

	vp_copy = vp;

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

void Graph::output_one_kDefectiveClique() {
	FILE *fout = Utility::open_file("kDefectiveClique.txt", "w");
	fprintf(fout, "%lu\n", kDefectiveClique.size());
	sort(kDefectiveClique.begin(), kDefectiveClique.end());
	for(ui i = 0;i < kDefectiveClique.size();i ++) fprintf(fout, " %u", kDefectiveClique[i]);
	fprintf(fout, "\n");
	fclose(fout);
}

void Graph::verify_kDefectiveClique() {
	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	FILE *fin = Utility::open_file("kDefectiveClique.txt", "r");

	ui max_kDefectiveClique_size = n, kDefectiveClique_n, idx = 0;
	char ok = 1;
	while(fscanf(fin, "%u", &kDefectiveClique_n) == 1) {
		++ idx;
		if(max_kDefectiveClique_size == n) {
			max_kDefectiveClique_size = kDefectiveClique_n;
			printf("max kDefectiveClique size: %u\n", max_kDefectiveClique_size);
		}
		if(kDefectiveClique_n != max_kDefectiveClique_size) printf("!!! WA kDefectiveClique size: %u!\n", kDefectiveClique_n);
		vector<ui> kDefectiveClique;
		for (ui i = 0; i < kDefectiveClique_n; i++) {
			ui tmp;
			fscanf(fin, "%u", &tmp);
			kDefectiveClique.pb(tmp);
		}

		for (ui i = 0; i < kDefectiveClique.size(); i++) {
			if (vis[kDefectiveClique[i]]) {
				printf("WA kDefectiveClique! Duplicate vertex: %u\n", idx);
				ok = 0;
				break;
			}
			vis[kDefectiveClique[i]] = 1;
		}
		ui cnt = 0;
		for(ui i = 0;i < kDefectiveClique.size();i ++) {
			ui d = 0;
			for(ui j = pstart[kDefectiveClique[i]];j < pstart[kDefectiveClique[i]+1];j ++) if(vis[edges[j]]) ++ d;
			assert(d <= kDefectiveClique.size()-1);
			cnt += kDefectiveClique.size()-1-d;
		}
		if(cnt > K*2) {
			ok = 0;
			printf("WA kDefectiveClique! Missing total %u edges\n", cnt/2);
		}
		for(ui i = 0;i < kDefectiveClique.size();i ++) vis[kDefectiveClique[i]] = 0;
	}
	if(ok) printf("Correct kDefectiveClique!\n");
	fclose(fin);

	delete[] vis;
}

void Graph::kDefectiveClique_degen() {
	Timer t;

	kDefectiveClique.clear();

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n-1);

	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);

	delete heap;
	delete[] vis;
	delete[] degree;
	delete[] core;
	delete[] peel_sequence;

	if(kDefectiveClique.size() < UB) printf("\tHeuristic kDefectiveClique Size: %lu, UB: %u, Time: %s (microseconds)\n", kDefectiveClique.size(), UB, Utility::integer_to_string(t.elapsed()).c_str());
	else printf("\tMaximum kDefectiveClique Size(add prune): %lu, Total Time: %s (microseconds)\n", kDefectiveClique.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::kDefectiveClique_exact() {
	Timer t;

	if(K == 0) {
		printf("For k=0, please invoke a maximum clique computation algorithm!\n");
		return ;
	}

	kDefectiveClique.clear();

	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	ui *degree_pre = new ui[n];
	char *vis = new char[n];
	ui *color_matrix = NULL;
	ListLinearHeap *heap = new ListLinearHeap(n, n-1);

	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);

	bool flag = 0;

	if(kDefectiveClique.size() < UB) {
		ui old_size = kDefectiveClique.size();
		ui *out_mapping = new ui[n];
		ui *rid = new ui[n];
		std::vector<ui> zero_point_array,empty_array;
		empty_array.clear();
		core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true, zero_point_array);
		color_shrink_graph(n, m, peel_sequence, out_mapping, rid, pstart, edges, degree,color_matrix, true);
		if(kDefectiveClique.size()>K+1){
			truss_shrink_graph(n, m, peel_sequence, out_mapping, rid, pstart, edges, degree,heap,true);
		}
		ego_degen(n, m, peel_sequence, pstart, edges, degree, rid, vis, heap, true);

		// std::cerr << kDefectiveClique.size() << " ? ? " << old_size << std::endl;

		if(kDefectiveClique.size() > old_size) {
			old_size = kDefectiveClique.size();
			for(ui i = 0;i < kDefectiveClique.size();i ++) {
				assert(kDefectiveClique[i] < n);
				kDefectiveClique[i] = out_mapping[kDefectiveClique[i]];
			}
			// std::cerr << kDefectiveClique.size() << " !!! " << old_size << std::endl;
			if(kDefectiveClique.size()>K+1){
				truss_shrink_graph(n, m, peel_sequence, out_mapping, rid, pstart, edges, degree, heap, true);
			}
		}
		ego_color_shrink_graph(n, m, peel_sequence, out_mapping, rid, pstart, edges, degree,color_matrix, true);

		Timer tt;

		if(UB > kDefectiveClique.size()) {
			ui preprocess_time = 0,preprocess_once;
			ui run_time = 0,run_once;
			kDefectiveClique_BB *kDefectiveClique_solver = new kDefectiveClique_BB();
			{
				vector<ui> ids;
				ui *pstart_new = new ui[n+1];
				ui *edges_new = new ui[m];
				ui *peel_sequence_rid = core;
				for(ui i = 0;i < n;i ++) peel_sequence_rid[peel_sequence[i]] = i;

				memset(vis, 0, sizeof(char)*n);
				
				kDefectiveClique_BB_matrix *kDefectiveClique_solver_m = new kDefectiveClique_BB_matrix();
				kDefectiveClique_solver->allocate_memory(n, K);
				kDefectiveClique_solver_m->allocate_memory(n, K);

				ui search_cnt = 0;
				double min_density = 1, total_density = 0;

				reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);
				Timer t0;

				// std::vector<ui> ans_array = {1,12,19,28,33,37,39,41,45,46,64,66,71,74,80,81,83,92,93,96,99,101,104,107,116,118,119,122,126,128,131,137,140,145,148,149,150,155,165,169,179,185,189,192,194 };

				int ret = 0;

				for(ui i = n;i&&kDefectiveClique.size() < UB;--i) {
					// if(!(i&8191)) std::cerr << i << " " << n << std::endl;
					// std::cerr << i << " " << ret << std::endl;
					ui u = peel_sequence[i-1];

					// auto idx = lower_bound(ans_array.begin(), ans_array.end(), out_mapping[u]);

					// if(idx != ans_array.end() && *idx == out_mapping[u])
					// ++ret;
					// if(ret == ans_array.size())
					// std::cerr << out_mapping[u] << std::endl;
					// else
					// continue;
					// std::cerr << i << " ??? " << u << std::endl;
					// ui test_idx = lower_bound(vp_copy.begin(),vp_copy.end(),mp(0u,110u)) - vp_copy.begin();
					// printf("%d %d %d %d %d\n",i,vp_copy.size(),test_idx,vp_copy[test_idx].first,vp_copy[test_idx].second);
					if(pend[u]-pstart[u]+K < kDefectiveClique.size()||n-i < kDefectiveClique.size()) continue;
					preprocess_once = t0.elapsed();

					#ifdef BITSET_TEST
					extract_bit_subgraph_with_prune(u, peel_sequence_rid, degree, ids, rid,pstart_new,edges_new,vis, pstart, pend, edges);
					#else
					extract_subgraph_with_prune(u, peel_sequence_rid, degree, ids, rid, vp, vis, pstart, pend, edges);
					#endif
					preprocess_time += t0.elapsed()-preprocess_once;

					if(ids.empty() || ids.size() <= kDefectiveClique.size()) continue;

					// if(!(i&511))
					// std::cerr << i << " " << n << " " << ids.size() << " " << vp.size() << " " << kDefectiveClique.size() << std::endl;

					double density = pstart_new[ids.size()]/(double)ids.size()/(ids.size()-1);
					++ search_cnt;
					total_density += density;
					if(density < min_density) min_density = density;

					// if(out_mapping[ids[0]] == 6)
					{
						preprocess_once = t0.elapsed();
						#ifdef BITSET_TEST
						kDefectiveClique_solver->load_graph(ids.size(), pstart_new, pstart_new+1, edges_new, 1, empty_array);
						#else
						kDefectiveClique_solver_m->load_graph(ids.size(), vp ,1);
						#endif
						preprocess_time += t0.elapsed()-preprocess_once;

						run_once = t0.elapsed();
						ui t_old_size = kDefectiveClique.size();
						#ifdef BITSET_TEST
						kDefectiveClique_solver->kDefectiveClique(K, UB, kDefectiveClique);
						#else
						kDefectiveClique_solver_m->kDefectiveClique(K, UB, kDefectiveClique);
						#endif
						run_time += t0.elapsed()-run_once;
						if(kDefectiveClique.size() > t_old_size) {
							// std::cerr << u << " " << i << " " << n << " " << ids.size() << " " << vp.size() << " " << kDefectiveClique.size() << std::endl;
							// std::cerr << i << " " << out_mapping[ids[0]] << " " << n << " " << ids.size() << " " << vp.size() << " " << kDefectiveClique.size() << std::endl;
							for(ui j = 0;j < kDefectiveClique.size();j ++) kDefectiveClique[j] = ids[kDefectiveClique[j]];
						}
					}

					// break;
				}
				// printf("total_ti:%lf\n",kDefectiveClique_solver_m->total_ti/1e6);
				// printf("2-jump branch_tot:%d,time:%d,flow_inside_ahead_time:%d, flow_init_time:%d, flow_run_time:%d, flow_color_time:%d, flow_ahead_time:%d, degen_time:%d, flow_jump_time:%d\n",kDefectiveClique_solver_m->branch_tot, t0.elapsed(), kDefectiveClique_solver_m->flow_inside_ahead_time, kDefectiveClique_solver_m->flow_init_time,kDefectiveClique_solver_m->flow_run_time,kDefectiveClique_solver_m->flow_color_time,kDefectiveClique_solver_m->flow_ahead_time,kDefectiveClique_solver_m->degen_time,kDefectiveClique_solver_m->flow_jump_time);
				delete[] pstart_new;
				delete[] edges_new;
				delete kDefectiveClique_solver_m;
			}
			// std::cerr << "hey!" << " " << kDefectiveClique.size() << std::endl;
			if(kDefectiveClique.size() <= K){
				preprocess_once = t.elapsed();
				kDefectiveClique_solver->load_graph(n, pstart, pstart+1, edges,0,zero_point_array);
				kDefectiveClique_solver->allocate_memory(n, K);
				preprocess_time += t.elapsed()-preprocess_once;
				run_once = t.elapsed();
				kDefectiveClique_solver->kDefectiveClique(K, UB, kDefectiveClique);
				run_time += t.elapsed()-run_once;
			}
			delete kDefectiveClique_solver;
			printf("preprocess_time: %u run_time: %u\n", preprocess_time, run_time);
		}


		if(kDefectiveClique.size() > old_size) {
			old_size = kDefectiveClique.size();
			for(ui i = 0;i < kDefectiveClique.size();i ++) {
				assert(kDefectiveClique[i] < n);
				kDefectiveClique[i] = out_mapping[kDefectiveClique[i]];
			}
		}




		ui kSize = kDefectiveClique.size(),lost_edges = 0;
		vp = vp_copy;
		FILE* f_out = fopen("kDefectiveClique.txt", "w");
		printf("%d %d\n", kSize,old_size);
		for(ui i = 0;i < kSize; ++i){
			ui u = kDefectiveClique[i];
			for(ui j = i+1;j < kSize; ++j){
				ui v = kDefectiveClique[j];
				ui idx = lower_bound(vp.begin(),vp.end(),mp(u,v)) - vp.begin();
				if(idx == vp.size() || vp[idx] != mp(u,v)) ++lost_edges;
			}
			fprintf(f_out,"%u ", u);
			// std::cerr << u << " ";
			printf("%u ", u);
		}
		fprintf(f_out,"\n");
		printf("\n");
		// std::cerr << std::endl;

		if(lost_edges > K) {
			printf("WA kDefectiveClique! Missing total %u edges(exact)\n", lost_edges);
			// std::cerr << "WA kDefectiveClique! Missing total " << lost_edges << " edges(exact)" << std::endl;
			flag = 1;
			FILE* f = fopen("kDefectiveCliqueTime.txt", "w");
			fprintf(f,"-1 %s (ms)\n", Utility::integer_to_string(t.elapsed()/1000ll).c_str());
			// std::cerr << "-1 " << Utility::integer_to_string(t.elapsed()/1000ll).c_str() << " (ms)" << std::endl;
		}
		else{
			printf("AC kDefectiveClique! Missing total %u edges(exact)\n", lost_edges);
		}

		delete[] out_mapping;
		delete[] rid;

		printf("*** Search time: %s\n", Utility::integer_to_string(tt.elapsed()).c_str());
	}

	delete heap;
	delete[] vis;
	delete[] degree;
	delete[] core;
	delete[] peel_sequence;
	delete[] color_matrix;
	printf("\tMaximum kDefectiveClique Size(add prune): %lu, Total Time: %s (microseconds)\n", kDefectiveClique.size(), Utility::integer_to_string(t.elapsed()).c_str());
	if(flag == 0){
		FILE* f = fopen("kDefectiveCliqueTime.txt", "w");
		fprintf(f,"%d %s (ms)\n",kDefectiveClique.size(), Utility::integer_to_string(t.elapsed()/1000ll).c_str());
		std::cerr << kDefectiveClique.size() << " " << Utility::integer_to_string(t.elapsed()/1000ll).c_str() << " (ms)" << std::endl;
	}
}

void Graph::reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges) {
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		for(ui j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[end ++] = edges[j];
	}
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		for(ui j = pstart[u];j < pend[u]&&rid[edges[j]] > rid[u];j ++) {
			ui v = edges[j];
			edges[pend[v] ++] = u;
			assert(pend[v] <= pstart[v+1]);
		}
	}
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend[i] == pstart[i+1]);
#endif
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		while(end < pstart[i+1]&&rid[edges[end]] > rid[i]) ++ end;
	}
}

// degeneracy-based k-defective clique
// return an upper bound of the maximum k-defective clique size
ui Graph::degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output) {
	Timer t;

	ui threshold = (kDefectiveClique.size() > K? kDefectiveClique.size()-K: 0); // all vertices with degree < threshold can be pruned

	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] < threshold) peel_sequence[queue_n ++] = i;
	for(ui i = 0;i < queue_n;i ++) {
		ui u = peel_sequence[i]; degree[u] = 0;
		for(ept j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j]] > 0) {
			if((degree[edges[j]] --) == threshold) peel_sequence[queue_n ++] = edges[j];
		}
	}
	ui UB = n;
	if(queue_n == n) UB = kDefectiveClique.size();

	ept total_edges = 0;
	memset(vis, 0, sizeof(char)*n);
	for(ui i = 0;i < n;i ++) {
		if(degree[i] >= threshold) {
			peel_sequence[queue_n + (new_size ++)] = i;
			total_edges += degree[i];
		}
		else {
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);

	if(new_size != 0) {
		heap->init(new_size, new_size-1, peel_sequence+queue_n, degree);
		ui max_core = 0;
		ui idx = n;
		UB = 0;
		for(ui i = 0;i < new_size;i ++) {
			ui u, key;
			heap->pop_min(u, key);
			if(key > max_core) max_core = key;
			core[u] = max_core;
			peel_sequence[queue_n + i] = u;

			ui t_UB = core[u] + K + 1;
			if(new_size - i < t_UB) t_UB = new_size - i;
			if(t_UB > UB) UB = t_UB;

			if(idx == n&&(total_edges + 2*K)/(new_size-i) >= new_size - i - 1){
				idx = i;
			}
			vis[u] = 1;

			for(ept j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 0) {
				total_edges -= 2;
				heap->decrement(edges[j], 1);
			}
		}

		if(output) printf("*** Degeneracy kDefectiveClique size: %u, max_core: %u, UB: %u, Time: %s (microseconds)\n", new_size-idx, max_core, UB, Utility::integer_to_string(t.elapsed()).c_str());

		if(new_size - idx > kDefectiveClique.size()) {
			kDefectiveClique.clear();
			for(ui i = idx;i < new_size;i ++) kDefectiveClique.pb(peel_sequence[queue_n + i]);
			if(!output) printf("Degen finds a kDefectiveClique of size: %u\n", new_size - idx);
		}
	}

	return UB;
}

void Graph::extract_bit_subgraph_with_prune(ui u, const ui *p_rid, ui *degree, vector<ui> &ids, ui *rid,ui *pstart_new, ui *edges_new, char *exists, ept *pstart, ept *pend, ui *edges) {
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(!exists[i]);
#endif

	ids.clear();
	ids.push_back(u); exists[u] = 1;
	for(ept i = pstart[u];i < pend[u];i ++) {
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]); exists[edges[i]] = 2;
	}
	assert(pend[u] >= pstart[u+1]||p_rid[edges[pend[u]]] < p_rid[u]);
	
	// Utility::print_array("ids1", ids.data(), 0, ids.size(), 0);

	ui *Q = rid;
	ui Q_n = 0;
	for(ui i = 1;i < ids.size();i ++) {
		ui v = ids[i];
		degree[v] = 0;
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] >= p_rid[u];j ++) {
			if(exists[edges[j]]) ++ degree[v];
		}
		if(degree[v]+K < kDefectiveClique.size()) Q[Q_n++] = v;
	}
	for(ui i = 0;i < Q_n;i ++) {
		ui v = Q[i];
		exists[v] = 3;
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] >= p_rid[u];j ++) if(exists[edges[j]] == 2) {
			if(degree[edges[j]]+K == kDefectiveClique.size()) Q[Q_n++] = edges[j];
			-- degree[edges[j]];
		}
	}
	assert(Q_n < ids.size());
	if(ids.size() - Q_n + K <= kDefectiveClique.size()) {
		for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0;
		ids.clear();
		return ;
	}
	
	ui old_size = ids.size();
	for(ui i = 1;i < old_size;i ++) if(exists[ids[i]] == 2) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] >= p_rid[u];j ++) {
			if(!exists[edges[j]]) {
				ids.push_back(edges[j]);
				exists[edges[j]] = 1;
				degree[edges[j]] = 1;
			}
			else ++ degree[edges[j]];
		}
	}

	ui new_size = 1;
	for(ui i = 1;i < old_size;i ++) {
		if(exists[ids[i]] == 3) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);
	for(ui i = old_size;i < ids.size();i ++) {
		if(degree[ids[i]] + K < kDefectiveClique.size()) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	ids.resize(new_size);
	ui edges_tot = 0;
	for(ui i = 0;i < ids.size();i ++) rid[ids[i]] = i,pstart_new[i] = 0;
	pstart_new[ids.size()] = 0;
	for(ui i = 0;i < ids.size();i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pend[v];j ++) if(exists[edges[j]]) {
			assert(rid[v] < ids.size()&&rid[edges[j]] < ids.size());
			++pstart_new[i+1];
			++pstart_new[rid[edges[j]]+1];
			++edges_tot;
		}
	}
	for(ui i = 0;i < ids.size();i ++) pstart_new[i+1] += pstart_new[i];
	for(ui i = 0;i < ids.size();i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pend[v];j ++) if(exists[edges[j]]) {
			assert(rid[v] < ids.size()&&rid[edges[j]] < ids.size());
			edges_new[--pstart_new[i+1]] = rid[edges[j]];
			edges_new[--pstart_new[rid[edges[j]]+1]] = i;
		}
	}
	for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0,pstart_new[i] = pstart_new[i+1];
	pstart_new[ids.size()] = edges_tot << 1;
}

void Graph::extract_subgraph_with_prune(ui u, const ui *p_rid, ui *degree, vector<ui> &ids, ui *rid, vector<pair<ui,ui> > &vp, char *exists, ept *pstart, ept *pend, ui *edges) {
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(!exists[i]);
#endif

	ids.clear(); vp.clear();
	ids.push_back(u); exists[u] = 1;
	for(ept i = pstart[u];i < pend[u];i ++) {
		assert(p_rid[edges[i]] > p_rid[u]);
		// std::cerr << p_rid[edges[i]] << " ?? " << p_rid[u] << std::endl;
		ids.push_back(edges[i]); exists[edges[i]] = 2;
	}
	assert(pend[u] >= pstart[u+1]||p_rid[edges[pend[u]]] < p_rid[u]);
	
	// Utility::print_array("ids1", ids.data(), 0, ids.size(), 0);

	ui *Q = rid;
	ui Q_n = 0;
	// ui S_n = 0,loss_edges = 0;
	for(ui i = 1;i < ids.size();i ++) {
		ui v = ids[i];
		degree[v] = 0;
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] >= p_rid[u];j ++) {
			if(exists[edges[j]]) ++ degree[v];
		}
		if(degree[v]+K < kDefectiveClique.size()) Q[Q_n++] = v;
	}
	for(ui i = 0;i < Q_n;i ++) {
		ui v = Q[i];
		exists[v] = 3;
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] >= p_rid[u];j ++) if(exists[edges[j]] == 2) {
			if(degree[edges[j]]+K == kDefectiveClique.size()) Q[Q_n++] = edges[j];
			-- degree[edges[j]];
		}
	}
	assert(Q_n < ids.size());
	if(ids.size() - Q_n + K <= kDefectiveClique.size()) {
		for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0;
		ids.clear();
		return ;
	}
	
	ui old_size = ids.size();
	for(ui i = 1;i < old_size;i ++) if(exists[ids[i]] == 2) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pstart[v+1]&&p_rid[edges[j]] >= p_rid[u];j ++) {
			if(!exists[edges[j]]) {
				ids.push_back(edges[j]);
				exists[edges[j]] = 1;
				degree[edges[j]] = 1;
			}
			else if(exists[edges[j]] == 1) ++ degree[edges[j]];
		}
		// std::cerr << v << " " << pstart[v] << " " << pstart[v+1] << " " << pstart[v+1] - pstart[v] << std::endl;
	}

	ui new_size = 1;
	for(ui i = 1;i < old_size;i ++) {
		if(exists[ids[i]] == 3) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);
	// std::cerr << loss_edges << std::endl;
	for(ui i = old_size;i < ids.size();i ++) {
		// std::cerr << degree[ids[i]] << " " << K << " ? " << kDefectiveClique.size() + loss_edges << " " << ids.size() << std::endl;
		if(degree[ids[i]] + K < kDefectiveClique.size()) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	ids.resize(new_size);

	for(ui i = 0;i < ids.size();i ++) rid[ids[i]] = i;
	for(ui i = 0;i < ids.size();i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pend[v];j ++) if(exists[edges[j]]) {
			assert(rid[v] < ids.size()&&rid[edges[j]] < ids.size());
			vp.push_back(make_pair(rid[v], rid[edges[j]]));
		}
	}
	for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0;
}

void Graph::extract_normal_subgraph_and_prune(vector<ui> u_array, ui *degree, vector<ui> &ids,ui *rid, ui *pstart_new, ui *edges_new, char *exists, ept *pstart, ept *pend, ui *edges) {
	for(ui i = 0;i < n;i ++) exists[i] = 0;
	ids.clear();
	for(ui i = 0;i < u_array.size();i ++)if(!exists[u_array[i]]){
		ui u = u_array[i];
		ids.push_back(u); exists[u] = 1;
	}

	for(ui i = 0;i < u_array.size();i ++){
		ui u = u_array[i];
		for(ept j = pstart[u];j < pstart[u+1];j ++)if(!exists[edges[j]]){
			ids.push_back(edges[j]); exists[edges[j]] = 2;
		}
	}

	for(ui i = 0;i < ids.size();++i){
		ui u = ids[i];
		degree[u] = 0;
		for(ept j = pstart[u];j < pstart[u+1];j ++) if(exists[edges[j]] != 0) ++ degree[u];
	}

	ui* Q = rid;
	ui Q_n = 0;

	for(ui i = 0;i < ids.size();++i){
		ui u = ids[i];
		if(degree[u] + K <= kDefectiveClique.size()) Q[Q_n++] = u;
	}

	for(ui i = 0;i < Q_n;++i){
		ui u = Q[i];
		exists[u] = 3;
		for(ept j = pstart[u];j < pstart[u+1];j ++) if(exists[edges[j]] == 2){
			if(degree[edges[j]] + K == kDefectiveClique.size()) Q[Q_n++] = edges[j];
			-- degree[edges[j]];
		}
	}

	ui new_size = u_array.size(),old_size = ids.size();
	for(ui i = u_array.size();i < old_size;i ++) {
		if(exists[ids[i]] == 3) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);
	ids.resize(new_size);

	ui edges_tot = 0;
	for(ui i = 0;i < ids.size();i ++) rid[ids[i]] = i,pstart_new[i] = 0;
	pstart_new[ids.size()] = 0;
	for(ui i = 0;i < ids.size();i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pstart[v+1];j ++) if(exists[edges[j]]) {
			++pstart_new[i+1];
			++edges_tot;
		}
	}
	for(ui i = 0;i < ids.size();i ++) pstart_new[i+1] += pstart_new[i];
	for(ui i = 0;i < ids.size();i ++) {
		ui v = ids[i];
		for(ept j = pstart[v];j < pstart[v+1];j ++) if(exists[edges[j]]) {
			edges_new[--pstart_new[i+1]] = rid[edges[j]];
		}
	}
	for(ui i = 0;i < ids.size();i ++) exists[ids[i]] = 0,pstart_new[i] = pstart_new[i+1];
	pstart_new[ids.size()] = edges_tot;
	
	return;
}

void Graph::ego_degen(ui n, ui m, ui *peel_sequence, ept *pstart, ui *edges, ui *degree, ui *rid, char *vis, ListLinearHeap *heap, bool output) {
	Timer t;
	if(pend == nullptr) pend = new ept[n+1];
	ui *pid = new ui[n+1];
	ui *r_pid = new ui[n+1];
	ui *pid_limit = new ui[n+1];
	ui *pid_low = new ui[n+1];
	ui *degree_pre = new ui[n+1];
	ui *edges_pid = new ui[m+1];
	ui *color = new ui[n+1];
	ui *color_degree = new ui[n+1];
	ui *color_matrix = nullptr;
	ui *color_tag = new ui[2*(n+1)];
	ui max_color = 0,old_vs_size = 0;
	orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);

	if(pend_buf == nullptr) pend_buf = new ept[n+1];
	if(edges_pointer == nullptr) edges_pointer = new ui[m];
	ui *pstart_s = pend_buf;
	ui *pend_s = rid;
	ui *edges_s = edges_pointer;
	unsigned long long MY_CNT = 0;

	vector<ui> Q;
	vector<ui> vs;
	vector<ui> vs_rid;
	memset(vis, 0, sizeof(char)*n);
	vs_rid.resize(n);
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		if(pend[u] - pstart[u] < kDefectiveClique.size()) continue;

		vs.clear();
		color_degree[u] = 0;
		for(ui j = pstart[u];j < pend[u];j ++) {
			color_degree[edges[j]] = 0;
			vs.push_back(edges[j]);
			vis[edges[j]] = 1;
			degree[edges[j]] = 0;
		}
		for(ui j = 0;j < vs.size();j ++) for(ui k = pstart[vs[j]];k < pend[vs[j]];k ++){
			if(vis[edges[k]]) {
				++ degree[vs[j]]; ++ degree[edges[k]];
			}
			++ MY_CNT;
		}
		pend_s[vs[0]] = pstart_s[vs[0]] = 0;
		for(ui j = 1;j < vs.size();j ++) pend_s[vs[j]] = pstart_s[vs[j]] = pstart_s[vs[j-1]] + degree[vs[j-1]];
		for(ui j = 0;j < vs.size();j ++) for(ui k = pstart[vs[j]];k < pend[vs[j]];k ++) if(vis[edges[k]]) {
			edges_s[pend_s[vs[j]]++] = edges[k];
			edges_s[pend_s[edges[k]]++] = vs[j];
		}

		ui threshold = (kDefectiveClique.size() > K+1? kDefectiveClique.size()-K-1: 0); // all vertices with degree < threshold can be pruned
		Q.clear();
		for(ui j = 0;j < vs.size();j ++) if(degree[vs[j]] < threshold) {
			Q.push_back(vs[j]);
			vis[vs[j]] = 0;
		}
		for(ui j = 0;j < Q.size();j ++) for(ui k = pstart_s[Q[j]];k < pend_s[Q[j]];k ++) if(vis[edges_s[k]]) {
			if( (degree[edges_s[k]]--) == threshold) {
				Q.push_back(edges_s[k]);
				vis[edges_s[k]] = 0;
			}
		}
		ui cnt = 0, total_edges = 0;
		for(ui j = 0;j < vs.size();j ++) if(vis[vs[j]]) {
			total_edges += degree[vs[j]];
			vs_rid[vs[j]] = cnt;
			vs[cnt++] = vs[j];
		}
		assert(cnt + Q.size() == vs.size());
		vs.resize(cnt);
		if(cnt == 0) continue;

		heap->init(vs.size(), vs.size()-1, vs.data(), degree);
		bool found = false;
		for(ui ii = 0;ii < vs.size();ii ++) {
			ui v, key;
			heap->pop_min(v, key);
			if(found) {
				kDefectiveClique.push_back(v);
				continue;
			}

			if(vs.size()-ii+1 <= kDefectiveClique.size()) break;

			if((total_edges + 2*K)/(vs.size()-ii) >= vs.size() - ii - 1) {
				kDefectiveClique.clear();
				kDefectiveClique.push_back(u);
				kDefectiveClique.push_back(v);
				found = true;
				continue;
			}

			vis[v] = 0;
			for(ept j = pstart_s[v],mat_id;j < pend_s[v];j ++) if(vis[edges_s[j]]) {
				total_edges -= 2;
				heap->decrement(edges_s[j], 1);
			}
		}
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 0;
	}
	for(ui i = 0;i < n;i ++) pend_buf[i] = pend[i];
	for(ui i = 0;i < n;i ++) for(ept j = pstart[i];j < pend[i];j ++) edges[pend_buf[edges[j]]++] = i;
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend_buf[i] == pstart[i+1]);
#endif

	if(output) printf("*** EGo-Degen kDefectiveClique size: %lu, Time: %s (microseconds)\n", kDefectiveClique.size(), Utility::integer_to_string(t.elapsed()).c_str());

	// std::cerr << "MY_CNT: " << MY_CNT << std::endl;

	delete[] pid;
	delete[] r_pid;
	delete[] pid_limit;
	delete[] pid_low;
	delete[] degree_pre;
	delete[] edges_pid;
	delete[] color;
	delete[] color_degree;
	delete[] color_tag;
	if(color_matrix != nullptr)
		delete[] color_matrix;
}

/*




*/

void Graph::color_shrink_graph(ui &n, ui &m, ui *peel_sequence, ui *out_mapping, ui *rid, ept *&pstart, ui *&edges, ui *degree, ui *color_matrix, bool output){

	Timer t;
	
	ui *color = new ui[n+1];
	ui *tag = new ui[2*(n+1)];
	ui max_color = 0;
	
	if(color_matrix != nullptr)
	delete[] color_matrix;
	color_matrix = new ui[n + 1];

	color_graph(n, m, peel_sequence, pstart, pstart+1, edges, color, max_color, tag);
	while(n && remove_and_color_shrink(n, m, out_mapping, peel_sequence, pstart, edges, rid, degree, color, max_color, color_matrix)){
		color_graph(n, m, peel_sequence, pstart, pstart+1, edges, color, max_color, tag);
	}

	
	
	delete[] color;
	delete[] tag;
	delete[] color_matrix;
	color_matrix = nullptr;

	// std::cerr << "After color shrink, n: " << n << " m: " << m/2 << std::endl;

	if(output) printf("*** After color shrink: n = %s, m = %s (undirected), Time: %s (microseconds)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str(),Utility::integer_to_string(t.elapsed()).c_str());
}


bool Graph::ego_color_shrink(ui &n, ui &m, ui *out_mapping, ui *peel_sequence, ept *pstart, ui *edges, ui *rid, ui *degree, ui *color, ui *color_matrix){
	//printf("begin\n");
	ui *Qv = new ui[n+1],head = 0,tail = 0;
	ui *tag = new ui[n+1];
	ui *vis = new ui[n+1];
	ui *pend = new ui[n+1];	
	ui bsz = kDefectiveClique.size();
	memset(vis,0,sizeof(ui)*(n+1));
	memset(Qv,0,sizeof(ui)*(n+1));
	memset(tag,0,sizeof(ui)*(n+1));
	memset(color,0,sizeof(ui)*(n+1));
	memset(color_matrix, 0, sizeof(ui)*(n+1));
	orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1] - pstart[i];
	ept removed_nodes = 0;
	ui ui_c = 0;
	for(ui i = 0;i < n;i ++){
		int color_cnt = 0;
		for(ui j = pstart[i];j < pstart[i+1];j ++) vis[edges[j]] = i+1;
		for(ui j = pstart[i],u,v;j < pstart[i+1];j ++){
			u = edges[j];
			tag[0] = 1;
			for(ui k = pstart[u];k < pend[u];k ++){
				v = edges[k];
				if(vis[v] == i+1){
					tag[color[v]] = 1;
				}
				++ui_c;
			}
			while(tag[color[u]]) ++ color[u];
			if(color_matrix[color[u]] != i+1)
				++color_cnt;
			color_matrix[color[u]] = i+1;
			tag[0] = 0;
			for(ui k = pstart[u];k < pend[u];k ++){
				v = edges[k];
				if(vis[v] == i+1){
					tag[color[v]] = 0;
				}
			}
		}
		if(color_cnt + K + 1 <= bsz){
			Qv[tail++] = i;
			++removed_nodes;
		}
	}

	// std::cerr << "ui_c:" << ui_c << std::endl;

	for(int i=head;i<tail;++i) tag[Qv[i]] = 1;

	while(head < tail){
		ui u = Qv[head++];
		tag[u] = 2;
		for(ept j = pstart[u];j < pstart[u+1];j ++)if(edges[j] < n && tag[edges[j]] != 2){
			ui v = edges[j];
			--degree[v];
			--degree[u];
			if(v < n && tag[v] == 0 && degree[v] + K + 1 == bsz){
				Qv[tail++] = v;
				tag[v] = 1;
				++removed_nodes;
			}
			edges[j] = n + v;
		}
	}

	// std::cerr << "finish one step" << std::endl;

	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		out_mapping[cnt] = out_mapping[i];
		rid[i] = cnt++;
	}
	cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[peel_sequence[i]] > 0) peel_sequence[cnt++] = rid[peel_sequence[i]];

	// std::cerr << "???" << std::endl;

	ui pos = 0; cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		// std::cerr << degree[i] << std::endl;
		ui start = pstart[i];
		pstart[cnt] = pos;
		for(ui j = start;j < pstart[i+1];j ++) if(edges[j] < n && degree[edges[j]] > 0) edges[pos ++] = rid[edges[j]];
		pos = degree[i] + pstart[cnt];
		++ cnt;
	}
	// std::cerr << " ... " << std::endl;
	// printf("*** ego color removed_nodes: %u, n: %u, cnt: %u\n", removed_nodes, n, cnt);
	// fflush(stdout);


	pstart[cnt] = m = pos;
	n = cnt;

	delete[] Qv;
	delete[] tag;
	delete[] vis;
	delete[] pend;

	if(removed_nodes <= n/10) {
		return false;
	}

	return true;
}

void Graph::ego_color_shrink_graph(ui &n, ui &m, ui *peel_sequence, ui *out_mapping, ui *rid, ept *&pstart, ui *&edges, ui *degree, ui *color_matrix, bool output){

	Timer t;
	
	ui *color = new ui[n+1];
	ui *tag = new ui[n+1];
	if(color_matrix != nullptr)
	delete[] color_matrix;
	color_matrix = new ui[n+1];
	while(n && ego_color_shrink(n, m, out_mapping, peel_sequence, pstart, edges, rid, degree, color, color_matrix));
	
	delete[] color;
	delete[] tag;
	delete[] color_matrix;
	color_matrix = nullptr;

	if(output) printf("*** After ego color shrink: n = %s, m = %s (undirected), Time: %s (microseconds)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str(),Utility::integer_to_string(t.elapsed()).c_str());
}

// in_mapping and out_mapping can be the same array
void Graph::core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output, std::vector<ui> &zero_point_array) {
	Timer t;
	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(core[i] + K >= kDefectiveClique.size()) {
		rid[i] = cnt;
		if(in_mapping == nullptr) out_mapping[cnt] = i;
		else out_mapping[cnt] = in_mapping[i];
		++ cnt;
	}

	if(cnt != n) {
		cnt = 0;
		ept pos = 0;
		for(ui i = 0;i < n;i ++) if(core[i] + K >= kDefectiveClique.size()) {
			ept t_start = pstart[i];
			pstart[cnt] = pos;
			for(ept j = t_start;j < pstart[i+1];j ++) if(core[edges[j]] + K >= kDefectiveClique.size()) {
				edges[pos ++] = rid[edges[j]];
			}
			++ cnt;
		}
		pstart[cnt] = pos;

		//printf("%u %u %u %u\n", n, cnt, core[peel_sequence[n-cnt-1]], core[peel_sequence[n-cnt]]);
		assert(core[peel_sequence[n-cnt-1]] == 0||core[peel_sequence[n-cnt-1]] + K < kDefectiveClique.size());
		assert(cnt == 0||core[peel_sequence[n-cnt]] + K >= kDefectiveClique.size());
		for(ui i = 0;i < cnt;i ++) {
			peel_sequence[i] = rid[peel_sequence[n-cnt+i]];
			core[i] = core[out_mapping[i]];
		}

		if(pos > 0&&pos < m/2) {
			ept *pstart_new = new ept[cnt+1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept)*(cnt+1));
			memcpy(edges_new, edges, sizeof(ui)*pos);
			delete[] pstart; pstart = pstart_new;
			delete[] edges; edges = edges_new;
		}

		n = cnt;
		m = pos;
	}


	if(output) printf("*** After core shrink: n = %s, m = %s (undirected), Time: %s (microseconds)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str(),Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::truss_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *out_mapping, ui *rid, ept *pstart, ui *edges, ui *degree, ListLinearHeap *heap, bool output) {
	Timer t;
	if(pend == nullptr) pend = new ept[n+1];
	if(tri_cnt == nullptr) tri_cnt = new ui[m];

	orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);

	// Timer t0;

	// std::cerr << "color_time:" << t0.elapsed() << std::endl;

	oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);

	while(n&&remove_and_shrink_oriented_tri(n, m, out_mapping, peel_sequence, pstart, pend, edges, tri_cnt, rid, degree)) {
		oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);
	}

	if(pend_buf == nullptr) pend_buf = new ept[n+1];
	if(edges_pointer == nullptr) edges_pointer = new ui[m];
	reorganize_oriented_graph(n, tri_cnt, pstart, pend, pend_buf, edges, edges_pointer, rid);

	if(Qe == nullptr) Qe = new ui[m];
	for(ui i = 0;i < n;i ++) {
		pend[i] = pstart[i+1];
		degree[i] = pstart[i+1] - pstart[i];
	}
	if(deleted == nullptr) deleted = new char[m];
	memset(deleted, 0, sizeof(char)*m);
	truss_peeling(Qe, tri_cnt, edges_pointer, deleted, degree, pstart, pend, edges);

	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		assert(degree[i]+K+1 > kDefectiveClique.size());
		out_mapping[cnt] = out_mapping[i];
		rid[i] = cnt++;
	}
	ui t_cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[peel_sequence[i]] > 0) peel_sequence[t_cnt++] = rid[peel_sequence[i]];
	assert(t_cnt == cnt);
	ui pos = 0; cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		ui start = pstart[i];
		pstart[cnt] = pos;
		for(ui j = start;j < pend[i];j ++) if(!deleted[j]) {
			assert(degree[edges[j]] > 0);
			edges[pos++] = rid[edges[j]];
		}
		++ cnt;
	}
	pstart[cnt] = m = pos;
	n = cnt;

	if(output) printf("*** After truss shrink: n = %s, m = %s (undirected), Time: %s (microseconds)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str(),Utility::integer_to_string(t.elapsed()).c_str());

}

// orient graph
void Graph::orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid) {
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	ui TOT_CNT = 0;
	for(ui i = 0;i < n;i ++) {
		ept &end = pend[i] = pstart[i];
		for(ept j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) std::swap(edges[end],edges[j]),++ end;
	}

#ifndef NDEBUG
	long long sum = 0;
	for(int i = 0;i < n;i ++) sum += pend[i] - pstart[i];
	assert(sum*2 == m);
#endif
}

void Graph::color_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges,ui *color, ui &max_color, ui *tag){
	memset(color,0,sizeof(ui)*n);
	memset(tag,0,sizeof(ui)*(2*n));
	max_color = 0;
	for(ui i = 0;i < n;i ++){
		for(ept j = pstart[i];j < pend[i];j ++){
			tag[color[edges[j]]] = 1;
		}
		color[i] = 1;
		while(tag[color[i]]) ++ color[i];
		for(ept j = pstart[i];j < pend[i];j ++){
			tag[color[edges[j]]] = 0;
		}
		max_color = std::max(max_color,color[i]);
	}
	++max_color;
 	return;
}

bool Graph::remove_and_color_shrink(ui &n, ui &m, ui *out_mapping, ui *peel_sequence, ept *pstart, ui *edges, ui *rid, ui *degree, ui *color, ui max_color, ui *color_matrix){
	memset(color_matrix, 0, sizeof(ui)*(max_color + 1));
	//printf("begin\n");
	ui *Qv = new ui[n+1],head = 0,tail = 0;
	ui *tag = new ui[n+1];
	ui bsz = kDefectiveClique.size();
	memset(Qv,0,sizeof(ui)*(n+1));
	memset(tag,0,sizeof(ui)*(n+1));
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1] - pstart[i];
	ept removed_nodes = 0;
	for(ui i = 0;i < n;i ++){
		int color_cnt = 0;
		for(ui j = pstart[i];j < pstart[i+1];j ++){
			if(color_matrix[color[edges[j]]] != i+1)
				++color_cnt;
			color_matrix[color[edges[j]]] = i+1;
		}
		if(color_cnt + K + 1 <= bsz){
			Qv[tail++] = i;
			tag[i] = 1;
			++removed_nodes;
		}
	}

	while(head < tail){
		ui u = Qv[head++];
		tag[u] = 2;
		for(ept j = pstart[u];j < pstart[u+1];j ++)if(edges[j] < n && tag[edges[j]] != 2){
			ui v = edges[j];
			--degree[v];
			--degree[u];
			if(v < n  && tag[v] == 0 && degree[v] + K + 1 == bsz){
				Qv[tail++] = v;
				tag[v] = 1;
				++removed_nodes;
			}
			edges[j] = n + v;
		}
	}

	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		out_mapping[cnt] = out_mapping[i];
		rid[i] = cnt++;
	}
	cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[peel_sequence[i]] > 0) peel_sequence[cnt++] = rid[peel_sequence[i]];

	ui pos = 0; cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		ui start = pstart[i];
		pstart[cnt] = pos;
		for(ui j = start;j < pstart[i+1];j ++) if(edges[j] < n && degree[edges[j]] > 0) edges[pos ++] = rid[edges[j]];
		pos = degree[i] + pstart[cnt];
		++ cnt;
	}
	// printf("*** color removed_nodes: %u, n: %u, cnt: %u\n", removed_nodes, n, cnt);


	pstart[cnt] = m = pos;
	n = cnt;

	delete[] Qv;
	delete[] tag;

	if(removed_nodes <= n/10) {
		return false;
	}

	return true;
}

// orient graph and triangle counting
void Graph::oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj) {
    Timer tri_compute;
	memset(adj, 0, sizeof(ui)*n);
	long long UI_CNT = 0,ui_cnt = 0,u2_cnt = 0,u3_cnt = 0;
	memset(tri_cnt, 0, sizeof(ui)*m);
	for(ui u = 0;u < n;u ++) {
		for(ept j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;
		ui_cnt += pend[u] - pstart[u];
		for(ept j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];
			for(ept k = pstart[v],u1,u2;k < pend[v];k ++){
				++UI_CNT;
				u1 = edges[k],u2 = adj[u1];
				if(u2) {
					++ ui_cnt;
					++ tri_cnt[j];
					++ tri_cnt[k];
					++ tri_cnt[u2-1];
				}
			}
		}

		for(ept j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}
	// std::cerr << "UI_CNT:" << UI_CNT << " " << ui_cnt << " " << u2_cnt << " " << u3_cnt << std::endl;
#ifndef NDEBUG
	//printf("*** Total number of triangles: %s\n", Utility::integer_to_string(cnt).c_str());
#endif
}

bool Graph::remove_and_shrink_oriented_tri(ui &n, ui &m, ui *out_mapping, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *rid, ui *degree) {
	//printf("begin\n");
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1] - pstart[i];
	ept removed_edges = 0;
	for(ui i = 0;i < n;i ++){
		for(ui j = pstart[i];j < pend[i];j ++) if(tri_cnt[j] + K + 2 <= kDefectiveClique.size()) {
			-- degree[i]; -- degree[edges[j]];
			++ removed_edges;
			if(degree[i] > n || degree[edges[j]] > n) {
				fprintf(stderr,"Error: %u %u %u %u %u %u\n", i, edges[j], degree[i], degree[edges[j]], tri_cnt[j],j);
				fflush(stderr);
			}
		}
	}

	printf("*** removed_edges: %u, m: %u\n", removed_edges, m/2);

	if(removed_edges <= m/4) {
		//printf("finish\n");
		return false;
	}

	//printf("here1\n");

	ui cnt = 0;
	for(ui i = 0;i < n;i ++){
		if(degree[i] > 0) {
			out_mapping[cnt] = out_mapping[i];
			rid[i] = cnt++;
		}
	}
	ui t_cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[peel_sequence[i]] > 0) peel_sequence[t_cnt++] = rid[peel_sequence[i]];
	assert(t_cnt == cnt);

	//printf("here2\n");

	ui pos = 0; cnt = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		ui start = pstart[i];
		pstart[cnt] = pos;
		for(ui j = start;j < pend[i];j ++){
			if(tri_cnt[j] + K + 2 > kDefectiveClique.size()){
				edges[pos ++] = rid[edges[j]];
			}
		}
		pend[cnt] = pos;
		assert(pos-pstart[cnt] <= degree[i]);
		pos += degree[i] - (pos-pstart[cnt]);
		++ cnt;
	}

	pstart[cnt] = m = pos;
	n = cnt;

	//printf("finish\n");
	return true;
}

// reorganize the adjacency lists
// and sort each adjacency list to be in increasing order
void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edges_pointer, ui *buf) {
	for(ui i = 0;i < n;i ++) pend2[i] = pend[i];
	for(ui i = 0;i < n;i ++) {
		for(ept j = pstart[i];j < pend[i];j ++) {
			ept &k = pend2[edges[j]];
			edges[k] = i; tri_cnt[k] = tri_cnt[j];
			++ k;
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend2[i] == pstart[i+1]);
#endif

	for(ui i = 0;i < n;i ++) {
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for(ui i = 0;i < n;i ++) {
		for(ept j = pend2[i];j < pstart[i+1];j ++) {
			ept &k = pend[edges[j]];
			edges[k] = i; tri_cnt[k] = tri_cnt[j];
			edges_pointer[k] = j; edges_pointer[j] = k;
			++ k;
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(pend[i] == pend2[i]);
#endif

	ept *ids = pend2;
	for(ui i = 0;i < n;i ++) {
		if(pend[i] == pstart[i]||pend[i] == pstart[i+1]) continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while(j < pend[i]||k < pstart[i+1]) {
			if(k >= pstart[i+1]||(j < pend[i]&&edges[j] < edges[k])) {
				ids[pos] = edges[j];
				buf[pos ++] = edges_pointer[j ++];
			}
			else {
				ids[pos] = edges[k];
				buf[pos ++] = edges_pointer[k ++];
			}
		}
		assert(pos+pstart[i] == pstart[i+1]);
		for(ept j = 0;j < pos;j ++) {
			ui idx = pstart[i]+j, k = buf[j];
			edges[idx] = ids[j];
			tri_cnt[idx] = tri_cnt[k];
			edges_pointer[idx] = k; edges_pointer[k] = idx;
		}
	}
}

void Graph::compact_neighbors(ui u, ui *tri_cnt, ui *edges_pointer, char *deleted, ept *pstart, ept *pend, ui *edges) {
	ui end = pstart[u];
	for(ui i = pstart[u];i < pend[u];i ++) if(!deleted[i]) {
		edges[end] = edges[i];
		tri_cnt[end] = tri_cnt[i];
		edges_pointer[end] = edges_pointer[i];
		edges_pointer[edges_pointer[end]] = end;
		deleted[end] = 0;
		++ end;
	}
	pend[u] = end;
}

char Graph::find(ui u, ui w, ept b, ept e, char *deleted, ept &idx, ui *edges) {
	if(b >= e) return 0;

	while(b+1 < e) {
		idx = b + (e-b)/2;
		if(edges[idx] > w) e = idx;
		else b = idx;
	}

	idx = b;
	if(edges[idx] == w&&!deleted[idx]) return 1;
	return 0;
}

// return the number of peeled edges
void Graph::truss_peeling(ui *Qe, ui *tri_cnt, ui *edges_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges) {
#ifndef NDEBUG
	char *exist = deleted;
	for(ui i = 0;i < n;i ++) {
		assert(pend[i] == pstart[i+1]);
		for(ui j = pstart[i];j < pstart[i+1];j ++) exist[edges[j]] = 1;
		for(ui j = pstart[i]+1;j < pstart[i+1];j ++) assert(edges[j] > edges[j-1]);
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			assert(edges_pointer[edges_pointer[j]] == j);
			assert(tri_cnt[j] == tri_cnt[edges_pointer[j]]);
			assert(edges[j] != i);
			ui cnt = 0, v = edges[j];
			for(ui k = pstart[v];k < pstart[v+1];k ++) if(exist[edges[k]]) ++ cnt;
			assert(cnt == tri_cnt[j]);
		}
		for(ui j = pstart[i];j < pstart[i+1];j ++) exist[edges[j]] = 0;
	}
	memset(deleted,0,sizeof(char)*n);
#endif
	unsigned long long TRI_CNT = 0;
	assert(kDefectiveClique.size() >= K+2);
	ui t_threshold = kDefectiveClique.size()-K-1;
	ept Qe_n = 0;
	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pend[i];j ++) if(tri_cnt[j] < t_threshold&&edges[j] > i) {
		Qe[Qe_n++] = i; Qe[Qe_n++] = edges[j];
	}
	for(ept j = 0;j < Qe_n;j += 2) {
		ui u = Qe[j], v = Qe[j+1], idx;
		find(u, v, pstart[u], pend[u], deleted, idx, edges);
		assert(edges[idx] == v);

		ui tri_n = tri_cnt[idx];
		deleted[idx] = deleted[edges_pointer[idx]] = 1;
		-- degree[u]; -- degree[v];
		if(pend[u]-pstart[u] > degree[u]*2) {
			compact_neighbors(u, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			//printf("degree[u]: %u, %u\n", degree[u], pend[u]-pstart[u]);
			assert(degree[u] == pend[u] - pstart[u]);
		}
		if(pend[v]-pstart[v] > degree[v]*2) {
			compact_neighbors(v, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			assert(degree[v] == pend[v] - pstart[v]);
		}
		if(pend[u]-pstart[u] < pend[v]-pstart[v]) swap(u,v);

		if(pend[u]-pstart[u] > (pend[v]-pstart[v])*2) { // binary search
		// if(false) {
			for(ept k = pstart[v];k < pend[v];k ++){
				++TRI_CNT;
				if(!deleted[k]) {
					if(tri_n&&find(u, edges[k], pstart[u], pend[u], deleted, idx, edges)) {
						-- tri_n;
						-- tri_cnt[edges_pointer[idx]];
						if( (tri_cnt[idx]--) == t_threshold) {
							Qe[Qe_n++] = u; Qe[Qe_n++] = edges[idx];
						}
						-- tri_cnt[edges_pointer[k]];
						if( (tri_cnt[k]--) == t_threshold) {
							Qe[Qe_n++] = v; Qe[Qe_n++] = edges[k];
						}
					}
				}
			}
		// }
		}
		else { // sorted_merge
			ept ii = pstart[u], jj = pstart[v];
			while(ii < pend[u]&&jj < pend[v]) {
				++ TRI_CNT;
				if(edges[ii] == edges[jj]) {
					if(!deleted[ii]&&!deleted[jj]) {
						-- tri_n;
						-- tri_cnt[edges_pointer[ii]];
						if( (tri_cnt[ii]--) == t_threshold) {
							Qe[Qe_n++] = u; Qe[Qe_n++] = edges[ii];
						}
						-- tri_cnt[edges_pointer[jj]];
						if( (tri_cnt[jj]--) == t_threshold) {
							Qe[Qe_n++] = v; Qe[Qe_n++] = edges[jj];
						}
					}

					++ ii;
					++ jj;
				}
				else if(edges[ii] < edges[jj]) ++ ii;
				else ++ jj;
			}
		}
		//if(tri_n != 0) printf("tri_n: %u\n", tri_n);
		assert(tri_n == 0);
	}

	// std::cerr << "TRI_CNT:" << TRI_CNT << " " << Qe_n << std::endl;
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pend[i];j ++) assert(deleted[j]||tri_cnt[j] >= t_threshold);
#endif
}

#ifndef _DINIC_SOLVER_
#define _DINIC_SOLVER_

#include <bits/stdc++.h>

#include "Utility.h"

class DinicSolver{
    ui *head;
    ui *ver;
    ui *nex;
    ui *e;
    int *cost;
    ui *deg;
    int *dep;
    ui *vis;
    ui tot;
    ui *edge_id;
    ui *max_deg;
    ui last_edge;
    std::deque<ui> Q;
    ui n , m, S, T, D, init_n , n1 , n2;
    ui look_ahead_prune;
    const int inf = 0x3f3f3f3f;
    inline void add(int x,int y,int u,int w){
        ver[++tot] = y,nex[tot] = head[x],head[x] = tot;e[tot] = u,cost[tot] = w;
        ver[++tot] = x,nex[tot] = head[y],head[y] = tot;e[tot] = 0,cost[tot] = -w;
        return;
    }
    bool bfs(){
        // std::cerr << "hey" << std::endl;
        for(ui i = 0;i < n;i ++) vis[i] = 0,dep[i] = inf;
        Q.push_back(T);
        vis[T] = 1;dep[T] = 0;
        while(!Q.empty()){
            ui x = Q.front();Q.pop_front();vis[x] = 0;
            // std::cerr << "X:" << x << std::endl;
            for(ui i = head[x];i;i = nex[i]){
                // std::cerr << i << std::endl;
                ui y = ver[i];
                // std::cerr << "i:" << i << " y:" << y << " " << e[i^1] << " " << dep[y] << " " << dep[x] << " " << cost[i] << std::endl;
                if(e[i^1] && dep[y] > dep[x] - cost[i]){
                    dep[y] = dep[x] - cost[i];
                    if(!vis[y]){
                        vis[y] = 1;
                        // std::cerr << L << " " << buf[L] << std::endl;
                        if(!Q.empty() && dep[Q.front()] < dep[y]) Q.push_front(y);
                        else
                        Q.push_back(y);
                    }
                }
                // std::cerr << "end!" << Q.size() << std::endl;
            }
        }
        // std::cerr << " endl?" << std::endl;
        return dep[S] < inf;
    }

    void spfa(int x){
        // std::cerr << "hey" << std::endl;
        for(ui i = 0;i < n;i ++) vis[i] = 0,dep[i] = inf;
        Q.push_back(x);
        vis[x] = 1;dep[x] = 0;
        while(!Q.empty()){
            ui x = Q.front();Q.pop_front();vis[x] = 0;
            for(ui i = head[x];i;i = nex[i]){
                // std::cerr << i << std::endl;
                ui y = ver[i];
                // std::cerr << "i:" << i << " y:" << y << " " << e[i^1] << " " << dep[y] << " " << dep[x] << " " << cost[i] << std::endl;
                if(e[i^1] && dep[y] > dep[x] - cost[i]){
                    dep[y] = dep[x] - cost[i];
                    if(!vis[y]){
                        vis[y] = 1;
                        // std::cerr << L << " " << buf[L] << std::endl;
                        if(!Q.empty() && dep[Q.front()] < dep[y]) Q.push_front(y);
                        else
                        Q.push_back(y);
                    }
                }
                // std::cerr << "end!" << Q.size() << std::endl;
            }
        }
        // std::cerr << " endl?" << std::endl;
        return;
    }
    int dfs(ui x,ui F){
        if(x == T || F == 0) return F;
        vis[x] = 1;
        int f,flow = 0;
        for(ui i = head[x];i;i = nex[i]){
            ui y = ver[i];
            if(!vis[y] && dep[y] == dep[x] - cost[i] && e[i]){
                f = dfs(y,std::min(F,e[i]));
                e[i] -= f;
                e[i^1] += f;
                flow += f;
                F -= f;
                if(F == 0) return flow;
            }
        }
        vis[x] = 0;
        return flow;
    }
    public:
    DinicSolver(ui _n,ui _look_ahead_prune){
        init_n = _n;
        head = new ui[6*init_n+5];
        ver = new ui[6*init_n+5];
        nex = new ui[6*init_n+5];
        e = new ui[6*init_n+5];
        cost = new int[6*init_n+5];
        deg = new ui[6*init_n+5];
        dep = new int[6*init_n+5];
        vis = new ui[6*init_n+5];
        edge_id = new ui[6*init_n+5];
        max_deg = new ui[6*init_n+5];
        tot = 0;
        look_ahead_prune = _look_ahead_prune;
        Q.clear();
    }
    ~DinicSolver(){
        if(head != nullptr){
            delete[] head;
            head = nullptr;
        }
        if(ver != nullptr){
            delete[] ver;
            ver = nullptr;
        }
        if(nex != nullptr){
            delete[] nex;
            nex = nullptr;
        }
        if(e != nullptr){
            delete[] e;
            e = nullptr;
        }
        if(cost != nullptr){
            delete[] cost;
            cost = nullptr;
        }
        if(deg != nullptr){
            delete[] deg;
            deg = nullptr;
        }
        if(dep != nullptr){
            delete[] dep;
            dep = nullptr;
        }
        if(vis != nullptr){
            delete[] vis;
            vis = nullptr;
        }
        if(edge_id != nullptr){
            delete[] edge_id;
            edge_id = nullptr;
        }
        if(max_deg != nullptr){
            delete[] max_deg;
            max_deg = nullptr;
        }
    }
    void init(ui _n1,ui _n2,ui _m,ui *from,ui *to,ui *val,ui *id, ui K,ui lim){
        tot = 1;
        // std::cerr << "n:" << n << std::endl;
        // std::cerr << _n1 << " n1+n2 " << _n2 << std::endl;
        n = _n1+_n2+3;m = _m;
        n1 = _n1;n2 = _n2;
        // std::cerr << "n:" << n << " " << _n1 << " " << _n2 << " " << m << std::endl;
        S = n-1,T = n-2,D = n-3;
        memset(head,0,sizeof(ui)*n);
        memset(deg,0,sizeof(ui)*n);
        memset(max_deg,0,sizeof(ui)*n);
        for(ui i = 0;i < m;++ i){
            ui u = from[id[i]],v = to[id[i]];
            // std::cerr << u << " -u-v- " << v << " " << n << std::endl;
            add(u,v,1,val[id[i]]);
            edge_id[tot] = edge_id[tot^1] = id[i];
            ++deg[u];
            ++deg[v];
            // std::cerr << u << " ?  " << v << " " << val[i] << " " << m << " " << n << std::endl;
        }
        for(ui i = 0;i < _n1;i ++){
            for(ui j = 0,s = 0;j < deg[i] && s < K;s+=j, j ++){
                add(D,i,1,j);
                // std::cerr << D << " ! " << i << " " << j << " " << deg[i] << std::endl;
            }
        }
        for(ui i = _n1;i < _n1 + _n2;i ++){
            for(ui j = 0,s = 0;j < deg[i] && s < K;s+=j, j ++){
                add(i,T,1,j);
                // std::cerr << i << " .. " << T << " " << j << std::endl;
            }
        }
        add(S,D,lim,0);
        last_edge = tot-1;
        // std::cerr << _n1 << " -n1-n2- " << _n2 << std::endl;
        // std::cerr << tot << std::endl;
        return;
    }
    int run(ui *rid, ui &reduce_val, ui K){
        int fans = 0;
        while(bfs()){
            vis[T] = 1;
            while(vis[T]){
                memset(vis,0,sizeof(ui)*n);
                ui f = dfs(S,inf);
                fans += f*dep[S];
                if(f) reduce_val = dep[S];
                if(fans > K) break;
                // std::cerr << "f:" << f << " " << dep[S] << " " << fans << std::endl;
            }
            if(fans > K) break;
        }

        for(ui i = 0;i < n1;++i){
            for(ui j = head[i];j;j = nex[j]){
                if(ver[j] != D && e[j] == 0){
                    rid[edge_id[j]] = 1;
                }
            }
        }

        if(look_ahead_prune){
            for(ui i = 0;i < n1;++i){
                spfa(i);
                for(ui j = head[i];j;j = nex[j]){
                    if(ver[j] != D && e[j] == 1){
                        if(dep[ver[j]] + cost[j] > 0){
                            rid[edge_id[j]] = 2;
                        }
                    }
                }
            }
        }
        // if(fans <= K){
        //     ++e[last_edge];
        //     if(bfs()){
        //         memset(vis,0,sizeof(ui)*n);
        //         ui f = dfs(S,inf);
        //         fans += f*dep[S];
        //         if(f) reduce_val = dep[S];
        //     }
        // }

        return fans;
    }
};

#endif
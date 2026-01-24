#ifndef _DINIC_SOLVER_
#define _DINIC_SOLVER_

#include <bits/stdc++.h>

#include "Utility.h"
#define fi first
#define se second
#define pii std::pair<int,int>

class DinicSolver{
    public:
    ui ahead_time;
    private:
    ui *head,*cur;
    ui *ver;
    ui *nex;
    ui *e,*dep,*dep2,*dep3;
    int *cost;
    ui *deg;
    int *h;
    ui *vis;
    ui tot,edge_tot;
    ui *edge_id;
    ui *max_deg;
    ui check_ID;
    ui *dfn,*low,*in,*id,*stk,top,timestamp,scc_num;
    ui last_edge;
    std::priority_queue<pii> V;
    ListLinearHeap *Q;
    std::deque<ui> P;
    ui n , m, S, T, D, init_n , n1 , n2;
    ui look_ahead_prune,flow_lim;
    ui INIT_LIMIT;
    ui start,target;
    const int inf = 0x3f3f3f3f;
    inline void add(int x,int y,int u,int w){
        ver[++tot] = y,nex[tot] = head[x],head[x] = tot;e[tot] = u,cost[tot] = w;
        ver[++tot] = x,nex[tot] = head[y],head[y] = tot;e[tot] = 0,cost[tot] = -w;
        return;
    }
    bool dijk(ui *dep,ui flag = 1,bool update = 1){
        if(flag) start = S,target = T;
        else     start = T,target = S;
        for(ui i = 0;i < n;i ++) vis[i] = 0,dep[i] = INIT_LIMIT;
        dep[start] = 0;
        Q->init(n,INIT_LIMIT,nullptr,dep);
        ui x,key;
        while(Q->pop_min(x,key)){
            if(vis[x]) continue;
            vis[x] = 1;
            for(ui i = head[x];i;i = nex[i]){
                ui y = ver[i];
                if(e[i] && dep[y] > dep[x] + h[x] + cost[i] - h[y]){
                    Q->decrement(y, dep[y] - (dep[x] + h[x] + cost[i] - h[y]));
                    dep[y] = dep[x] + h[x] + cost[i] - h[y];
                }
            }
        }
        if(update)
        for(ui i = 0;i < n;i ++) h[i] += dep[i];
        return dep[target] < INIT_LIMIT;
    }

    void inv_dijk(ui *dep,ui flag = 1){
        if(flag) start = T,target = S;
        else     start = D,target = T;
        for(ui i = 0;i < n;i ++) vis[i] = 0,dep[i] = INIT_LIMIT;
        dep[start] = 0;
        Q->init(n,INIT_LIMIT,nullptr,dep);
        ui x,key,bfs_cnt = 0;
        while(Q->pop_min(x,key)){
            if(vis[x]) continue;
            vis[x] = 1;
            for(ui i = head[x];i;i = nex[i]){
                ui y = ver[i];
                if(e[i^1] && dep[y] > dep[x] - (h[x] + cost[i] - h[y])){
                    Q->decrement(y, dep[y] - (dep[x] - (h[x] + cost[i] - h[y])));
                    dep[y] = dep[x] - (h[x] + cost[i] - h[y]);
                }
            }
        }
        return;
    }

    void tarjan(int x){
        // std::cerr << x << " " << init_n << std::endl;
        dfn[x]=low[x]=++timestamp;//赋予一个新的时间戳
        stk[++top]=x;//入栈
        in[x]=1;
        int y;
        for(int i=head[x];i;i=nex[i])if(e[i] && h[x] + cost[i] == h[y=ver[i]]){
            if(!dfn[y]){
                tarjan(y);
                low[x]=std::min(low[x],low[y]);
            }else if(in[y]) low[x]=std::min(low[x],dfn[y]);
        }
        if(dfn[x]==low[x]){//是强连通分量的第一个节点
            ++scc_num;
            int z;
            do{
                z=stk[top--];
                in[z]=0;
                id[z]=scc_num;
            }while(z!=x);
        }
        return;
    }

    
    void tarjan_no_cost(int x){
        dfn[x]=low[x]=++timestamp;//赋予一个新的时间戳
        stk[++top]=x;//入栈
        in[x]=1;
        int y;
        for(int i=head[x];i;i=nex[i])if(e[i]){
            y = ver[i];
            if(!dfn[y]){
                tarjan_no_cost(y);
                low[x]=std::min(low[x],low[y]);
            }else if(in[y]) low[x]=std::min(low[x],dfn[y]);
        }
        if(dfn[x]==low[x]){//是强连通分量的第一个节点
            ++scc_num;
            int z;
            do{
                z=stk[top--];
                in[z]=0;
                id[z]=scc_num;
            }while(z!=x);
        }
        return;
    }

    ui dfs_cnt;

    int dfs(ui x,ui F){
        ++dfs_cnt;
        // std::cerr << x << "x:F" << " " << F << std::endl;
        if(x == T || F == 0) return F;
        vis[x] = 1;
        int f,flow = 0;
        for(ui &i = cur[x];i;i = nex[i]){
            ui y = ver[i];
            // std::cerr << dep[y] << " " << dep[x] + cost[i] << " " << x << " " << y << " " << i << std::endl;
            if(!vis[y] && h[y] == h[x] + cost[i] && e[i]){
                f = dfs(y,std::min(F,e[i]));
                e[i] -= f;
                e[i^1] += f;
                flow += f;
                F -= f;
                if(F == 0) break;
            }
        }
        vis[x] = 0;
        return flow;
    }
    public:
    DinicSolver(ui _n,ui _look_ahead_prune){
        init_n = _n;
        head = new ui[6*init_n+5];
        cur = new ui[6*init_n+5];
        ver = new ui[6*init_n+5];
        nex = new ui[6*init_n+5];
        e = new ui[6*init_n+5];
        h = new int[6*init_n+5];
        cost = new int[6*init_n+5];
        deg = new ui[6*init_n+5];
        dep = new ui[6*init_n+5];
        dep2 = new ui[6*init_n+5];
        dep3 = new ui[6*init_n+5];
        vis = new ui[6*init_n+5];
        edge_id = new ui[6*init_n+5];
        max_deg = new ui[6*init_n+5];
        dfn = new ui[6*init_n+5];
        low = new ui[6*init_n+5];
        in = new ui[6*init_n+5];
        id = new ui[6*init_n+5];
        stk = new ui[6*init_n+5];
        tot = timestamp = scc_num = top = dfs_cnt = ahead_time = 0;
        look_ahead_prune = _look_ahead_prune;
        INIT_LIMIT = 6*init_n+5;
        Q = new ListLinearHeap(6*init_n+5,6*init_n+5);
    }
    ~DinicSolver(){
        if(head != nullptr){
            delete[] head;
            head = nullptr;
        }
        if(cur != nullptr){
            delete[] cur;
            cur = nullptr;
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
        if(h != nullptr){
            delete[] h;
            h = nullptr;
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
        if(dep2 != nullptr){
            delete[] dep2;
            dep2 = nullptr;
        }
        if(dep3 != nullptr){
            delete[] dep3;
            dep3 = nullptr;
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
        if(dfn != nullptr){
            delete[] dfn;
            dfn = nullptr;
        }
        if(low != nullptr){
            delete[] low;
            low = nullptr;
        }
        if(in != nullptr){
            delete[] in;
            in = nullptr;
        }
        if(id != nullptr){
            delete[] id;
            id = nullptr;
        }
        if(stk != nullptr){
            delete[] stk;
            stk = nullptr;
        }
    }
    void init(ui _n1, ui _n2, ui _m, ui *from, ui *to, ui *val, ui *id, ui K, ui lim, ui check_id){
        check_ID = check_id;
        tot = 1;
        // std::cerr << "n:" << n << std::endl;
        // std::cerr << _n1 << " n1+n2 " << _n2 << std::endl;
        n = _n1+_n2+3;m = _m;
        n1 = _n1;n2 = _n2;
        // std::cerr << "n:" << n << " " << _n1 << " " << _n2 << " " << m << std::endl;
        S = n-1,D = n-2,T = n-3;
        memset(head,0,sizeof(ui)*n);
        memset(deg,0,sizeof(ui)*n);
        memset(max_deg,0,sizeof(ui)*n);
        memset(h,0,sizeof(int)*n);
        memset(h + n1,0x3f,sizeof(int)*(n2+1));
        for(ui i = 0;i < m;++ i){
            ui u = from[id[i]],v = to[id[i]];
            add(u,v,1,val[id[i]]);
            edge_id[tot] = edge_id[tot^1] = id[i];
            ++deg[u];
            ++deg[v];
            h[v] = std::min(h[v],(int)val[id[i]]);
            h[T] = std::min(h[T],(int)val[id[i]]);
        }
        edge_tot = tot;
        for(ui i = 0;i < _n1;i ++){
            for(ui j = 0,s = 0,bef_s = 0;j < deg[i] && bef_s <= K+1;bef_s=s,s+=j,j ++){
                add(D,i,1,j);
            }
        }
        for(ui i = _n1;i < _n1 + _n2;i ++){
            for(ui j = 0,s = 0,bef_s = 0;j < deg[i] && bef_s <= K+1;bef_s=s,s+=j,j ++){
                add(i,T,1,j);
            }
        }
        add(S,D,lim+1,0);
        last_edge = tot-1;
        flow_lim = lim;
        return;
    }
    int run(ui *rid, ui *net_loss, ui *color_buf_1, ui *color_buf_2, ui &reduce_val, ui K){
        int fans = 0;
        //std::cerr << "--T--T--" << std::endl;
        do{
            reduce_val = h[T];
            if(flow_lim == 0) break;
            vis[T] = 1;
            while(vis[T]){
                memset(vis,0,sizeof(ui)*n);
                memcpy(cur,head,sizeof(ui)*n);
                dfs_cnt = 0;
                ui f = dfs(S,flow_lim);
                flow_lim -= f;
                fans += f*h[T];
                if(f) reduce_val = h[T];
                if(fans > K) break;
            }
            if(fans > K) break;
        }while(dijk(dep));
        
        if(fans + reduce_val > K) return fans;

        for(ui i = 0;i < n1;++i){
            for(ui j = head[i];j;j = nex[j]){
                if(ver[j] != D && e[j] == 0){
                    rid[edge_id[j]] = 1;
                }
            }
        }

        if(look_ahead_prune){
            inv_dijk(dep);
            top = timestamp = scc_num = 0;
            memset(dfn,0,sizeof(ui)*n);
            for(ui i = 0;i < n;++i){
                if(dfn[i] == 0)
                tarjan(i);
            }
            
            for(ui i = 0;i < n1;++i){
                net_loss[i] = 0;
                for(ui j = head[i];j;j = nex[j]){
                    if(ver[j] != D && e[j] == 1){
                        net_loss[edge_id[j]] = dep[ver[j]] + h[T] - h[ver[j]] + h[i] - h[D] + cost[j];
                        if(id[ver[j]] == id[i] && h[i] + cost[j] == h[ver[j]])
                        rid[edge_id[j]] = 2;
                        else
                        rid[edge_id[j]] = 0;
                    }
                    else if(ver[j] != D && e[j] == 0){
                        if(h[i] + cost[j] != h[ver[j]])
                        rid[edge_id[j]] = 3;
                    }
                    else if(ver[j] == D && e[j] == 0 && id[i] == id[ver[j]] && h[i] + cost[j] == h[ver[j]]){
                        color_buf_1[i] = 1;
                    }
                }
            }
            
            for(ui i = n1;i < n1+n2; ++i)
                for(ui j = head[i];j;j = nex[j])
                    if(ver[j] == T && e[j] == 1 && id[i] == id[ver[j]] && h[i] + cost[j] == h[ver[j]]){
                        color_buf_2[i-n1] = 1;
                        break;
                    }
        }
        return fans;
    }
};

#endif
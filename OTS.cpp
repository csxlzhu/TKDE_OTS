#include <bits/stdc++.h>

using namespace std;
typedef pair <int, int> P;
const int N = 500050;
const int Maxk = 30;
vector <double> f[N][Maxk];
vector <vector <int> > S[N][Maxk];
//vector <int> ans[N][Maxk];
vector <int> anc[N];
vector <int> I;
vector <P> vec[N];
map <int, int> mp_id;
int id[N], feq[N] = {}, fa[N], dep[N], K;
double g[N][Maxk] = {};
int n, num, cnt = 0, root = 0;

struct Graph{
	vector <int> vec[N];
	int cnt = 0, dfn[N], dep[N], p[N][30];
	inline int lca(int x, int y){
		if (dep[x] < dep[y]) swap(x, y);
		for (int i = log2(dep[x]); i >= 0; i--)
            if (dep[x] - (1 << i) >= dep[y])
                x = p[x][i];
		if (x == y) return y;
		for (int i = log2(dep[x]); i >= 0; i--)
            if (p[x][i] ^ p[y][i])
                x = p[x][i], y = p[y][i];
		return p[x][0];
	}
	inline void ST(){
		int k = log2(n);
		for (int j = 1; j <= k; j++)
            for (int i = 1; i <= n; i++)
                if (p[i][j - 1] != -1) p[i][j] = p[p[i][j - 1]][j - 1];
	}
	inline void dfs(int u, int depth){
		dfn[u] = ++cnt;
		for (int i = 0; i < vec[u].size(); i++)
            //if(vec[u][i] != fa[u])
                p[vec[u][i]][0] = u,
                dfs(vec[u][i], dep[vec[u][i]] = depth + 1);
	}
	inline void init(int root){
		cnt = 0;
		p[root][0] = -1, dep[root] = 1;
		dfs(root, dep[root]);
        ST();
	}
	inline void Clear(int n){
        memset(p,-1,sizeof(p));
        for (int i = 0; i <= n; i++) vec[i].clear();
	}
}G;
bool cmp(int x, int y){
    return G.dfn[x] < G.dfn[y];
}
void input(){
    //cout << cnt++ << endl;
    //FILE *fin, *fout;
    //fin = fopen("tree-MED.txt", "r");
    //fin = fopen("test_tree.txt", "r");
    //freopen("data.in", "r", stdin);
    //freopen("dp.out", "w", stdout);
    scanf("%d", &n);
    memset(fa, -1, sizeof fa);
    mp_id.clear();
    G.Clear(n);
    cnt = 0;
    for (int i = 1; i < n; i++)
	{
	    int x, y;
		scanf("%d %d", &x, &y);
		//if (y == 49) printf("%d %d\n", x, y);
		if (mp_id.find(x) == mp_id.end()) mp_id[x] = cnt++, id[mp_id[x]] = x;
		if (mp_id.find(y) == mp_id.end()) mp_id[y] = cnt++, id[mp_id[y]] = y;
		G.vec[mp_id[x]].push_back(mp_id[y]);
		fa[mp_id[y]] = mp_id[x];
		//cout << mp_id[x] << " " << mp_id[y] << endl;
	}
	//cout << n << " " << cnt << endl;
	for (int i = 0; i < n; i++)
        if (fa[i] == -1){
                root = i; break;
        }
    G.init(root);
    //cout << root << endl;
	//fclose(fin);
	//fin = fopen("test_feq.txt", "r");
	//fin = fopen("freq2.txt", "r");
	//fin = fopen("freq3.txt", "r");
	//fin = fopen("feq-patient.txt", "r");
	//fin = fopen("feq-nurse.txt", "r");
	memset(feq, 0, sizeof feq);
	int T;
	scanf("%d", &T);
	//int sum = 0;
	for (int i = 0; i < T; i++)
	{
	    int x, y;
		scanf("%d %d", &x, &y);
		//if (mp_id[x] == 0 && y != 0) cout << "fuck" << endl;
		if (y != 0 || mp_id[x] == root) {
            feq[mp_id[x]] = y;
            I.push_back(mp_id[x]);
		}
		//sum += y;
	}
	sort(I.begin(), I.end(), cmp);
}
void add(int u, int v){
    if (u == v) return ;
    vec[u].push_back(P(v, G.dep[v] - G.dep[u]));
    fa[v] = u;
    //cout << u << " " << v << " " << G.dep[v] - G.dep[u] << endl;
}
int stk[N] = {};
void get_vtree(){
    for (int i = 0; i <= n; i++)
        vec[i].clear(), anc[i].clear();
    sort(I.begin(), I.end(), cmp);
    memset(stk, 0, sizeof stk);
    int cnt = 0, top = 0;
    stk[0] = root;
    stk[++top] = root;
    for (int i = 0; i < I.size(); i++){
        int now = I[i], lca = G.lca(now, stk[top]);
        //cout << now << " " << stk[top] << " " << lca << endl;
        while(G.dep[lca] < G.dep[stk[top - 1]]){
            add(stk[top - 1], stk[top]);
            top--;
        }
        add(lca, stk[top]);
        stk[top] = lca;
        if (stk[top] != now) stk[++top] = now;
    }
    while(--top) add(stk[top], stk[top + 1]);
}
void get_anc(int x, int len){
    dep[x] = len;
    //cout << x << " " << len << endl;
    anc[x].push_back(-1);
    if (fa[x] != -1){
        for (int i = 1; i < anc[fa[x]].size(); i++)
            anc[x].push_back(anc[fa[x]][i]);
    }
    anc[x].push_back(x);
    for (int i = 0; i < vec[x].size(); i++){
        int y = vec[x][i].first;
        get_anc(y, len + vec[x][i].second);
    }
}

//void get_val(){
//    for (int i = 0; i < n; i++)
//        if (feq[i])
//        for (int j = 1; j < anc[i].size(); j++){
//            //cout << i << " " << anc[i][j] << endl;
//            val[anc[i][j]][i] = 1.0 / (dep[i] - dep[anc[i][j]] + 1) * feq[i];
//        }
//
//}

double val(int x, int y){
    //printf("%d %d\n", dep[y], dep[x]);
    return 1.0 / (dep[y] - dep[x] + 1) * feq[y];
}

void pre(){
    input();
    //cout << "ok" << endl;
    get_vtree();
    //cout << "ok" << endl;
    get_anc(root, 0);
    //cout << "ok" << endl;
    //get_val();
    //cout << "ok" << endl;
}
void dfs(int x, int K){
    for (int i = 0; i < vec[x].size(); i++){
        int y = vec[x][i].first;
        dfs(y, K);
    }
    for (int j = 0; j < anc[x].size(); j++)
        for (int k = 0; k <= K; k++)
            f[x][k].push_back(0), S[x][k].push_back(vector<int>());
    for (int j = 0; j < anc[x].size() - 1; j++)
    for (int i = 0; i < vec[x].size(); i++){
        int y = vec[x][i].first, z = anc[x][j];
        //if (y == fa[x]) continue;
        double tmp[Maxk] = {};
        vector <int> tmpS[Maxk];
        for (int k1 = 0; k1 <= K; k1++) tmp[k1] = f[x][k1][j];
        for (int k1 = 0; k1 <= K; k1++){
            //double tmp = 0;
            //if (z != -1) f[x][k1 - k2][j] = val[z][x];
            for (int k2 = 0; k2 <= k1; k2++){
                if (f[x][k1 - k2][j] + f[y][k2][j] >= tmp[k1])
                    tmp[k1] = max(tmp[k1], f[x][k1 - k2][j] + f[y][k2][j]),
                    tmpS[k1].clear(),
                    tmpS[k1].insert(tmpS[k1].end(), S[x][k1 - k2][j].begin(), S[x][k1 - k2][j].end()),
                    tmpS[k1].insert(tmpS[k1].end(), S[y][k2][j].begin(), S[y][k2][j].end());
            }
            //if (z != -1) f[x][k1][j] += val[z][x];
            //cout << f[x][k1][j] << endl;
        }
        for (int k1 = 0; k1 <= K; k1++) f[x][k1][j] = tmp[k1], S[x][k1][j] = tmpS[k1], tmpS[k1].clear();
        //cout << x << " " << k1 << " " << z << " " << f[x][k1][j] << endl;
    }
    int j = anc[x].size() - 1;
    for (int i = 0; i < vec[x].size(); i++){
        int y = vec[x][i].first, z = anc[x][j];
        //if (y == fa[x]) continue;
        double tmp[Maxk] = {};
        vector <int> tmpS[Maxk];
        for (int k1 = 1; k1 <= K; k1++) tmp[k1] = f[x][k1][j];
        for (int k1 = 1; k1 <= K; k1++){
            //double tmp = 0;
            //if (z != -1) f[x][k1 - k2][j] = val[z][x];
            for (int k2 = 0; k2 < k1; k2++){
                if (f[x][k1 - k2][j] + f[y][k2][j] >= tmp[k1])
                    tmp[k1] = max(tmp[k1], f[x][k1 - k2][j] + f[y][k2][j]),
                    tmpS[k1].clear(),
                    tmpS[k1].insert(tmpS[k1].end(), S[x][k1 - k2][j].begin(), S[x][k1 - k2][j].end()),
                    tmpS[k1].insert(tmpS[k1].end(), S[y][k2][j].begin(), S[y][k2][j].end());
                    //cout << k1 << " " << tmpS[k1].size() << endl;
            }
            //if (z != -1) f[x][k1][j] += val[z][x];
            //cout << f[x][k1][j] << endl;
        }
        for (int k1 = 0; k1 <= K; k1++) f[x][k1][j] = tmp[k1], S[x][k1][j] = tmpS[k1], tmpS[k1].clear();
        //cout << x << " " << k1 << " " << z << " " << f[x][k1][j] << endl;
    }
    for (int j = 0; j < anc[x].size(); j++)
        for (int k1 = 0; k1 <= K; k1++){
            int z = anc[x][j];
            if (z != -1 && !(z == x && k1 == 0)) f[x][k1][j] += val(z, x);
    }
    j = anc[x].size() - 1;
    for (int k = 1; k <= K; k++) g[x][k] = f[x][k][j], S[x][k][j].push_back(id[x]);
    for (int j = 0; j < anc[x].size(); j++)
    for (int k = 0; k <= K; k++)
        if (f[x][k][j] <= g[x][k])
            S[x][k][j] = S[x][k][anc[x].size() - 1],
            f[x][k][j] = max(f[x][k][j], g[x][k]);
}
double solve(int K){
    for (int i = 0; i <= n; i++)
        for (int k = 0; k <= K; k++)
            f[i][k].clear(), S[i][k].clear();
    dfs(root, K);
    double ans = 0;
    vector <int> Ans;
    //cout << anc[root].size() << endl;
    for (int k = 0; k <= K; k++){
        if (ans < f[root][k][1]){
            Ans = S[root][k][1], ans = max(ans, f[root][k][1]);
        }
        if (ans < f[root][k][0]){
            Ans = S[root][k][0], ans = max(ans, f[root][k][0]);
        }
    }
    sort(Ans.begin(), Ans.end());
    for (int i = 0; i < Ans.size(); i++)
        cout << Ans[i] << endl;
    return ans;
}
int main(int argc, char *argv[])
{
    K = atoi(argv[3]);
    char input_name[1000], output_name[1000];
    sprintf(input_name, "%s", argv[1]);
    sprintf(output_name, "%s", argv[2]);
    freopen(input_name, "r", stdin);
    freopen(output_name, "w", stdout);
    clock_t t1;
	pre();
    //double t = 0;
    //t1 = clock();
	printf("%.3f\n", solve(K));
	//t = (clock() - t1) * 1.0 / CLOCKS_PER_SEC;
    //printf("%.3f\n", log(t) / log(10));
    return 0;
}

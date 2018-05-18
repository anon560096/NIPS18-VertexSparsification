// Last update May 16, 2018

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <vector>
using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int, int> ii;
typedef vector<ii> vii;
typedef long long ll;

struct Edge {
  int u, v;
  ll cap, flow;
  Edge() {}
  Edge(int u, int v, ll cap): u(u), v(v), cap(cap), flow(0) {}
};

struct Dinic {
  int N;
  vector<Edge> E;
  vector<vector<int>> g;
  vector<int> d, pt;
  
  Dinic(int N): N(N), E(0), g(N), d(N), pt(N) {}

  void AddEdge(int u, int v, ll cap) {
    if (u != v) {
      E.emplace_back(u, v, cap);
      g[u].emplace_back(E.size() - 1);
      E.emplace_back(v, u, 0);
      g[v].emplace_back(E.size() - 1);
    }
  }

  bool BFS(int S, int T) {
    queue<int> q({S});
    fill(d.begin(), d.end(), N + 1);
    d[S] = 0;
    while(!q.empty()) {
      int u = q.front(); q.pop();
      if (u == T) break;
      for (int k: g[u]) {
        Edge &e = E[k];
        if (e.flow < e.cap && d[e.v] > d[e.u] + 1) {
          d[e.v] = d[e.u] + 1;
          q.emplace(e.v);
        }
      }
    }
    return d[T] != N + 1;
  }

  ll DFS(int u, int T, ll flow = -1) {
    if (u == T || flow == 0) return flow;
    for (int &i = pt[u]; i < g[u].size(); ++i) {
      Edge &e = E[g[u][i]];
      Edge &oe = E[g[u][i]^1];
      if (d[e.v] == d[e.u] + 1) {
        ll amt = e.cap - e.flow;
        if (flow != -1 && amt > flow) amt = flow;
        if (ll pushed = DFS(e.v, T, amt)) {
          e.flow += pushed;
          oe.flow -= pushed;
          return pushed;
        }
      }
    }
    return 0;
  }

  ll MaxFlow(int S, int T) {
    ll total = 0;
    while (BFS(S, T)) {
      fill(pt.begin(), pt.end(), 0);
      while (ll flow = DFS(S, T))
        total += flow;
    }
    return total;
  }

  vvi Residual() {
    vvi adj(N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < g[i].size(); j++) {
        Edge e = E[g[i][j]];
        if (e.cap == 0) continue;
        if (e.flow < e.cap) adj[e.u].push_back(e.v);
      }
    }
    return adj;
  }
};

class UnionFind {
private:
  vi p, rank, setSize;
  int numSets;
public:
  UnionFind(int N) {
    setSize.assign(N, 1); numSets = N; rank.assign(N, 0);
    p.assign(N, 0); for (int i = 0; i < N; i++) p[i] = i; }
  int findSet(int i) { return (p[i] == i) ? i : (p[i] = findSet(p[i])); }
  bool isSameSet(int i, int j) { return findSet(i) == findSet(j); }
  void unionSet(int i, int j) { 
    if (!isSameSet(i, j)) { numSets--; 
    int x = findSet(i), y = findSet(j);
    if (rank[x] > rank[y]) { p[y] = x; setSize[x] += setSize[y]; }
    else                   { p[x] = y; setSize[y] += setSize[x];
                             if (rank[x] == rank[y]) rank[y]++; } } }
  int numDisjointSets() { return numSets; }
  int sizeOfSet(int i) { return setSize[findSet(i)]; }
};

int n, m, num_terminals;
UnionFind cluster(0);
vector<map<int, int>> adj, original_adj;
vi true_label;
set<int> terminals;
vector<int> terminal_list;
vector<map<int, int>> current_graph;
set<int> independent_set;
int approx_power = 1;

void PrintAdj() {
  for (int i = 0; i < n; i++) {
    cout << i << ":";
    for (auto p : adj[i]) {
      cout << " (" << p.first << ", " << p.second << ")";
    }
    cout << endl;
  } 
}

int RandomNeighbor(int node) {
  int degree = 0;
  for (auto key_pair : adj[node]) {
    degree += key_pair.second;
  }
  double r = ((double)random()/(double)RAND_MAX); // bad idea, but fine for now
  int running_sum = 0;
  int ans = -1;
  for (auto key_pair : adj[node]) {
    running_sum += key_pair.second;
    if (r*degree <= running_sum) {
      ans = key_pair.first;
      break;
    }
  }
  assert(ans != -1);
  return ans;
}

void Contract(int u, int v) {
  // Contract u and v (keep v's index)
  //cout << u << " " << v << endl;
  cluster.unionSet(u, v);
  // Build new neighborhood for v
  map<int, int> new_neighorhood;
  for (auto key_pair : adj[u]) {
    if (key_pair.first != v) {
      new_neighorhood[key_pair.first] += key_pair.second;
    }
  }
  for (auto key_pair : adj[v]) {
    if (key_pair.first != u) {
      new_neighorhood[key_pair.first] += key_pair.second;
    }
  }
  // Delete edges pointing into u and v
  for (auto key_pair : adj[u]) {
    int neighbor = key_pair.first;
    adj[neighbor].erase(u);
  }
  for (auto key_pair : adj[v]) {
    int neighbor = key_pair.first;
    adj[neighbor].erase(v);
  }
  // Delete edges pointing out of u and v
  adj[u].clear();
  adj[v].clear();
  // Rebuild 
  for (auto key_pair : new_neighorhood) {
    int neighbor = key_pair.first;
    adj[v][neighbor] = key_pair.second;
    adj[neighbor][v] = key_pair.second;
  }
}

void ContractGraph(vi &permutation) {
  //PrintAdj();
  //cout << endl;
  current_graph = adj;
  for (int i = 0; i < n; i++) {
    int node = permutation[i];
    if (terminals.count(node)) continue;
    // Check to see if still an independent set.
    bool reset_current_graph = false;
    for (auto key_pair : current_graph[node]) {
      int neighbor = key_pair.first;
      if (independent_set.count(neighbor)) {
        independent_set.clear();
        reset_current_graph = true;
        approx_power++;
        break;
      }
    }
    independent_set.insert(node);
    int neighbor = RandomNeighbor(node);
    Contract(node, neighbor);
    if (reset_current_graph) current_graph = adj;
    //PrintAdj();
    //cout << endl;
  }
}

int MaxFlow(int s, int t) {
  Dinic G = Dinic(n);
  for (int i = 0; i < n; i++) {
    for (auto key_pair : adj[i]) {
      int neighbor = key_pair.first;
      int capacity = key_pair.second;
      G.AddEdge(i, neighbor, capacity);
    }
  }
  return G.MaxFlow(s, t);
}

int best_cut_s = -1;
int best_cut_t = -1;
int GlobalMinCut() {
  vi candidates;
  for (int i = 0; i < n; i++) {
    if (adj[i].size()) candidates.push_back(i);
  }
  int min_cut = 1 << 30;
  for (int i = 1; i < candidates.size(); i++) {
    int s = candidates[0];
    int t = candidates[i];
    //cout << s << " " << t << endl;
    int flow_val = MaxFlow(s, t);
    if (flow_val < min_cut) {
      min_cut = flow_val;
      best_cut_t = t;
    }
  }
  best_cut_s = candidates[0];
  MaxFlow(candidates[0], best_cut_t);
  return min_cut;
}

vii GetRandomTerminalPairs(int k) {
  vii terminal_pairs;
  for (int i = 0; i < k; i++) {
    vi tmp_terminals = terminal_list;
    random_shuffle(tmp_terminals.begin(), tmp_terminals.end());
    terminal_pairs.push_back(ii(tmp_terminals[0], tmp_terminals[1]));
  }
  return terminal_pairs;
}

void RecoverPartition() {
  // Uses global s and t from earlier.
  // NOTE: No longer consider this approach for stochastic block model.
  cout << best_cut_s << " " << best_cut_t << endl;
  Dinic G = Dinic(n);
  for (int i = 0; i < n; i++) {
    for (auto key_pair : adj[i]) {
      int neighbor = key_pair.first;
      int capacity = key_pair.second;
      G.AddEdge(i, neighbor, capacity);
    }
  }
  int flow_val = G.MaxFlow(best_cut_s, best_cut_t);

  vvi residual = G.Residual();
  vi ans;
  set<int> seen;
  queue<int> q;
  q.push(best_cut_s);
  seen.insert(best_cut_s);
  while (!q.empty()) {
    int cur_node = q.front();
    q.pop();
    ans.push_back(cur_node);
    cout << " - " << cur_node << endl;
    for (auto neighbor : residual[cur_node]) {
      if (!seen.count(neighbor)) {
        seen.insert(neighbor);
        q.push(neighbor);
      }
    }
  }
  cout << ans.size() << endl;
  // Find edges in cut
  set<ii> edges_in_cut;
}

// Comput isolating cut
pair<int, set<ii>> IsolatingCut(int s, set<int> T) {
  int super_sink = n;
  Dinic G = Dinic(n + 1);
  // Construct modified graph
  for (int i = 0; i < n; i++) {
    for (auto key_pair : adj[i]) {
      int neighbor = key_pair.first;
      int capacity = key_pair.second;
      int from = i;
      int to = neighbor;
      if (T.count(from)) from = super_sink;
      if (T.count(to)) to = super_sink;
      G.AddEdge(from, to, capacity);
    }
  }
  int flow_val = G.MaxFlow(s, super_sink);

  // Get edges in cut
  vvi residual = G.Residual();
  set<int> ans;
  set<int> seen;
  queue<int> q;
  q.push(s);
  seen.insert(s);
  while (!q.empty()) {
    int cur_node = q.front();
    q.pop();
    ans.insert(cur_node);
    for (auto neighbor : residual[cur_node]) {
      if (!seen.count(neighbor)) {
        seen.insert(neighbor);
        q.push(neighbor);
      }
    }
  }
  set<ii> edges_in_cut;
  for (int i = 0; i < n; i++) {
    for (auto key_pair : adj[i]) {
      int neighbor = key_pair.first;
      if (neighbor < i) continue;
      int capacity = key_pair.second;
      if (ans.count(i) != ans.count(neighbor)) {
        //cout << " - " << i << " " << neighbor << endl;
        edges_in_cut.insert(make_pair(i, neighbor));
      }
    }
  }
  return make_pair(flow_val, edges_in_cut);
}

double CutQualityMean(const vector<int> &exact, const vector<int> &approx) {
  int trials = exact.size();
  double total = 0;
  for (int i = 0; i < trials; i++) {
    total += (double)approx[i]/exact[i];
  }
  return total/trials;
}

double CutQualityStd(const vector<int> &exact, const vector<int> &approx) {
  int trials = exact.size();
  double total = 0;
  for (int i = 0; i < trials; i++) {
    total += pow((double)approx[i]/exact[i], 2);
  }
  double variance = total/trials - CutQualityMean(exact, approx);
  return sqrt(variance);
}

double CutQualityMax(const vector<int> &exact, const vector<int> &approx) {
  double ans = 0;
  for (int i = 0; i < exact.size(); i++) {
    ans = max(ans, (double)approx[i]/exact[i]);
  }
  return ans;
}

void PrintStats(const vector<int> &exact, const vector<int> &approx) {
  int trials = exact.size();
  cout << "NumberOfTerminals: " << terminal_list.size() << endl;
  cout << "NumberOfTrials: " << trials << endl;
  cout << "CutQuality_MEAN: " << CutQualityMean(exact, approx) << endl;
  cout << "CutQuality_STD: " << CutQualityStd(exact, approx) << endl;
  cout << "CutQuality_MAX: " << CutQualityMax(exact, approx) << endl;
}

ll k_way_cut(vi k_terminals) {
  int k = k_terminals.size();
  // Compute isolating cut for each of the k terminals
  vector< pair<int, set<ii>> > all_isolating_cuts(k);
  for (int i = 0; i < k; i++) {
    int s = k_terminals[i];
    set<int> T;
    for (int j = 0; j < k; j++) {
      if (i == j) continue;
      T.insert(k_terminals[j]);
    }
    all_isolating_cuts[i] = IsolatingCut(s, T);
    //cout << " - " << i << ": " << s << " -> " << all_isolating_cuts[i].first << endl;
  }
  sort(all_isolating_cuts.begin(), all_isolating_cuts.end());

  set<ii> multiway_cut;
  for (int i = 0; i < k - 1; i++) {
    for (const auto &edge : all_isolating_cuts[i].second) {
      multiway_cut.insert(edge);
    }
  }
  ll total_cut_weight = 0;
  for (const auto &edge : multiway_cut) {
    int capacity = adj[edge.first][edge.second];
    total_cut_weight += capacity;
  }
  return total_cut_weight;
}

double k_way_cut_experiment(int k) {
  int trials = 50;
//  cout << "running " << trials << " trials with k = " << k << endl;
  // Using isolating cut heuristic now.
  vi exact_min_cut(trials), approx_min_cut(trials);
  for (int trial = 0; trial < trials; trial++) {
//    if ((trial + 1) % 10 == 0) {
//      cout << " - trial: " << trial + 1 << endl;
//    }
    adj = original_adj;

    // Compute random k terminals
    vi terminal_list_copy = terminal_list;
    random_shuffle(terminal_list_copy.begin(), terminal_list_copy.end());
    vi k_terminals(k);
    for (int j = 0; j < k; j++) k_terminals[j] = terminal_list_copy[j];

    exact_min_cut[trial] = k_way_cut(k_terminals);

    // Now perform random contractions down to terminals.
    cluster = UnionFind(n);
    vi permutation(n);
    for (int i = 0; i < n; i++) permutation[i] = i;
    random_shuffle(permutation.begin(), permutation.end());
    ContractGraph(permutation);

    approx_min_cut[trial] = k_way_cut(k_terminals);

//    cout << "trial: " << trial;
//    cout << " exact value: " << exact_min_cut[trial];
//    cout << "  approx value: " << approx_min_cut[trial];
//    cout << " --> ratio: " << (double)approx_min_cut[trial]/exact_min_cut[trial] << endl;
  }
  //cout << " " << CutQualityMean(exact_min_cut, approx_min_cut);
  //PrintStats(exact_min_cut, approx_min_cut);
  return CutQualityMean(exact_min_cut, approx_min_cut);
}

vi ParseLine(string line) {
  vi ret;
  stringstream ss(line);
  int tmp;
  while (ss >> tmp) ret.push_back(tmp);
  return ret;
}

void ReadGraph() {
  map<int, int> indexer;
  int num_nodes = 0;
  int num_edges = 0;
  string line;
  getline(cin, line); // skip first line
  map<ii, int> edges;
  while (getline(cin, line)) {
    vi vals = ParseLine(line);
    if (vals.size() == 2) vals.push_back(1);
    int u = vals[0];
    int v = vals[1];
    int w = vals[2];
    edges[ii(u, v)] = w;
    num_edges++;
    if (!indexer.count(u)) {
      indexer[u] = num_nodes;
      num_nodes++;
    }
    if (!indexer.count(v)) {
      indexer[v] = num_nodes;
      num_nodes++;
    }
  }
  n = num_nodes;
  m = num_edges;
  cout << "READING INPUT" << endl;
  cout << "graph with " << n << " vertices and " << m << " edges" << endl;
  adj = vector<map<int, int>>(n);
  for (auto key_pair : edges) {
    int u = key_pair.first.first;
    int v = key_pair.first.second;
    int w = key_pair.second;
    adj[u][v] = w;
    adj[v][u] = w;
  }
}

void GenerateRandomTerminals(int num_terminals) {
  terminals.clear();
  terminal_list.clear();
  vi v(n);
  for (int i = 0; i < n; i++) v[i] = i;
  random_shuffle(v.begin(), v.end());
  for (int i = 0; i < num_terminals; i++) {
    terminals.insert(v[i]);
    terminal_list.push_back(v[i]);
  }
}

int main() {
  // Init
  srand(time(0));
  ReadGraph();
  original_adj = adj;

  // Init k-cut vals
  vi k_vals = {2, 4, 8};
  //vi terminal_vals = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  vi terminal_vals = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200};
  //vi terminal_vals = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
  
  for (auto num_terminals : terminal_vals) {
    GenerateRandomTerminals(num_terminals);
    vector<double> buffer;
    buffer.push_back(num_terminals);
    for (auto k : k_vals) {
      //cout << "k: " << k << endl;
      buffer.push_back(k_way_cut_experiment(k));
    }
    cout << num_terminals << " terminals: ";
    for (int i = 0; i < k_vals.size(); i++) {
      cout << "mean " << k_vals[i] << "-way cut quality:" << buffer[i + 1] << ",";
    }
    cout << endl;
  }
  cout << endl;

  int factor = 1;
  while (1) {
    factor *= 2;
    int num_terminals = n/factor;
    if (num_terminals < 8) break;
    GenerateRandomTerminals(num_terminals);
    vector<double> buffer;
    buffer.push_back(factor);
    for (auto k : k_vals) {
      //cout << "k: " << k << endl;
      buffer.push_back(k_way_cut_experiment(k));
    }
    cout << num_terminals << " terminals: ";
    for (int i = 0; i < k_vals.size(); i++) {
      cout << "mean " << k_vals[i] << "-way cut quality:" << buffer[i + 1] << ",";
    }
    cout << endl;
  }
  return 0;
}

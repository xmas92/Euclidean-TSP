//
//  main.cpp
//  Euclidean TSP
//
//  Created by Axel Boldt-Christmas on 21/10/15.
//  Copyright © 2015 Axelnet. All rights reserved.
//

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>
#include <set>
#include <cassert>
#include <random>

std::vector<std::vector<int32_t> > C;
std::vector<std::pair<float,float> > V;
std::vector<int16_t> T;
int16_t N;
int32_t AvrC;

struct edge_t {
    int16_t _n1, _n2;
    edge_t(int16_t n1, int16_t n2) :_n1((n1<n2)?n1:n2), _n2((n1<n2)?n2:n1) {}
    int16_t other(int16_t o) const {
        //assert(_n1==o || _n2 == o);
        return (_n1 == o)?_n2:_n1;
    }
};
inline bool operator==(const edge_t& lhs, const edge_t& rhs){
    return (lhs._n1 == rhs._n1 && lhs._n2 == rhs._n2);
}
inline bool operator!=(const edge_t& lhs, const edge_t& rhs){return !(lhs == rhs);}
inline bool operator<(const edge_t& lhs, const edge_t& rhs){
    if (lhs._n1 == rhs._n1)
        return lhs._n2 < rhs._n2;
    return lhs._n1 < rhs._n1;
}

std::ostream& operator<<(std::ostream& os, const edge_t& e) {
    os << "(" << e._n1 << "--" << e._n2 << ")";
    return os;
}


std::ostream& operator<<(std::ostream& os, const std::array<int16_t,3>& n) {
    os << "{" << n[0] << "," << n[1] << ","<< n[2] << ")";
    return os;
}

struct node_t {
    int16_t _in, _out, _order;
    node_t(int16_t in = -1, int16_t out = -1, int16_t order = -1)
    : _in(in) , _out(out), _order(order) {}
    void swap(int16_t o = -1) { std::swap(_in, _out); _order = o; }
};

std::ostream& operator<<(std::ostream& os, const node_t& n) {
    os << "(" << n._in << "->x->" << n._out << ")@" << n._order;
    return os;
}

struct cycle_t {
    std::vector<node_t> _nodes;
    int16_t _n;
    cycle_t(int16_t n = 0): _n(n), _nodes(std::vector<node_t>(n)) {}
    void fixOrder() {
        int i, o;
        i = o = 0;
        do {
            _nodes[i]._order = o++;
            i = _nodes[i]._out;
        } while(i != 0);
    }
    int64_t distance() {
        int64_t d = 0;
        for (int i = 0; i < _n; i++) {
            d += C[_nodes[i]._out][i];
        }
        return d;
    }
    int32_t distance(int16_t n) {
        return C[n][_nodes[n]._out];
    }
    void order(int16_t from, int16_t to, int16_t order) {
        while (from != to) {
            _nodes[from]._order = (order++ % _n);
            from = _nodes[from]._out;
        }
        _nodes[from]._order = (order++ % _n);
    }
    void reverse(int16_t from, int16_t to, int16_t order) {
        while (from != to) {
            _nodes[from].swap(order++ % _n);
            from = _nodes[from]._out;
        }
        _nodes[from].swap(order++ % _n);
    }
    int32_t optGain2(int16_t n1,int16_t n2) {
        return distance(n1)+distance(n2)-C[n1][n2]-C[_nodes[n1]._out][_nodes[n2]._out];
    }
    int32_t optGain3a(std::array<int16_t,3> n, std::array<int16_t,3> no) {
        return distance(n[0])+distance(n[1])+distance(n[2])
        -C[n[0]][no[1]]
        -C[n[2]][no[0]]
        -C[n[1]][no[2]];
    }
    int32_t optGain3b(std::array<int16_t,3> n, std::array<int16_t,3> no) {
        return distance(n[0])+distance(n[1])+distance(n[2])
        -C[n[0]][no[1]]
        -C[n[2]][n[1]]
        -C[no[0]][no[2]];
    }
    void swap2(int16_t n1,int16_t n2) {
        int16_t n1out = _nodes[n1]._out;
        int16_t n2out = _nodes[n2]._out;
        if ((_n-_nodes[n1out]._order+_nodes[n2]._order)%_n > _n/2) {
            std::swap(n1, n2);
            std::swap(n1out, n2out);
        }
        reverse(n2, n1out, _nodes[n1]._order+1);
        _nodes[n1]._out = n2; _nodes[n2]._in = n1;
        _nodes[n1out]._out = n2out; _nodes[n2out]._in = n1out;
    }
    void swap3a(std::array<int16_t,3> n, std::array<int16_t,3> no)  {
        _nodes[n[0]]._out = no[1]; _nodes[no[1]]._in = n[0]; // E1
        order(no[1], n[2], _nodes[n[0]]._order + 1); // Update Order
        _nodes[n[2]]._out = no[0]; _nodes[no[0]]._in = n[2]; // E2
        order(no[0], n[1], _nodes[n[2]]._order + 1); // Update Order
        _nodes[n[1]]._out = no[2]; _nodes[no[2]]._in = n[1]; // E3
        //fixOrder();
        //check();
    }
    void swap3b(std::array<int16_t,3> n, std::array<int16_t,3> no)  {
        _nodes[n[0]]._out = no[1]; _nodes[no[1]]._in = n[0]; // E1
        order(no[1], n[2], _nodes[n[0]]._order + 1); // Update Order
        reverse(n[1], no[0], _nodes[n[2]]._order+1); // Reverse Edge
        _nodes[n[2]]._out = n[1]; _nodes[n[1]]._in = n[2]; // E2
        _nodes[no[0]]._out = no[2]; _nodes[no[2]]._in = no[0]; // E3
        //fixOrder();
        //check();
    }
    void check() {
        int16_t i = 0, c = 0;
        do {
            c++;
            assert(_nodes[_nodes[i]._out]._order == (_nodes[i]._order + 1) % _n);
            assert(i == _nodes[_nodes[i]._out]._in);
            i = _nodes[i]._out;
            assert(c <= _n);
        } while (i != 0);
        assert(c == _n);
    }
    void print() {
        int16_t i = 0;
        do {
            std::cout << i << std::endl;
            i = _nodes[i]._out;
        } while (i != 0);
    }
    inline int64_t twoOpt(int16_t n1,int16_t n2) {
        if (n1 == _nodes[n2]._in || n1 == n2 || n1 == _nodes[n2]._in) return 0;
        int ret = optGain2(n1, n2);
        if (ret <= 0)
            return 0;
        swap2(n1,n2);
        //std::cout << "{" << n1 << "," << n2 << "}" << std::endl;
        return ret;
    }
    inline int64_t threeOpt(int16_t n1,int16_t n2,int16_t n3) {
        if (n1 == _nodes[n2]._in || n1 == n2 || n1 == _nodes[n2]._in) return 0;
        if (n1 == _nodes[n3]._in || n1 == n3 || n1 == _nodes[n3]._in) return 0;
        if (n2 == _nodes[n3]._in || n2 == n3 || n2 == _nodes[n3]._in) return 0;
        int ret;
        std::array<int16_t,3> n = {{n1,n2,n3}};
        std::sort(n.begin(),n.end(),
                  [&](int16_t &a, int16_t &b) { return _nodes[a]._order < _nodes[b]._order; });
        std::array<int16_t,3> no = {{_nodes[n[0]]._out,_nodes[n[1]]._out,_nodes[n[2]]._out}};
        ret = optGain3a(n, no);
        if (ret > 0) {
            swap3a(n, no);
            //std::cout << n << " : " << ret << std::endl;
            return ret;
        }
        for (int i = 0; i < 3; i++) {
            ret = optGain3b(n,no);
            if (ret > 0) {
                swap3b(n,no);
                //std::cout << n << " : " << ret << std::endl;
                return ret;
            }
            std::rotate(n.begin(), n.begin()+1, n.end());
            std::rotate(no.begin(), no.begin()+1, no.end());
        }
        return 0;
    }
};

struct graph_t {
    std::vector<std::set<edge_t>> _edges;
    int16_t _n;
    graph_t(int16_t n = 0) : _n(n), _edges(std::vector<std::set<edge_t>>(n)) {}
    void add_edge(edge_t e) {
        auto a = _edges[e._n1].insert(e);
        auto b = _edges[e._n2].insert(e);
        //assert(a.second == b.second);
    }
    void remove_edge(edge_t e) {
        _edges[e._n1].erase(e);
        _edges[e._n2].erase(e);
    }
    int16_t degree(int16_t n) {
        return _edges[n].size();
    }
    int16_t end(int16_t i) {
        int16_t n = (*_edges[i].begin()).other(i);
        while (degree(n) == 2) {
            std::array<int16_t, 2> neig;
            size_t j = 0;
            for (auto &e : _edges[n]) {
                neig[j++] = e.other(n);
            }
            if (neig[0] == i) {
                i = n;
                n = neig[1];
            } else {
                i = n;
                n = neig[0];
            }
        }
        return n;
    }
    cycle_t cycle() {
        cycle_t ret(_n);
        int16_t p = -1, i = 0;
        do {
            //assert(degree(i) == 2);
            std::array<int16_t, 2> neig;
            size_t j = 0;
            for (auto &e : _edges[i]) {
                neig[j++] = e.other(i);
            }
            if (neig[0] == p) {
                p = i;
                i = neig[1];
            } else {
                p = i;
                i = neig[0];
            }
            ret._nodes[p]._out = i;
            ret._nodes[i]._in = p;
        } while (i != 0);
        ret.fixOrder();
        return ret;
    }
};

cycle_t Cyc;

auto Deadline = std::chrono::system_clock::now()+std::chrono::seconds(2);

void ReverseTour(int16_t s, int16_t e) {
    int16_t c = s;
    int16_t n = T[c];
    int16_t nn;
    do {
        nn = T[n];
        T[n] = c;
        c = n;
        n = nn;
    } while (c != e);
    
}

void ReadInput() {
    std::cin >> N;
    C = std::vector<std::vector<int32_t>>(N,std::vector<int32_t>(N,0));
    T = std::vector<int16_t>(N,-1);
    for (int i = 0; i < N; i++) {
        float x,y;
        std::cin >> x >> y;
        V.push_back(std::make_pair(x, y));
    }
}

void CalculateC() {
    int64_t t = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            t += C[j][i] = C[i][j] = std::round(std::sqrt(std::pow(V[i].first-V[j].first, 2)
                                           + std::pow(V[i].second-V[j].second, 2)));
        }
    }
    AvrC = static_cast<int32_t>(t/(N*(N-1)/2));
}

void WriteOutput() {
    int16_t i = 0;
    do {
        std::cout << i << std::endl;
        i = T[i];
    } while (i != 0);
}

int64_t Distance() {
    int16_t c = 0;
    int64_t t = 0;
    do {
        t += C[c][T[c]];
        c = T[c];
    } while (c != 0);
    return t;
}

bool IsTour() {
    int16_t c = 0;
    int16_t t = 0;
    do {
        t++;
        c = T[c];
    } while (c != 0);
    return t==N;
}

int64_t TwoOpt(int16_t s) {
    int64_t change = 0;
    int16_t o = -1;
    for (int64_t i = 0; i < N; i++) {
        if (i == s || i == T[s]) continue;
        if (C[s][T[s]]+C[i][T[i]]-C[s][i]-C[T[s]][T[i]] > change) {
            change = C[s][T[s]]+C[i][T[i]]-C[s][i]-C[T[s]][T[i]];
            o = i;
        }
    }
    if (o != -1) {
        int16_t pto = T[o];
        ReverseTour(T[s], o);
        T[T[s]] = pto;
        T[s] = o;
        //assert(IsTour());
    }
    return change;
}

bool Order(int16_t s1, int16_t s2, int16_t s3) {
    while (s1 != s3) {
        s1 = T[s1];
        if (s1 == s2)
            return true;
    }
    return false;
}

bool ThreeOpt(int16_t s1, int16_t s2) {
    if (s1 == s2 || s1 == T[s2] || s2 == T[s1]) return false;
    for (int64_t s3 = T[T[s2]]; s3 != s1; s3 = T[s3]) {
        //assert(Order(s1, s2, s3));
        if (s3 == s1 || s3 == T[s1] || s3 == s2 || s3 == T[s2]) continue;
        if (C[s1][T[s1]]+C[s2][T[s2]]+C[s3][T[s3]] > C[s1][T[s2]]+C[s3][T[s1]]+C[s2][T[s3]]) {
            int16_t pt1 = T[s1], pt3 = T[s3];
            T[s1] = T[s2];
            T[s3] = pt1;
            T[s2] = pt3;
            //assert(IsTour());
            return true;
        }
        if (C[s1][T[s1]]+C[s2][T[s2]]+C[s3][T[s3]] > C[s1][s3]+C[T[s2]][T[s1]]+C[s2][T[s3]]) {
            int16_t pt1 = T[s1], pt3 = T[s3];
            T[s1] = s3;
            ReverseTour(T[s2], s3);
            T[T[s2]] = pt1;
            T[s2] = pt3;
            //assert(IsTour());
            return true;
        }
    }
    return false;
}

void LocalCycOpt() {
    int64_t change = 0;
    int16_t s1 = 0, s2 = 0;
    int64_t counter = 0;
    while ((counter++ % (int)N != 0) || (Deadline-std::chrono::system_clock::now()) > std::chrono::milliseconds(100)) {
        change += Cyc.twoOpt(s1, s2++);
        if (s2 == N) {
            s1++;
            s2 = 0;
            if (s1 == N) {
                if (change == 0) break;
                change = s1 = 0;
            }
        }
    }
    std::vector<std::vector<bool>> nn(N,std::vector<bool>(N,false));
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            nn[i][j] = nn[j][i] = (C[i][j] < AvrC*0.25);
        }
    }
    s1 = 0, s2 = 0;
    while ((counter < (int)N) || (Deadline-std::chrono::system_clock::now()) > std::chrono::milliseconds(100)) {
        if ((counter >= (int)N)) {
            counter = 0;
        }
        if (nn[s1][s2]) {
            for (int16_t s3 = 0; s3 < N; s3++) {
                if (nn[s1][s3] && nn[s2][s3]) {
                    change += Cyc.threeOpt(s1, s2, s3);
                    counter++;
                }
            }
        }
        if (++s2 == N) {
            s1++;
            s2 = 0;
            if (s1 == N) {
                if (change == 0) break;
                change = s1 = 0;
            }
        }
    }
}

void LocalOpt() {
    int16_t s1 = 0, s2 = 0;
    int64_t change = 0;
    while ((Deadline-std::chrono::system_clock::now()) > std::chrono::milliseconds(100)) {
        change += TwoOpt(s1++);
        if (s1 == N) {
            if (change == 0) break;
            change = s1 = 0;
        }
    }
    s1 = 0;
    s2 = T[T[s1]];
    std::cout << Distance() << std::endl;
    while ((Deadline-std::chrono::system_clock::now()) > std::chrono::milliseconds(100)) {
        if (ThreeOpt(s1, s2)) {
            //s1 = T[s1];
            s2 = T[T[s1]];
        }
        else
            s2 = T[s2];
        if (s2 == s1) {
            s1 = T[s1];
            s2 = T[T[s1]];
        }
    }
}

void GreedyTour() {
    graph_t Graph(N);
    std::list<edge_t> edges;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            edges.push_back(edge_t(i,j));
        }
    }
    edges.sort([&](const edge_t &a, const edge_t &b){ return C[a._n1][a._n2] < C[b._n1][b._n2];});
    int16_t n = 0;
    while (n != N-1) {
        edge_t e = edges.front();
        edges.pop_front();
        if ((Graph.degree(e._n1) < 2 && Graph.degree(e._n2) == 0) ||
            (Graph.degree(e._n2) < 2 && Graph.degree(e._n1) == 0)) {
            Graph.add_edge(e);
            n++;
        }
        else if (Graph.degree(e._n1) == 1 &&
                 Graph.degree(e._n2) == 1 &&
                 Graph.end(e._n1) != e._n2) {
            Graph.add_edge(e);
            n++;
        }
    }
    int i = 0;
    while (Graph.degree(i) == 2) {i++;}
    int j = i+1;
    while (Graph.degree(j) == 2) {j++;}
    Graph.add_edge(edge_t(i,j));
    Cyc = Graph.cycle();
}

void NNTour() {
    int16_t c = 0;
    Cyc = cycle_t(N);
    std::vector<bool> visited(N, false);
    for(;;) {
        visited[c] = true;
        int32_t minC = std::numeric_limits<int32_t>::max();
        int16_t minV = -1;
        for (int16_t i = 1; i < N; i++) {
            if (visited[i]) continue;
            if (C[c][i] < minC) {
                minC = C[c][i];
                minV = i;
            }
        }
        if (minV == -1) {
            Cyc._nodes[c]._out = 0;
            Cyc._nodes[0]._in = c;
            break;
        }
        Cyc._nodes[c]._out = minV;
        Cyc._nodes[minV]._in = c;
        c = minV;
    }
    Cyc.fixOrder();
}

void CWTour() {
    graph_t Graph(N);
    int16_t h = -1;
    int64_t minAC = std::numeric_limits<int64_t>::max();
    for (int i = 0; i < N; i++) {
        int64_t AC = std::accumulate(C[i].begin(), C[i].end(), 0);
        if (AC < minAC) {
            h = i;
            minAC = AC;
        }
    }
    auto savings = [h](const int16_t &i,const int16_t &j)->int32_t {
        return C[h][i]+C[h][j] - C[i][j];
    };
    int16_t Vs = N-1;
    std::list<edge_t> edges;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            if (j == h) continue;
            edges.push_back(edge_t(i,j));
        }
    }
    edges.sort([&](const edge_t &a, const edge_t &b){
        return savings(a._n1,a._n2) > savings(b._n1,b._n2); });
    while (Vs > 2) {
        edge_t e = edges.front();
        edges.pop_front();
        if ((Graph.degree(e._n1) < 2 && Graph.degree(e._n2) == 0) ||
            (Graph.degree(e._n2) < 2 && Graph.degree(e._n1) == 0)) {
            Graph.add_edge(e);
            if (Graph.degree(e._n1) == 2 || Graph.degree(e._n2) == 2) Vs--;
        }
        else if (Graph.degree(e._n1) == 1 &&
                 Graph.degree(e._n2) == 1 &&
                 Graph.end(e._n1) != e._n2) {
            Graph.add_edge(e);
            Vs -= 2;
        }
    }
    for (int i = 0; i < N; i++) {
        if (i == h) continue;
        if (Graph.degree(i) == 1) Graph.add_edge(edge_t(i,h));
    }
    Cyc = Graph.cycle();
}

void InitialTour() {
    int16_t c = 0;
    std::vector<bool> visited(N, false);
    for(;;) {
        visited[c] = true;
        int32_t minC = std::numeric_limits<int32_t>::max();
        int16_t minV = -1;
        for (int16_t i = 1; i < N; i++) {
            if (visited[i]) continue;
            if (C[c][i] < minC) {
                minC = C[c][i];
                minV = i;
            }
        }
        if (minV == -1) {
            T[c] = 0;
            break;
        }
        T[c] = minV;
        c = minV;
    }
}

void RandomInput() {
    N = 1000;
    C = std::vector<std::vector<int32_t>>(N,std::vector<int32_t>(N,0));
    T = std::vector<int16_t>(N,-1);
    for (int i = 0; i < N; i++) {
        float x,y;
        std::random_device rd;
        std::uniform_int_distribution<float> dist(.0f, 1e6);
        x = dist(rd);
        y = dist(rd);
        V.push_back(std::make_pair(x, y));
    }
}

int main(int argc, const char * argv[]) {
    RandomInput();
    //ReadInput();
    //Deadline = std::chrono::system_clock::now()+std::chrono::seconds(2);
    CalculateC();
    Deadline = std::chrono::system_clock::now()+std::chrono::seconds(2);
    NNTour();
    std::cout << "NNeigh: " << Cyc.distance() << std::endl;
    LocalCycOpt();
    std::cout << "locOpt: " << Cyc.distance() << std::endl;
    Deadline = std::chrono::system_clock::now()+std::chrono::seconds(2);
    GreedyTour();
    std::cout << "Greedy: " << Cyc.distance() << std::endl;
    LocalCycOpt();
    std::cout << "locOpt: " << Cyc.distance() << std::endl;
    Deadline = std::chrono::system_clock::now()+std::chrono::seconds(2);
    CWTour();
    std::cout << "CWrigh: " << Cyc.distance() << std::endl;
    LocalCycOpt();
    std::cout << "locOpt: " << Cyc.distance() << std::endl;
    //Cyc.print();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(Deadline-std::chrono::system_clock::now()).count() << "ms" << std::endl;
}

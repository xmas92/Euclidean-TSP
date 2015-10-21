//
//  main.cpp
//  Euclidean TSP
//
//  Created by Axel Boldt-Christmas on 21/10/15.
//  Copyright © 2015 Axelnet. All rights reserved.
//

#include <iostream>
#include <vector>
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
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = std::round(std::sqrt(std::pow(V[i].first-V[j].first, 2)
                                           + std::pow(V[i].second-V[j].second, 2)));
        }
    }
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

std::pair<int16_t, int16_t> Edge(int16_t x, int16_t y){
    if (x < y)
        return std::make_pair(x, y);
    return std::make_pair(y, x);
}

void LKOpt(int16_t s) {
    std::set<std::pair<int16_t, int16_t>> X, Y;
    std::vector<int16_t> optT(T);
    int32_t GOpt = 0;
    int32_t G = 0;
    int32_t g;
    int32_t gOpt;
    int16_t ln = s;
    int16_t f = T[ln];
    int16_t n, nf;
    int16_t lpn = -1;
    std::pair<int16_t, int16_t> x;
    int32_t yOptC;
    int32_t xC;
    do {
        n = -1;
        x = Edge(ln, f);
        xC = C[ln][f];
        if (X.count(x) != 0) break;
        for (int16_t pn = T[f]; n == -1 && pn != s; pn = T[pn]) {
            g = xC - C[f][pn];
            if (!(X.count(Edge(f,pn))==0 &&
                  G+g > 0 &&
                  Y.count(Edge(lpn, pn))==0&&
                  T[pn] != 0 &&
                  pn != T[f])) {
                lpn = pn;
                continue;
            }
            n = pn;
        }
        if (n != -1) {
            X.insert(x);
            Y.insert(Edge(f, n));
            yOptC = C[f][s];
            gOpt = G + (xC - yOptC);
            if (gOpt > GOpt) {
                GOpt = gOpt;
                optT = T;
                optT[s] = f;
            }
            G += xC - C[f][n];
            ReverseTour(f, lpn);
            nf = lpn;
            T[f] = n;
            ln = n;
            f = nf;
        }
    } while (n != -1);
    T = optT;
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
    N = 200;
    C = std::vector<std::vector<int32_t>>(N,std::vector<int32_t>(N,0));
    T = std::vector<int16_t>(N,-1);
    for (int i = 0; i < N; i++) {
        float x,y;
        std::random_device rd;
        std::uniform_int_distribution<float> dist(.0f, 10e6);
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
    InitialTour();
    std::cout << Distance() << std::endl;
    LocalOpt();
    std::cout << Distance() << std::endl;
    WriteOutput();
}

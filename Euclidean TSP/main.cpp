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

void LocalOpt() {
    int16_t s = 0;
    while ((Deadline-std::chrono::system_clock::now()) > std::chrono::milliseconds(100)) {
        LKOpt(s);
        ++s %= N;
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

int main(int argc, const char * argv[]) {
    ReadInput();
    CalculateC();
    InitialTour();
    LocalOpt();
    WriteOutput();
}

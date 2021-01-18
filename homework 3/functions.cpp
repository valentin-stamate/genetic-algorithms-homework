#include "./functions.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>     
#include <time.h>     
#include <chrono> 
#include <random>
#include <math.h>

using namespace std::chrono;
using namespace std;
mt19937 generator(0);

void readGraph(char* source, int G[500][500], int* l) {
    ifstream in;
    in.open(source);

    int n;
    in>>n;
    (*l) = n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in>>G[i][j];
        }
    }
}

void toNormalTime(long long int time, ofstream* out) {
    int hours = time / 3600000;
    time -= hours * 3600000;

    int min = time / 60000;
    time -= min * 60000;

    int sec = time / 1000;
    int msec = time % 1000;

    (*out)<<hours<<"h "<<min<<"min "<<sec<<"s "<<msec<<"ms\n";
}

long long getCurrentTime() {
    milliseconds ms = duration_cast< milliseconds >(
        system_clock::now().time_since_epoch()
    );
    return ms.count();
}

Member::Member() {}

void Member::init() {
    for (int i = 0; i < this->n; i++) {
        gene[i] = generateRandom(0, this->n - 1 - i);
    }
}

void Member::copy(const Member& m) {
    this->n = m.n;
    this->cost = m.cost;
    this->fitness = m.fitness;
    for (int i = 0; i < this->n; i++) {
        gene[i] = m.gene[i];
    }
}

void Member::randomMutation() {
    int r = generateRandom(0, this->n - 1);

    gene[r]++;
    int mod = this->n - 1 - r;
    if (mod == 0) {
        gene[r] = 0;;
    } else {
        gene[r] %= mod;
    }
}

void Member::mutate(int o, int pop) {
    float max_pb = MAX_PROB - M_PROB;
    float m_pb = (o / (1.0 * pop)) * max_pb + M_PROB;

    for (int i = 0; i < this->n; i++) {
        float r = generateFloatRandom(0, 1);

        if (r < m_pb) {
            gene[i]++;
            int mod = this->n - 1 - i;
            if (mod == 0) {
                gene[i] = 0;;
            } else {
                gene[i] %= mod;
            }
        }
    }
}

void Member::mutate() {

    for (int i = 0; i < this->n; i++) {
        float r = generateFloatRandom(0, 1);

        if (r < M_PROB) {
            gene[i]++;
            int mod = this->n - 1 - i;
            if (mod == 0) {
                gene[i] = 0;;
            } else {
                gene[i] %= mod;
            }
        }
    }
}

void Member::secondMutation() {
    int p = generateRandom(1, 20);// let's say

    for (int k = 0; k < p; k++) {
        
        int i = generateRandom(0, this->n - 1);
        
        gene[i]++;
        int mod = this->n - 1 - i;
        if (mod == 0) {
            gene[i] = 0;;
        } else {
            gene[i] %= mod;
        }
    }

}

void Member::mutateGreedy(int G[500][500]) {
    int newGene[500];

    for (int i = 0; i < this->n; i++) {
        newGene[i] = gene[i];
    }

    for (int i = 0; i < this->n; i++) {
        float r = generateFloatRandom(0, 1);

        if (r < M_G_PROB) {
            newGene[i]++;
            int mod = this->n - 1 - i;
            if (mod == 0) {
                newGene[i] = 0;;
            } else {
                newGene[i] %= mod;
            }

            int newCost = score(newGene, this->n, G);

            if (newCost < cost) {
                gene[i] = newGene[i];
                cost = newCost;
            } else {
                newGene[i] = gene[i];
            }

        }
    }
}

void Member::randomProb() {
    rProb = generateFloatRandom(0, 1);
}

void Member::strongMutation() {
    for (int i = 0; i < this->n; i++) {
        float r = generateFloatRandom(0, 1);

        if (r < 0.4) {

            int m = generateRandom(1, 4);

            gene[i]+=m;
            int mod = this->n - 1 - i;
            if (mod == 0) {
                gene[i] = 0;;
            } else {
                gene[i] %= mod;
            }
        }
    }
}

void smallCrossover(Member* p[2], Member* c[2]) {

    int s = generateRandom(1, p[0]->n - 2);

    for (int i = 0; i < s; i++) {
        c[0]->gene[i] = p[0]->gene[i];
        c[1]->gene[i] = p[1]->gene[i];
    }

    for (int i = s; i < p[0]->n; i++) {
        c[0]->gene[i] = p[1]->gene[i];
        c[1]->gene[i] = p[0]->gene[i];
    }

}

void sSmallCrossover(Member* p[2], Member* c[2]) {

    int s = generateRandom(1, p[0]->n - 2);

    int n = p[0]->n;

    for (int i = 0; i < n; i++) {
        int r = generateFloatRandom(0, 1);
        if (r < M_S_PROB) {
            c[0]->gene[i] = p[0]->gene[i];
        } else {
            c[0]->gene[i] = p[1]->gene[i];
        }

        r = generateFloatRandom(0, 1);
        if (r < M_S_PROB) {
            c[1]->gene[i] = p[0]->gene[i];
        } else {
            c[1]->gene[i] = p[1]->gene[i];
        }
    }
}


int generateRandom(int x, int y) {
 
    int r = rand() % (y - x + 1) + x;

    return r;
}

float generateFloatRandom(float x, float y) {
    uniform_real_distribution<float> distribution(x, y); // [x, y]
    long long ct = getCurrentTime() % 10000;

    float rn = distribution(generator) * ct;
    float dif = y - x;

    while (rn >= dif) {
        rn -= dif;
    }

    rn += x;

    return rn;
}

float fitness(int sc, int minScore) {
    float score = sc - minScore + 10;
    score = 1.0 / score + 2;

    float f = (pow(5.0, score) - 25) * 100;

    return f;
}

float saFitness(int sc) {
    float fit = (1.0 * sc) / 1000;
    
    return 1 / fit;
}


int score(int gene[500],int n, int G[500][500]) {
    int s = 0;

    int *cycle = decodeGene(gene, n);

    for (int i = 0; i < n - 1; i++) {
        int c = G[cycle[i]][cycle[i + 1]];
        s += c;
    }

    int x = cycle[n - 1];
    int y = cycle[0];

    s += G[x][y];

    return s;
}


int cycle[500];
int* decodeGene(int gene[500], int n) {
    vector<int> list;

    for (int i = 0; i < n; i++) {
        list.push_back(i);
    }

    for (int i = 0; i < n; i++) {
        cycle[i] = list[gene[i]];
        list.erase(list.begin() + gene[i]);
    }

    return cycle;
}

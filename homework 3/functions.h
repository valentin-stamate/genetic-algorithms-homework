#pragma once

#include <fstream>
using namespace std;

#define M_PROB 0.001
#define MAX_PROB 0.02

#define M_S_PROB 0.8

#define M_G_PROB 0.1

#define M_C_PROB 0.5
#define M_POP 100
#define M_GEN 1000

#define MAX_POP 400

#define SA_IT 100000
#define M_ET 10
#define S_SIZE 30

class Member {
    public:
    
    int n;
    
    int gene[500];
    float fitness;
    int cost;

    float rProb;

    Member();
    void copy(const Member&);
    void init();
    void mutate(int, int);
    void mutate();
    void secondMutation();
    void mutateGreedy(int[500][500]);
    void mutateGreedy();
    void strongMutation();
    void randomMutation();
    void randomProb();
};

struct byFitness {
    bool operator() (const Member &a, const Member &b) {
        return a.fitness > b.fitness;
    }
};

struct byProb {
    bool operator() (const Member &a, const Member &b) {
        return a.rProb < b.rProb;
    }
};

int generateRandom(int, int);

void toNormalTime(long long int, ofstream*);
long long getCurrentTime();
void readGraph(char*, int[500][500], int*);
float generateFloatRandom(float, float);



void smallCrossover(Member*[2], Member*[2]);
void sSmallCrossover(Member*[2], Member*[2]);

float fitness(int, int);
float saFitness(int);
int score(int[500], int, int[500][500]);
int* decodeGene(int[500], int);


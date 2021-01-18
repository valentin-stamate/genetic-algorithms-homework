#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdlib.h>     
#include <time.h>   
#include "functions.h"
using namespace std;

int G[500][500];
int dim;

void initPop();

void select();
void mutate();
void crossover();
void evaluate();
void sortPopulation();

int getMinScore();

float simulateEvolution();
float simulatedAnnealing();

void showPopulation();
void showMember(Member*);

Member *population;
int nPop = M_POP;

ifstream in;
ofstream out;

int main() {

    char input[100];
    char output[100];

    int inT;
    cout<<"Enter the input: ";
    cin>>inT;

    sprintf(input, "inputs/input_%d", inT);
    sprintf(output, "output_%d", inT);

    in.open(input);
    out.open(output);

    srand(time(NULL));
    
    population = new Member[MAX_POP];

    readGraph(input, G, &dim);

    // initPop();

    float ct_1 = getCurrentTime();
    // for (int i = 1; i <= S_SIZE; i++) {
    //     printf("Repeat: %d\n", i);
    //     int r = simulateEvolution();
    //     // out<<r<<" ";
    // }
    float ct_2 = getCurrentTime();
    // out<<"\n";
    // toNormalTime(ct_2 - ct_1, &out);

    ct_1 = getCurrentTime();
    for (int i = 1; i <= 1; i++) {
        printf("Repeat: %d\n", i);
        int r = simulatedAnnealing();
        out<<r<<" ";
    }
    ct_2 = getCurrentTime();
    out<<"\n";
    toNormalTime(ct_2 - ct_1, &out);
    return 0;
}

void initPop() {
    for (int i = 0; i < MAX_POP; i++) {
        population[i].n = dim;
        population[i].init();
    }
}

float simulatedAnnealing() {

    Member m;
    Member nbor;

    int finalCost = 2000000000;

    for (int i = 1; i <= SA_IT; i++) {
        m.n = dim;
        m.init();

        m.cost = score(m.gene, dim, G);
        // float cFit = fitness(m.cost);

        float t = 100;
        while (t > 10e-8) {
            nbor = m;
            nbor.randomMutation();
            nbor.cost = score(nbor.gene, dim, G);
            nbor.fitness = saFitness(nbor.cost);
            
            if (nbor.fitness < m.fitness) {
                m = nbor;
            } else if(generateFloatRandom(0, 1 - 0.00001) < exp( -abs(m.fitness - nbor.fitness) / t )) {
                m = nbor;
            }

            t *= 0.90;
        }

        finalCost = min(finalCost, m.cost);
        out<<finalCost<<", ";
        printf("Sa: it %d current fitness %f and cost: %d\n", i,m.fitness , finalCost);
    }

    return finalCost;

}

float simulateEvolution() {

    evaluate();
    
    for (int i = 1; i <= M_GEN; i++) {
        sortPopulation();
        printf("Generation %d with best found %f with value %d\n", i, population[0].fitness, population[0].cost);
        out<<population[0].cost<<", ";
        select();
        // showPopulation();

        if (i % 10 == 0) {
            for (int j = 0; j < nPop; j++) {
                population[j].mutate(j, nPop);
            }
        }

        mutate();
        crossover();
        evaluate();
    }

    sortPopulation();

    return population[0].cost;
}

void evaluate() {

    int minScore = getMinScore();

    for (int i = 0; i < nPop; i++) {
        population[i].cost = score(population[i].gene, dim, G);
        population[i].fitness = fitness(population[i].cost, minScore);
    }
}

void sortPopulation() {
    sort(population, population + nPop, byFitness());
}

void select() {
    // float TOTAL_SCORE = 0;

    // for (int i = 0; i < nPop; i++) {
    //     TOTAL_SCORE += population[i].fitness;
    // }

    // float selectionWheel[500];

    // selectionWheel[0] = 0;
    // for (int i = 1; i <= nPop; i++) {
    //     selectionWheel[i] = population[i - 1].fitness / TOTAL_SCORE;
    // }

    // for (int i = 1; i <= nPop; i++) {
    //     selectionWheel[i] += selectionWheel[i - 1];
    // }

    // int oldPop = nPop;

    // for (int i = 0; i < M_POP; i++) {
    //     float r = generateFloatRandom(0.000001, 1);

    //     int s;
    //     for (s = 0; s < oldPop; s++) {
    //         if (selectionWheel[s] < r && r <= selectionWheel[s + 1]) {
    //             break;
    //         }
    //     }

    //     population[i].copy(population[s]);

    // }

    nPop = M_POP;
}

void mutate() {
    for (int i = M_ET; i < nPop; i++) {
        population[i].mutateGreedy(G);
    }
    for (int i = nPop - 20; i < nPop; i++) {
        population[i].strongMutation();
    }

    // for (int i = 0; i < nPop; i++) {
    //     population[i].mutate();
    // }
}

void crossover() {
    for (int i = 0; i < nPop; i++) {
        population[i].randomProb();
    }

    sort(population, population + nPop, byProb());
    
    int oldPop = nPop;

    for (int i = 0; i < oldPop; i+=2) {

        if (population[i].rProb > M_C_PROB) {
            break;
        }

        Member *p[2];
        Member *c[2];

        p[0] = population + i;
        p[1] = population + i + 1;

        c[0] = population + nPop++;
        c[1] = population + nPop++;

        smallCrossover(p, c);

        // c[0] = population + nPop++;
        // c[1] = population + nPop++;

        // sSmallCrossover(p, c);

    }
}

void showPopulation() {
    for (int i = nPop - 1; i >= 0; i--) {
        showMember(population + i);
    }
}

void showMember(Member* m) {
    // printf("Gene: ");
    // for (int i = 0; i < m->n; i++) {
    //     printf("%d ", m->gene[i]);
    // }
    // printf("\n");

    printf("Cost: %d Fitness: %f\n", m->cost, m->fitness);

    // printf("\n");
}

int getMinScore() {
    int minScore = 999999;

    for (int i = 0; i < nPop; i++) {
        minScore = min(minScore, population[i].cost);
    }

    return minScore;
}

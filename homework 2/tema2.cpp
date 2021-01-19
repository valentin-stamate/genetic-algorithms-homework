#define _USE_MATH_DEFINES
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <random>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>

using namespace std;
using namespace std::chrono;


#define BITSTRING_MAX_LENGTH 1000
#define MAX_COMPONENTS 50

#define PRECIZION 5
#define GLOGAL_ITERATIONS 1
#define POPULATIONS_SIZE 100
#define MAX_GENERATIONS 1000
#define CROSSOVER_PROBABILITY 0.2
#define MUTATION_PROBABILITY 0.001

int COMPONENTS;
int FUNCTION;

ofstream out;

int BITPOINT_LEN, BITSTRING_LENGTH;
// Codul reprezinta o generalizare. Pentru rezultate, am spart codul si am rulat in paralele functiile cu diferite imputuri : 2, 5, 10, 15, 30

//random_device device;
mt19937 generator(0);

int generateRandom(int, int);
float generateFloatRandom(float, float);
long long getCurrentTime();

float getFitness(int*, float);

class Member {
    private:
    
    public:
    int gene[BITSTRING_MAX_LENGTH];
    float fitness;
    float randomProb;

    Member() {}

    Member(const Member &m) {
        for (int i = 0; i < BITSTRING_LENGTH; i++) {
            this->gene[i] = m.gene[i];
        }
        this->fitness = m.fitness;
        this->randomProb = m.randomProb;
    }

    void mutate() {
        for (int i = 0; i < BITSTRING_LENGTH; i++) {
            float r = generateFloatRandom(0, 1);
            if (r < MUTATION_PROBABILITY) {
                this->gene[i] = 1 - this->gene[i];
            }
        }  
    }

    void assignRandomProb() {
        this->randomProb = generateFloatRandom(0, 1);
    }
};

struct byProb {
    bool operator() (const Member &a, const Member &b) {
        return a.randomProb < b.randomProb;
    }
};
struct byFitness {
    bool operator() (const Member &a, const Member &b) {
        return a.fitness < b.fitness;
    }
};

void showPopulation(Member *, int&);
void toNormalTime(long long int);
void bitStringToPoint(int*, float*);
// setup
pair<float, float> interval;

pair<float, float> deJongInt(-5.12, 5.12);
float DeJong(float* X) {
    float val = 0;

    for (int i = 0; i < COMPONENTS; i++) {
        val += (X[i] * X[i]);
    }

    return val;
}

pair<float, float> schwefelInt(-500, 500);
float Schwefel(float* X) {
    float val = 0;

    for (int i = 0; i < COMPONENTS; i++) {
        val += ( - X[i] * sin(sqrt(abs(X[i]))) );
    }

    return val;
}

pair<float, float> rastriginInt(-5.12, 5.12);
float Rastrigin(float* X) {
    float val = 10 * COMPONENTS;
    for (int i = 0; i < COMPONENTS; i++) {
        float part = X[i] * X[i] - 10.0 * cos(2.0 * M_PI * X[i]);
        val += part;
    }
    return val;
}

pair<float, float> michalewiczInt(0, 3.14);
float Michalewicz(float* X) {
    float val = 0;
    for (int i = 0; i < COMPONENTS; i++) {
        val -= (sin(X[i]) * pow(sin((i * X[i] * X[i]) / 3.1415), 20));
    }
    return val;
}

// end setup

float (*Function)(float*);

void init();
void generateBitString(int*);
void printBitString(int*);
float getPointFromBitString(int*, int);
float getPointValue(int*);
void getNeighbour(int*, int*, int, int);
void showPoint(int*);

float simulateEvolution();
void SmallCrossOver(const Member, const Member, const Member[2]);
void CrossOver(Member*, int&);
void generatePopulation(Member*, int&);
void Select(Member*, int&);
void Evaluate(Member*, int&);
void Mutate(Member*, int&);

int main() {
    printf("\nEnter a number from 1 to 4 then press enter, same for dimensions.\n");
    cout<<"Function: 1: DeJong, 2: Swefel, 3: Rastrigin, 4: Michalewich\n";
    cin>>FUNCTION;

    char fileName[30];
    strcpy(fileName, "output_");
    fileName[7] = '0' + FUNCTION;
    fileName[8] = '_';

    cout<<"Give the nr of dimensions: ";
    cin>>COMPONENTS;
    
    if (COMPONENTS >= 10) {
        fileName[9] = '0' + COMPONENTS / 10;
    }
    else {
        fileName[9] = '0' + COMPONENTS;
    }
    fileName[10] = '\0';

    out.open(fileName);

    // FUNCTION = 3;
    // COMPONENTS = 5;

    init();

    printf("Component number %d\n", COMPONENTS);
    long long int t1 = getCurrentTime();
    for (int i = 0; i < GLOGAL_ITERATIONS; i++) {
        printf("Iteration %d\n", i);
        float result= simulateEvolution();
        out<<result<<" ";
    }
    out<<"\n";
    long long int t2 = getCurrentTime();
    toNormalTime(t2 - t1);

    return 0;
}


float simulateEvolution() {

    Member population[200];
    int n = 0;
    generatePopulation(population, n);

    Evaluate(population, n);

    for (int generation = 1; generation <= MAX_GENERATIONS; generation++) {
        // printf("Generation: %d\n", generation);
        Select(population, n);

        Mutate(population, n);

        CrossOver(population, n);

        // sort(population, population + n, byFitness());

        // float lowestFitnessValue = getPointValue(population[0].gene);

        // float highestFitnessValue = getPointValue(population[n - 1].gene);

        // printf("L: %f - H: %f\n", lowestFitnessValue, highestFitnessValue);

        Evaluate(population, n);

    }

    sort(population, population + n, byFitness());

    // float lowestFitnessValue = getPointValue(population[0].gene);

    float highestFitnessValue = getPointValue(population[n - 1].gene);

    // printf("L: %f - H: %f\n", lowestFitnessValue, highestFitnessValue);

    return highestFitnessValue;

}



void Select(Member *population, int &n) {
    float TOTAL_SCORE = 0;

    for (int i = 0; i < n; i++) {
        TOTAL_SCORE += population[i].fitness;
    }

    float selectionWheel[POPULATIONS_SIZE * 2];

    selectionWheel[0] = 0;
    for (int i = 1; i <= n; i++) {
        selectionWheel[i] = population[i - 1].fitness / TOTAL_SCORE;
    }

    for (int i = 1; i <= n; i++) {
        selectionWheel[i] += selectionWheel[i - 1];
    }

    Member oldPopulation[200];
    int oldN = n;


    for (int i = 0; i < oldN; i++) {
        oldPopulation[i] = population[i];
    }

    for (int i = 0; i < POPULATIONS_SIZE; i++) {
        float r = generateFloatRandom(0.000001, 1);

        int s;
        for (s = 0; s < n; s++) {
            if (selectionWheel[s] < r && r <= selectionWheel[s + 1]) {
                break;
            }
        }

        population[i] = oldPopulation[s];

    }

    n = POPULATIONS_SIZE;

}

void generatePopulation(Member *population, int &n) {
    for (int i = 0; i < POPULATIONS_SIZE; i++) {
        Member m = Member();
        generateBitString(m.gene);
        population[i] = m;
    }
    n = POPULATIONS_SIZE;
}

void SmallCrossOver(const Member a, const Member b, Member children[2]) {

    int pivot = generateRandom(1, BITSTRING_LENGTH - 2);
    
    for (int i = 0; i <= pivot; i++) {
        children[0].gene[i] = a.gene[i];
        children[1].gene[i] = b.gene[i];
    }

    for (int i = pivot + 1; i < BITSTRING_LENGTH; i++) {
        children[0].gene[i] = b.gene[i];
        children[1].gene[i] = a.gene[i];
    }
    
}

void CrossOver(Member *population, int &n) {
    for (int i = 0; i < n; i++) {
        population[i].assignRandomProb();
    }

    sort(population, population + n, byProb());
    
    int crossoverLine = 0;

    while (population[crossoverLine].randomProb < CROSSOVER_PROBABILITY && crossoverLine < n) {
        crossoverLine++;
    }
    crossoverLine--;

    if (crossoverLine % 2 == 1) {
        crossoverLine += (generateRandom(1, 10) % 2 ? 1 : -1);
    }

    for (int i = 0; i < crossoverLine; i+=2) {
        Member children[2];
        SmallCrossOver(population[i], population[i + 1], children);

        // children[0].mutate();
        // children[1].mutate();

        population[n++] = children[0];
        population[n++] = children[1];
    }
}

void Mutate(Member *population, int &n) {
    for (int i = 0; i < n; i++) {
        population[i].mutate();
    }
}

void Evaluate(Member *population, int &n) {

    float currentMin = 1000000000;

    for (int i = 0; i < n; i++) {
        float value = getPointValue(population[i].gene);
        
        if (currentMin > value) {
            currentMin = value;
        }
    }

    for (int i = 0; i < n; i++) {
        population[i].fitness = getFitness(population[i].gene, currentMin);
    }
}

//
void init() {

    switch (FUNCTION) {
    case 1:
        Function = DeJong;
        interval = deJongInt;
        break;
    case 2:
        Function = Schwefel;
        interval = schwefelInt;
        break;
    case 3:
        Function = Rastrigin;
        interval = rastriginInt;
        break;
    case 4:
        Function = Michalewicz;
        interval = michalewiczInt;
        break;
    default:
        break;
    }

    BITPOINT_LEN = ceil( log2(pow(10, PRECIZION) * (interval.second - interval.first)));
    BITSTRING_LENGTH = COMPONENTS * BITPOINT_LEN;
}

float getPointFromBitString(int* bitString, int position) {
    float n = 0;

    for (int i = BITPOINT_LEN * (position + 1) - 1; i >= BITPOINT_LEN * position; i--) {
        n = n * 2 + bitString[i];
    }

    float S = n / (2 << (BITPOINT_LEN - 1));
    S *= (interval.second - interval.first);
    S += interval.first;
    
    return S;
}

void generateBitString(int* bitString) {

    for (int i = 0; i < BITSTRING_LENGTH; i++) {
        int nr = generateRandom(1, 100);
        bitString[i] = nr % 2;
    }
}

void printBitString(int* bitString) {
    for (int i = 0; i < BITSTRING_LENGTH; i++) {
        if (i % BITPOINT_LEN == 0)
            cout<<" ";
        cout<<bitString[i]<<" ";
    }
    cout<<"\n";
}

float getPointValue(int *bitString) {
    
    float X[MAX_COMPONENTS];
    bitStringToPoint(bitString, X);

    // for (int i = 0; i < COMPONENTS; i++) {
    //     cout<<X[i]<<" ";
    // }
    // cout<<"\n";

    return Function(X);
}

float getFitness(int* bitString, float currentMin) {
    float fit = getPointValue(bitString);

    if (currentMin < 0) {
        fit += (-currentMin);
    }

    fit += 100;

    fit = 1 / fit;

    return fit;
}

int generateRandom(int x, int y) {
    uniform_int_distribution<int> distribution(x, y); // [x, y]
    long long ct = getCurrentTime() % 10000;
    int rn = distribution(generator) * ct;
    
    int dif = y - x;
    return (rn % dif) + x;
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

void toNormalTime(long long int time) {
    int hours = time / 3600000;
    time -= hours * 3600000;

    int min = time / 60000;
    time -= min * 60000;

    int sec = time / 1000;
    int msec = time % 1000;

    out<<hours<<"h "<<min<<"min "<<sec<<"s "<<msec<<"ms\n";
}

long long getCurrentTime() {
    milliseconds ms = duration_cast< milliseconds >(
        system_clock::now().time_since_epoch()
    );
    return ms.count();
}

void showPopulation(Member *population, int &n) {

    for (int i = 0; i < n; i++) {
        cout<<"This is "<<i<<"'th member\n";
        // printBitString(population[i].gene);
        showPoint(population[i].gene);
        cout<<"Function value is: "<<getPointValue(population[i].gene)<<"\n";
        cout<<"Function fitness: "<<population[i].fitness<<"\n";
        cout<<"\n";
    }

}

void bitStringToPoint(int* bitString, float* X) {

    for (int i = 0; i < COMPONENTS; i++) {
        X[i] = 0;
    }

    for (int i = BITSTRING_LENGTH - 1; i >= 0; i--) {
        X[i / BITPOINT_LEN] = X[i / BITPOINT_LEN] * 2 + bitString[i];
    }

    for (int i = 0; i < COMPONENTS; i++) {
        float S = X[i] / (2 << (BITPOINT_LEN - 1));
        S *= (interval.second - interval.first);
        S += interval.first;

        X[i] = S;
    }
}

void showPoint(int* bitString) {
    float X[MAX_COMPONENTS];

    bitStringToPoint(bitString, X);

    for (int i = 0; i < COMPONENTS; i++) {
        cout<<X[i]<<" ";
    }   
    cout<<"\n";

}

#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <algorithm>
#include <math.h>
#include <fstream>

using namespace std;
using namespace std::chrono;

#define LMAX 1000
#define NMAX 50

// Codul reprezinta o generalizare. Pentru rezultate, am spart codul si am rulat in paralele functiile cu diferite imputuri : 2, 5, 10, 15, 30
ofstream out("output.out");

struct result {
    float minimumValue;
    float X[NMAX];
};

const int pre = 5;
const int GlobalIterations = 30;
const int InnerIterations = 100000; //

int generateRandom(int, int);
float generateFloatRandom(float, float);
long long getCurrentTime();

//random_device device;
mt19937 generator(0);


// setup
int N = 30;
pair<float, float> interval;

pair<float, float> deJongInt(-5.12, 5.12);
float DeJong(float* X) {
    float val = 0;

    for (int i = 0; i < N; i++) {
        val += (X[i] * X[i]);
    }

    return val;
}

pair<float, float> schwefelInt(-500, 500);
float Schwefel(float* X) {
    float val = 0;

    for (int i = 0; i < N; i++) {
        val += ( - X[i] * sin(sqrt(abs(X[i]))) );
    }

    return val;
}

pair<float, float> rastriginInt(-5.12, 5.12);
float Rastrigin(float* X) {
    float val = 10 * N;
    for (int i = 0; i < N; i++) {
        val += (X[i] * X[i] - 10 * cos(2.0 * 3.1415 * X[i]));
    }
    return val;
}

pair<float, float> michalewiczInt(0, 3.14);
float Michalewicz(float* X) {
    float val = 0;
    for (int i = 0; i < N; i++) {
        val -= (sin(X[i]) * pow(sin((i * X[i] * X[i]) / 3.1415), 20));
    }
    return val;
}

// end setup

float (*Function)(float*);

int l, L;

void init() {
    l = ceil( log2(pow(10, pre) * (interval.second - interval.first)));
    L = N * l;
}

void generateBitString(int*);
void printBitString(int*);
float getPointFromBitString(int*, int);
void getNeighbour(int*, int*, int, int);
float fitness(int*);
result iteratedHillClimbing(bool);
result iteratedSimulatedAnnealing();

void HillClimbingFI() {
    cout<<"Starting Hill Climbing - FirstImprovement\n";
    out<<"Hill Climbing - FirstImprovement\n";

    int t1 = getCurrentTime();

    result RESULT = iteratedHillClimbing(true);
    out<<RESULT.minimumValue<<" ";
    for (int i = 1; i < GlobalIterations; i++) {
        result newRes = iteratedHillClimbing(true);

        out<<newRes.minimumValue<<" ";

        if (newRes.minimumValue < RESULT.minimumValue) {
            RESULT.minimumValue = newRes.minimumValue;
            
            // for (int j = 0; j < N; j++) {
            //     RESULT.X[j] = newRes.X[j];
            // }
        }
    }

    // cout<<"Minimum value found: "<<RESULT.minimumValue<<"\nPoint: X = ";
    // for (int i = 0; i < N; i++) {
    //     cout<<RESULT.X[i]<<" ";
    // }

    cout<<"\n";
    int t2 = getCurrentTime();
    
    out<<"\n"<<(t2 - t1)<<"\n\n";
    cout<<"Time(milis): "<<t2 - t1<<"\n\n\n";

} 

void HillClimbingBI() {
    cout<<"Starting Hill Climbing - BestImprovement\n";
    out<<"Hill Climbing - BestImprovement\n";

    int t1 = getCurrentTime();
    result RESULT = iteratedHillClimbing(false);
    out<<RESULT.minimumValue<<" ";
    for (int i = 1; i < GlobalIterations; i++) {
        result newRes = iteratedHillClimbing(false);
        out<<newRes.minimumValue<<" ";
        if (newRes.minimumValue < RESULT.minimumValue) {
            RESULT = newRes;
            // for (int j = 0; j < N; j++) {
            //     RESULT.X[j] = newRes.X[j];
            // }
        }
    }  
    
    // cout<<"Minimum value found: "<<RESULT.minimumValue<<"\nPoint: X = ";
    // for (int i = 0; i < N; i++) {
    //     cout<<RESULT.X[i]<<" ";
    // }
    cout<<"\n";
    int t2 = getCurrentTime();

    out<<"\n"<<(t2 - t1)<<"\n\n";
    cout<<"Time(milis): "<<t2 - t1<<"\n\n\n";
}  

void SimulatedAnnealing() {
    cout<<"Starting Simulated Annealing\n";
    out<<"Simulated Annealing\n";

    int t1 = getCurrentTime();
    result RESULT = iteratedSimulatedAnnealing();
    out<<RESULT.minimumValue<<" ";
    for (int i = 1; i < GlobalIterations; i++) {
        result newRes = iteratedSimulatedAnnealing();
        out<<newRes.minimumValue<<" ";
        if (newRes.minimumValue < RESULT.minimumValue) {
            RESULT = newRes;
            // for (int j = 0; j < N; j++) {
            //     RESULT.X[j] = newRes.X[j];
            // }
        }
    }

    // cout<<"Minimum value found: "<<RESULT.minimumValue<<"\nPoint: ";
    // for (int i = 0; i < N; i++) {
    //     cout<<RESULT.X[i]<<" ";
    // }
    cout<<"\n";
    int t2 = getCurrentTime();

    out<<"\n"<<(t2 - t1)<<"\n\n";
    cout<<"Time(milis): "<<t2 - t1;
    cout<<"\n";
}


void DeJongsProcess() {
    Function = &DeJong;
    interval = deJongInt;
    init();

    cout<<"De Jong's Function\n";
    cout<<"Minimum value 0 in point x_i = 0 i=1:n \n";
    
    out<<"De Jong's Function: minimum value 0\n";
    HillClimbingFI();

    HillClimbingBI();
    
    SimulatedAnnealing();

    out<<"\n\n";
    cout<<"\n\n\n";
}

void SchwefelProcess() {
    Function = &Schwefel;
    interval = schwefelInt;
    init();

    cout<<"Schwefel's Function\n";
    cout<<"Minimum value "<<- N * 418.9829<<" in point x_i = 420.9687 i=1:n \n";

    out<<"Schwefel's Function: minimum value "<<- N * 418.9829<<"\n";
    HillClimbingFI();

    HillClimbingBI();
    
    SimulatedAnnealing();

    out<<"\n\n";
    cout<<"\n\n\n";
}

void RastriginProcess() {
    Function = &Rastrigin;
    interval = rastriginInt;
    init();

    cout<<"Rastrigin's Function\n";
    cout<<"Minimum value 0 in point x_i = 0 i=1:n \n";

    out<<"Rastrigin's Function: minimum value 0\n";
    HillClimbingFI();

    HillClimbingBI();
    
    SimulatedAnnealing();

    out<<"\n\n";
    cout<<"\n\n\n";
}

void MichalewiczProcess() {
    Function = &Michalewicz;
    interval = michalewiczInt;
    init();

    cout<<"Michalewicz's Function\n";
    cout<<"Minimum value ?? in point x_i = ?? i=1:n \n\n";

    out<<"Michalewicz's Function: minimum value ??\n";
    HillClimbingFI();

    HillClimbingBI();
    
    SimulatedAnnealing();

    cout<<"\n\n\n";
}

int main() {
    cout<<"Components: "<<N<<"\n\n";
    out<<"Components: "<<N<<"\n";

    DeJongsProcess();

    SchwefelProcess();

    RastriginProcess();

    MichalewiczProcess();

    return 0;
}

result iteratedHillClimbing(bool firstImprovement) {
    const int ITERATIONS = InnerIterations;

    int currentSolution[LMAX];
    generateBitString(currentSolution);
    float finalValue = fitness(currentSolution);
    int finalPoint[LMAX];


    int neighbour[LMAX];
    for (int i = 0; i < ITERATIONS; i++) {
        
        generateBitString(currentSolution);
        float currentValue = fitness(currentSolution);

        bool done = false;

        while (!done) {

            bool found = false;
            for (int i = 0; i < N * 2; i++) {

                getNeighbour(currentSolution, neighbour, i / 2, i % 2);
                
                float newValue = fitness(neighbour);
                if (newValue < currentValue) {

                    currentValue = newValue;

                    for (int j = 0; j < L; j++) {
                        currentSolution[j] = neighbour[j];
                    }

                    found = true;
                    if (firstImprovement) {
                        break;
                    }

                }

            }

            done = !found;
        }

        if (currentValue < finalValue) {
            finalValue = currentValue;
            for (int i = 0; i < L; i++) {
                finalPoint[i] = currentSolution[i];
            }
        }
        
    }

    
    result res;
    res.minimumValue = finalValue;
    for (int i = 0; i < N; i++) {
       // res.X[i] = getPointFromBitString(currentSolution, i);
    }

    return res;

}

result iteratedSimulatedAnnealing() {
    float t = 100;
    const int ITERATIONS = InnerIterations;
    int currentSolution[LMAX];

    generateBitString(currentSolution);
    float totalValue = fitness(currentSolution);

    for (int i = 0; i < ITERATIONS; i++) {
        generateBitString(currentSolution);
    
        float currentValue = fitness(currentSolution);

        int neighbour[LMAX];

        bool done = false;
        while (t > 10e-8) {

            for (int i = 0; i < L; i++) {
                neighbour[i] = currentSolution[i];
            }
            
            int randomPosition = generateRandom(0, L - 1);
            neighbour[randomPosition] = 1 - neighbour[randomPosition];

            float newValue = fitness(neighbour);

            if (newValue < currentValue) {
                for (int j = 0; j < L; j++) {
                    currentSolution[j] = neighbour[j];
                }
                currentValue = newValue;
            } else if (generateFloatRandom(0, 1 - 0.00001) < exp( -abs(currentValue - newValue) / t ) ) {
                
                for (int j = 0; j < L; j++) {
                    currentSolution[j] = neighbour[j];
                }
                currentValue = newValue;
            }

            t *= 0.99;
        }

        if (currentValue < totalValue) {
            totalValue = currentValue;
        }

    }

    result res;
    res.minimumValue = totalValue;
    // for (int i = 0; i < N; i++) {
        // res.X[i] = getPointFromBitString(currentSolution, i);
    //}

    return res;

}

void getNeighbour(int* point, int* neighbour, int position, int type) {

    for (int i = 0; i < L; i++) {
        neighbour[i] = point[i];
    }

    if (type == 1) { // incremenet

        int i;
        for (i = (position + 1) * l - 1; i >= position * l && neighbour[i] == 1; i--) {
            neighbour[i] = 0;
        }

        neighbour[i] = 1;

    } else { // decrement
        int i;
        for (i = (position + 1) * l - 1; i >= position * l && neighbour[i] == 0; i--) {}
        neighbour[i] = 0;
        
        for (int j = i + 1; j <= (position + 1) * l - 1; j++) {
            neighbour[j] = 1;
        }

    }
}

float getPointFromBitString(int* bitString, int position) {
    float n = 0;

    for (int i = l * (position + 1) - 1; i >= l * position; i--) {
        n = n * 2 + bitString[i];
    }

    float S = n / (2 << (l - 1));
    S *= (interval.second - interval.first);
    S += interval.first;
    
    return S;
}

void generateBitString(int* bitString) {

    for (int i = 0; i < L; i++) {
        int nr = generateRandom(1, 10);
        bitString[i] = nr % 2;
    }
}

void printBitString(int* bitString) {
    for (int i = 0; i < L; i++) {
        if (i % l == 0)
            cout<<" ";
        cout<<bitString[i]<<" ";
    }
    cout<<"\n";
}

float fitness(int* bitString) {
    float X[NMAX];
    for (int i = 0; i < N; i++) {
        X[i] = 0;
    }

    for (int i = L - 1; i >= 0; i--) {
        X[i / l] = X[i / l] * 2 + bitString[i];
    }

    for (int i = 0; i < N; i++) {
        float S = X[i] / (2 << (l - 1));
        S *= (interval.second - interval.first);
        S += interval.first;

        X[i] = S;
    }

    return Function(X);
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
    return distribution(generator);
}

long long getCurrentTime() {
    milliseconds ms = duration_cast< milliseconds >(
        system_clock::now().time_since_epoch()
    );
    return ms.count();
}


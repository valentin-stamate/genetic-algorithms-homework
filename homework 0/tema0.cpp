#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <algorithm>
using namespace std;
using namespace std::chrono;

float eps = 0.01;
float gRandom(float, float);
long long getCurrentTime();

//random_device device;
mt19937 generator(0);

int n = 5;
pair<float, float> interval[30];
// Beale Function
float fun(float val[30]) {
    float p1 = 1.5 - val[0] + val[0] * val[2];
    float p2 = 2.25 - val[0] + val[0] * val[1] * val[1];
    float p3 = 2.625 - val[0] + val[0] * val[1] * val[1] * val[1];
    return p1 * p1 + p2 * p2 + p3 * p3;
}


void determ(int, float *, float *, float *);
void heuris();

void solve();

int main() {

    interval[0] = pair<float, float>(-4.5, 4.5);
    interval[1] = pair<float, float>(-4.5, 4.5);



    cout<<"Real result: 0 and x = (3, 0.5)\n\n";
    
    solve();

}

void solve() {

    long long start, end;
    start = getCurrentTime();

    float X[30], XSOL[30];
    for (int i = 0; i < n; i++) {
        X[i] = gRandom(interval[i].first, interval[i].second);
    }
    float fMin = fun(X);

    cout<<"Starting deterministic method:\n";
    
    determ(0, X, &fMin, XSOL);

    cout<<"Algorithm result: "<<fMin<<" for point X = (";
    for (int i = 0; i < n - 1; i++) {
        cout<<XSOL[i]<<", ";
    }
    cout<<XSOL[n - 1]<<")\n";

    end = getCurrentTime();

    cout<<"Time(milis): "<<(end - start);


    cout<<"\n\n";


    start = getCurrentTime();
    heuris();
    end = getCurrentTime();

    cout<<"Time(milis): "<<(end - start);

}

void determ(int depth, float *X, float *fMin, float *XSOL) {
    if (depth == n) {
        float Y = fun(X);
        if (*fMin > Y) {
            *fMin = Y;
            for (int j = 0; j < n; j++) {
                XSOL[j] = X[j];
            }
        }
        return;
    }

    for (float i = interval[depth].first; i <= interval[depth].second; i+=eps) {
        X[depth] = i;
        determ(depth + 1, X, fMin, XSOL);
    }

}

void heuris() {
    cout<<"Starting heuristic method:\n";

    float fMin;

    float X[30];
    float XSOL[30];
    for (int i = 0; i < n; i++) {
        
        X[i] = gRandom(interval[i].first, interval[i].second);
    }

    fMin = fun(X);

    int steps = 1000;

    while (--steps != 0) {

        for (int i = 0; i < n; i++) {
            X[i] = gRandom(interval[i].first, interval[i].second);

            float Y = fun(X);
            if (fMin > Y) {
                fMin = Y;
                for (int j = 0; j < n; j++) {
                    XSOL[j] = X[j];
                }
            }

        }

    }
    
    cout<<"Algorithm result: "<<fMin<<" for point X = (";
    for (int i = 0; i < n - 1; i++) {
        cout<<XSOL[i]<<", ";
    }
    cout<<XSOL[n - 1]<<")\n";
    cout<<")\n";

}

float gRandom(float x, float y) {
    uniform_real_distribution<float> distribution(x, y);
    return (int)(distribution(generator) * 100) * 1.0 / 100;
}

long long getCurrentTime() {
    milliseconds ms = duration_cast< milliseconds >(
        system_clock::now().time_since_epoch()
    );
    return ms.count();
}


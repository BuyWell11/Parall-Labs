#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#define DT 0.05
#define E 1e-4
#define THREAD_COUNT 8

int bodiesCount, timeSteps;
double *masses, GravConstant;

class Point2D {
public:
    double x;
    double y;

    Point2D() : x(0.0), y(0.0) {}

    Point2D(double _x, double _y) : x(_x), y(_y) {}

    void setNewData(double _x, double _y){
        this->x = _x;
        this->y = _y;
    }
};

class Body {
public:
    Point2D point;
    double m;
    int step = 0;
    int number;
    std::vector<double> forces;
    Point2D speed;
    Body(double x, double y, double Vx, double Vy, double m, int number) {
        this->point = Point2D(x, y);
        this->speed = Point2D(Vx, Vy);
        this->m = m;
        this->number = number;
        this->forces.resize(bodiesCount);
        pthread_rwlock_init(&rwlock, nullptr);
    }
    ~Body(){
        pthread_rwlock_destroy(&rwlock);
    }
private:
    pthread_rwlock_t rwlock;
};

std::vector<Body> bodies;

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}

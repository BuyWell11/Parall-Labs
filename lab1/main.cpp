#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#define DT 0.05
#define E 1e-4
#define THREAD_COUNT 8

int bodiesCount, timeSteps;
double GravConstant;

class Point2D {
public:
    double x;
    double y;

    Point2D() : x(0.0), y(0.0) {}

    Point2D(double _x, double _y) : x(_x), y(_y) {}

    void setNewData(double _x, double _y) {
        this->x = _x;
        this->y = _y;
    }
};

Point2D addVectors(Point2D vector1, Point2D vector2) {
    Point2D resultVector = {vector1.x + vector2.x, vector1.y + vector2.y};
    return resultVector;
}
        
Point2D subtractVectors(Point2D vector1, Point2D vector2) {
    Point2D resultVector = {vector1.x - vector2.x, vector1.y - vector2.y};
    return resultVector;
}
        
Point2D scaleVector(double constant, Point2D vector) {
    Point2D resultVector = {vector.x * constant, vector.y * constant};
    return resultVector;
}
        
double mod(Point2D vector) {
    return sqrt(vector.x * vector.x + vector.y * vector.y);
}

class Body {
public:
    Point2D point;
    Point2D force;
    Point2D speed;
    double m;
    int step = 0;
    int number;
    Body(double x, double y, double Vx, double Vy, double m, int number) {
        this->point = Point2D(x, y);
        this->speed = Point2D(Vx, Vy);
        this->m = m;
        this->number = number;
        this->force = Point2D(0,0);
        pthread_rwlock_init(&rwlock, nullptr);
    }
    
    ~Body(){
        pthread_rwlock_destroy(&rwlock);
    }
    
    
    void calculateForce(Body body) {
        double denominator = pow(mod(subtractVectors(this->point, body.point)), 3);
        
        if (denominator < E) {
            denominator = E;
        }
        
        this->force = addVectors(this->force, scaleVector(GravConstant * body.m / denominator, subtractVectors(body.point, this->point)));
    }
    
    void calculatePosition() {
        this->point = addVectors(this->point, scaleVector(DT, this->speed));
    }
    
    void calculateSpeed() {
        this->speed = addVectors(this->speed, scaleVector(DT, this->force));
    }
    
private:
    pthread_rwlock_t rwlock;
};

std::vector<Body> bodies;

bool initiateSystem(const std::string &filename) {
    std::ifstream inputFile(filename);
    double m, x, y, Vx, Vy;

    if (!inputFile) {
        std::cout << "Could not open the file" << std::endl;
        return false;
    }

    inputFile >> GravConstant >> bodiesCount >> timeSteps;
    for (int i = 0; i < bodiesCount; i++) {
        inputFile >> m;
        inputFile >> x >> y;
        inputFile >> Vx >> Vy;
        bodies.emplace_back(x, y, Vx, Vy, m, i);
    }

    if (inputFile.fail()) {
        std::cerr << "Could not read data from the file." << std::endl;
        return false;
    }

    inputFile.close();
    return true;
}

int main(int argC, char *argV[]) {
    pthread_t threads[THREAD_COUNT];
    if (argC != 2)
        printf("Usage : %s <file name containing system configuration data>", argV[0]);
    else {
        if (initiateSystem(argV[1])) {
            std::cout << "Body   :     x              y           vx              vy   " << std::endl;
            if(bodies.size() <= THREAD_COUNT){
                for(int i = 0; i < bodies.size(); i++)
            }
        }
        return 0;
    }
}

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

    static Point2D addVectors(Point2D vector1, Point2D vector2){
        return Point2D(vector1.x + vector2.x, vector1.y + vector2.y);
    }

    static Point2D subtractVectors(Point2D vector1, Point2D vector2){
        return Point2D(vector1.x - vector2.x, vector1.y - vector2.y);
    }

    static Point2D scaleVector(double constant, Point2D vector){
        return Point2D(vector.x * constant, vector.y * constant);
    }

    static double mod(Point2D vector){
        return sqrt(vector.x * vector.x + vector.y * vector.y);
    }
};

class Body {
public:
    Point2D point;
    Point2D force;
    Point2D speed;
    std::vector<Point2D> forces;
    double m;
    int step = 0;
    int number;
    Body(double x, double y, double Vx, double Vy, double m, int number) {
        this->point = Point2D(x, y);
        this->speed = Point2D(Vx, Vy);
        this->m = m;
        this->number = number;
        forces.resize(bodiesCount);

        pthread_rwlock_init(&rwlock, nullptr);
    }
    
    void calculateForce(Body& body) {
        double denominator = pow(Point2D::mod(Point2D::subtractVectors(getPoint(), body.getPoint())), 3);
        
        if (denominator < E) {
            denominator = E;
        }
        
        forces[body.number] = Point2D::addVectors(forces[this->number], Point2D::scaleVector(GravConstant * body.m / denominator, Point2D::subtractVectors(body.getPoint(), getPoint())));
    }
    
    void calculateForceSum() {
        for (int i = 0; i < bodiesCount; i++) {
            this->force = Point2D::addVectors(this->force, this->forces[i]);
        }
    }
    
    void calculatePosition() {
        this->point = Point2D::addVectors(getPoint(), Point2D::scaleVector(DT, this->speed));
    }
    
    void calculateSpeed() {
        this->speed = Point2D::addVectors(this->speed, Point2D::scaleVector(DT, this->force));
    }

    Point2D getPoint(){
        return this->point;
    }

    Point2D getForce(int index){
        return forces[index];
    }

    ~Body(){
        pthread_rwlock_destroy(&rwlock);
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

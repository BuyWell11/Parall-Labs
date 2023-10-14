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

    void addVector(Point2D& vector){
        this->x += vector.x;
        this->y += vector.y;
    }

    void subtractVector(Point2D& vector){
        this->x -= vector.x;
        this->y -= vector.y;
    }

    void scaleVector(double constant){
        this->x *= constant;
        this->y *= constant;
    }

    double mod(Point2D& vector){
        return sqrt(this->x * vector.x + this->y * vector.y);
    }

    void copyRevers(Point2D point) {
        this->x = point.x * -1;
        this->y = point.y * -1;
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
    Body(double x, double y, double Vx, double Vy, double m, int number) : point(Point2D(x, y)), speed(Point2D(Vx, Vy)), m(m), number(number) {
        forces.resize(bodiesCount);

        pthread_rwlock_init(&rwlock, nullptr);
    }
    
    void calculateForce(Body& body) {
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

    Point2D getPoint() {
        return this->point;
    }


    Point2D getForce(int index){
        return forces[index];
    }

    void copyForce(Body body) {
        this->forces[body.number].copyRevers(body.getForce(this->number));
    }

    bool isNotNullForce(int index) {
        if(forces[index].x != 0 || forces[index].y != 0){
            return true;
        }
        return false;
    }

    void writeForce(Body body, int index){
        //тут функция подсчёта как переданная точка действует на точку метод которой вызван и сохраняет в массив с forces[index]
        //по сути колкулейт форсес можешь сюда вствить, а сам метод удалить и добавить потом саммари метод который по массиву сил проходит и суммирует
        //ну и точка берётся только через метод гет поинт чтобы эффекта гонки не было
    }

    ~Body(){
        pthread_rwlock_destroy(&rwlock);
    }

private:
    pthread_rwlock_t rwlock;
};

struct ThreadData{
    std::vector<Body> partOfBodies;
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

void* ThreadFunction(void* arg){
    auto* data = static_cast<ThreadData*>(arg);

    for(Body& partBody: data->partOfBodies){
        for(Body& body: bodies){
            if(partBody.number == body.number){
                continue;
            }

            if(body.isNotNullForce(partBody.number)){
                partBody.copyForce(body);
            }
            else{
                
            }

        }
    }
    return NULL;
}

int main(int argC, char *argV[]) {
    pthread_t threads[THREAD_COUNT];


    if (argC != 2)
        printf("Usage : %s <file name containing system configuration data>", argV[0]);
    else {
        if (initiateSystem(argV[1])) {


            int num_threads_to_create = (bodiesCount < THREAD_COUNT) ? bodiesCount : THREAD_COUNT;
            int bodies_per_tread = bodiesCount / num_threads_to_create;

            std::vector<ThreadData> thread_data(num_threads_to_create);

            int start_index = 0;

            std::cout << "Body   :     x              y           vx              vy   " << std::endl;

            for(int i = 0; i < num_threads_to_create; i++){
                int slice_size = i < bodiesCount % num_threads_to_create ? bodies_per_tread + 1 : bodies_per_tread;
                int end_index = start_index + slice_size;

                thread_data[i].partOfBodies.assign(bodies.begin() + start_index, bodies.begin() + end_index);
            }
        }
        return 0;
    }
}

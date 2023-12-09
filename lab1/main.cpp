#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <pthread.h>
#include "timer.h"

#define DT 0.05
#define E 1e-4
#define THREAD_COUNT 8

pthread_barrier_t barrier_forces;
pthread_barrier_t barrier_calculate;
pthread_barrier_t barrier_thread;
int bodiesCount, timeSteps;
double GravConstant;


class PointForPrint {
public:
    double x;
    double y;
    double Vx;
    double Vy;

    PointForPrint(double _x, double _y, double _Vx, double _Vy) : x(_x), y(_y), Vx(_Vx), Vy(_Vy) {}
};

std::vector<std::vector<PointForPrint>> cycleVector;

class Point2D {
public:
    double x;
    double y;
    bool isNotNull;

    Point2D() : x(0), y(0), isNotNull(false) {}

    Point2D(double _x, double _y) : x(_x), y(_y), isNotNull(false) {}

    static Point2D addVectors(Point2D vector1, Point2D vector2) {
        Point2D point2D = Point2D(vector1.x + vector2.x, vector1.y + vector2.y);
        return point2D;
    }

    static Point2D subtractVectors(Point2D vector1, Point2D vector2) {
        Point2D point2D = Point2D(vector1.x - vector2.x, vector1.y - vector2.y);
        return point2D;
    }

    static Point2D scaleVector(double constant, Point2D vector) {
        Point2D point2D = Point2D(vector.x * constant, vector.y * constant);
        return point2D;
    }

    static double mod(Point2D vector) {
        double temp = sqrt(vector.x * vector.x + vector.y * vector.y);
        return temp;
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

    Body(double x, double y, double Vx, double Vy, double m, int number) : point(Point2D(x, y)), speed(Point2D(Vx, Vy)),
                                                                           m(m), number(number) {
        forces.resize(bodiesCount);
    }

    void calculateForce(Body &body) {
        Point2D firstPoint = this->point;
        Point2D secondPoint = body.point;

        double denominator = pow(Point2D::mod(Point2D::subtractVectors(firstPoint, secondPoint)),3);

        if (denominator < E) {
            denominator = E;
        }

        forces[body.number] = Point2D::addVectors(forces[this->number],
                                                  Point2D::scaleVector(GravConstant * body.m / denominator,
                                                                       Point2D::subtractVectors(secondPoint,
                                                                                                firstPoint)));

        body.forces[this->number].copyRevers(forces[body.number]);
        forces[body.number].isNotNull = true;
    }

    void calculateForceSum() {
        for (int i = 0; i < bodiesCount; i++) {
            this->force = Point2D::addVectors(this->force, this->forces[i]);
            this->forces[i].isNotNull = false;
        }
    }

    void calculatePosition() {
        this->point = Point2D::addVectors(this->point, Point2D::scaleVector(DT, this->speed));
    }

    void calculateSpeed() {
        this->speed = Point2D::addVectors(this->speed, Point2D::scaleVector(DT, this->force));
    }
};

struct ThreadData {
    int start;
    int end;
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
    cycleVector.resize(bodiesCount);
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

void *ThreadFunction(void *arg) {
    auto *data = static_cast<ThreadData *>(arg);
    int start = data->start;
    int end = data->end;

    for (int j = 0; j < timeSteps; j++) {
        for (int i = start; i < end; i++) {
            Body &partBody = bodies[i];
            for (Body &body: bodies) {
                if (partBody.number < body.number) {
                    partBody.calculateForce(body);
                }
            }
        }

        pthread_barrier_wait(&barrier_forces);
        for (int i = start; i < end; i++) {
            Body &partBody = bodies[i];
            partBody.calculateForceSum();
            partBody.calculatePosition();
            partBody.calculateSpeed();
            partBody.step++;
            cycleVector[i].emplace_back(partBody.point.x, partBody.point.y, partBody.speed.x, partBody.speed.y);
        }
        pthread_barrier_wait(&barrier_calculate);
    }
    pthread_barrier_wait(&barrier_thread);
    return NULL;
}

void writeOutput() {
    std::ofstream outputFile("output.txt", std::ios::out);
    if (outputFile.is_open()) {
        for (int i = 0; i < timeSteps; i++) {
            outputFile << "Cycle " << i + 1 << std::endl;
            for (int j = 0; j < bodiesCount; j++) {
                outputFile << "Body " << j << " : " << cycleVector[j][i].x << "\t" << cycleVector[j][i].y << "\t"
                           << cycleVector[j][i].Vx << "\t" << cycleVector[j][i].Vy << std::endl;
            }
        }
        outputFile.close();
    } else {
        std::cerr << "Error while writing output" << std::endl;
    }
}

int main(int argC, char *argV[]) {
    double start, end;

    if (argC != 2)
        printf("Usage : %s <file name containing system configuration data>", argV[0]);
    else {
        if (initiateSystem(argV[1])) {
            pthread_t threads[THREAD_COUNT];

            int num_threads_to_create = (bodiesCount < THREAD_COUNT) ? bodiesCount : THREAD_COUNT;
            int bodies_per_tread = bodiesCount / num_threads_to_create;

            std::vector<ThreadData> thread_data(num_threads_to_create);

            int start_index = 0;

            for (int i = 0; i < num_threads_to_create; i++) {

                int slice_size = i < bodiesCount % num_threads_to_create ? bodies_per_tread + 1 : bodies_per_tread;
                int end_index = start_index + slice_size;
                thread_data[i].start = start_index;
                thread_data[i].end = end_index;
                start_index = end_index;
            }

            pthread_barrier_init(&barrier_forces, NULL, num_threads_to_create);
            pthread_barrier_init(&barrier_calculate, NULL, num_threads_to_create);
            pthread_barrier_init(&barrier_thread, NULL, num_threads_to_create);


            GET_TIME(start);

            for (int j = 0; j < num_threads_to_create; j++) {
                pthread_create(&threads[j], nullptr, ThreadFunction, &thread_data[j]);
            }

            for (int j = 0; j < num_threads_to_create; j++) {
                pthread_join(threads[j], nullptr);
            }

            GET_TIME(end);

            std::cout << "Time spent: " << end - start << " sec" << std::endl;
            writeOutput();
            pthread_barrier_destroy(&barrier_calculate);
            pthread_barrier_destroy(&barrier_forces);
            pthread_barrier_destroy(&barrier_thread);
        }
        return 0;
    }
}

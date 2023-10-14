#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <pthread.h>

#define DT 0.05
#define E 1e-4
#define THREAD_COUNT 8

int bodiesCount, timeSteps, curStep = 1;
double GravConstant;

class Point2D {
public:
    double x;
    double y;
    bool used;

    Point2D() : x(0), y(0), used(false) {}

    Point2D(double _x, double _y) : x(_x), y(_y), used(false) {}

    static Point2D addVectors(Point2D vector1, Point2D vector2) {
        return Point2D(vector1.x + vector2.x, vector1.y + vector2.y);
    }

    static Point2D subtractVectors(Point2D vector1, Point2D vector2) {
        return Point2D(vector1.x - vector2.x, vector1.y - vector2.y);
    }

    static Point2D scaleVector(double constant, Point2D vector) {
        return Point2D(vector.x * constant, vector.y * constant);
    }

    static double mod(Point2D vector) {
        return sqrt(vector.x * vector.x + vector.y * vector.y);
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

        pthread_rwlock_init(&rwlock, nullptr);
    }

    void calculateForce(Body &body) {
        pthread_rwlock_wrlock(&rwlock);
        double denominator = pow(Point2D::mod(Point2D::subtractVectors(getPoint(), body.getPoint())), 3);

        if (denominator < E) {
            denominator = E;
        }

        forces[body.number] = Point2D::addVectors(forces[this->number],
                                                  Point2D::scaleVector(GravConstant * body.m / denominator,
                                                                       Point2D::subtractVectors(body.getPoint(),
                                                                                                getPoint())));
        forces[body.number].used = true;
        pthread_rwlock_unlock(&rwlock);
    }

    void calculateForceSum() {
        pthread_rwlock_wrlock(&rwlock);
        for (int i = 0; i < bodiesCount; i++) {
            this->force = Point2D::addVectors(this->force, this->forces[i]);
        }
        for (Point2D& point: this->forces) {
            point.used = false;
        }
        this->forces.resize(bodiesCount);
        pthread_rwlock_unlock(&rwlock);
    }

    void calculatePosition() {
        this->point = Point2D::addVectors(getPoint(), Point2D::scaleVector(DT, this->speed));
    }

    void calculateSpeed() {
        this->speed = Point2D::addVectors(this->speed, Point2D::scaleVector(DT, this->force));
    }

    Point2D getPoint() {
        return this->point;
    }


    Point2D getForce(int index) {
        return forces[index];
    }


    void copyForceRevers(Body body) {
        pthread_rwlock_wrlock(&rwlock);
        this->forces[body.number].copyRevers(body.getForce(this->number));
        pthread_rwlock_unlock(&rwlock);
    }

    bool isNotNullForce(int index) {
        pthread_rwlock_rdlock(&rwlock);
        if (forces[index].used) {
            pthread_rwlock_unlock(&rwlock);
            return true;
        }
        pthread_rwlock_unlock(&rwlock);
        return false;
    }

    bool isAllForcesCalculated() {
        int countOfNullForces = 0;
        for (Point2D &f: forces) {
            if (f.x == 0 && f.y == 0) {
                countOfNullForces++;
                if (countOfNullForces > 1) {
                    return false;
                }
            }
        }
        return true;
    }

    void plusStep() {
        this->step++;
    }

    int getStep() {
        pthread_rwlock_rdlock(&rwlock);
        pthread_rwlock_unlock(&rwlock);
        return this->step;
    }

    ~Body() {
        pthread_rwlock_destroy(&rwlock);
    }

private:
    pthread_rwlock_t rwlock;
};

struct ThreadData {
    int start;
    int end;
    int thread_num;
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

void *ThreadFunction(void *arg) {
    auto *data = static_cast<ThreadData *>(arg);
    int start = data->start;
    int end = data->end;
    int thread_num = data->thread_num;
    std::cout << "Thread " << thread_num << " started" << std::endl;

    for (int i = start; i < end; i++) {
        Body &partBody = bodies[i];
        for (Body &body: bodies) {
            if (partBody.number == body.number) {
                std::cout << "Thread " << thread_num << " : same body " << partBody.number << " and " << body.number
                          << std::endl;
                continue;
            }

            if (body.isNotNullForce(partBody.number)) {
                std::cout << "Thread " << thread_num << " : copy force from " << body.number << " to "
                          << partBody.number << std::endl;
                partBody.copyForceRevers(body);
            } else {
                std::cout << "Thread " << thread_num << " : calculate force for " << partBody.number << " from "
                          << body.number << std::endl;
                partBody.calculateForce(body);
            }
        }
    }
    for (int i = start; i < end; i++) {
        Body &partBody = bodies[i];
        while (!partBody.isAllForcesCalculated()) {
            continue;
        }
        partBody.calculateForceSum();
        partBody.calculatePosition();
        partBody.calculateSpeed();
        partBody.plusStep();
    }
    return NULL;
}

bool checkStep() {
    for (Body &body: bodies) {
        if (body.getStep() != curStep) {
            return false;
        }
    }
    return true;
}

int main(int argC, char *argV[]) {

    if (argC != 2)
        printf("Usage : %s <file name containing system configuration data>", argV[0]);
    else {
        if (initiateSystem(argV[1])) {


            int num_threads_to_create = (bodiesCount < THREAD_COUNT) ? bodiesCount : THREAD_COUNT;
            int bodies_per_tread = bodiesCount / num_threads_to_create;

            std::vector<ThreadData> thread_data(num_threads_to_create);

            int start_index = 0;

            std::cout << "Body   :     x              y           vx              vy   " << std::endl;

            for (int i = 0; i < num_threads_to_create; i++) {

                int slice_size = i < bodiesCount % num_threads_to_create ? bodies_per_tread + 1 : bodies_per_tread;
                int end_index = start_index + slice_size;
                thread_data[i].start = start_index;
                thread_data[i].end = end_index;
                thread_data[i].thread_num = i;
                start_index = end_index;
            }

            for (int i = 0; i < timeSteps; i++) {
                pthread_t threads[THREAD_COUNT];

                for (int j = 0; j < num_threads_to_create; j++) {
                    pthread_create(&threads[j], nullptr, ThreadFunction, &thread_data[j]);
                }

                for (int j = 0; j < num_threads_to_create; j++) {
                    pthread_join(threads[j], nullptr);
                }
                while (!checkStep()) {
                }
                std::cout << "Cycle " << curStep << std::endl;
                for (Body &body: bodies) {
                    std::cout << "Body " << body.number << " : " << body.point.x << "\t" << body.point.y << "\t"
                              << body.speed.x << "\t" << body.speed.y << std::endl;
                }
                curStep++;
            }
        }
        return 0;
    }
}

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <pthread.h>
#include <unistd.h>

#define DT 0.05
#define E 1e-4
#define THREAD_COUNT 8

int bodiesCount, timeSteps, curStep = 1;
double GravConstant;

class Point2D {
public:
    double x;
    double y;
    bool isNotNull;

    Point2D() : x(0), y(0), isNotNull(false) {}

    Point2D(double _x, double _y) : x(_x), y(_y), isNotNull(false) {}

    static Point2D addVectors(Point2D vector1, Point2D vector2) {
        pthread_mutex_lock(&mutex);
        Point2D point2D = Point2D(vector1.x + vector2.x, vector1.y + vector2.y);
        pthread_mutex_unlock(&mutex);
        return point2D;
    }

    static Point2D subtractVectors(Point2D vector1, Point2D vector2) {
        pthread_mutex_lock(&mutex);
        Point2D point2D = Point2D(vector1.x - vector2.x, vector1.y - vector2.y);
        pthread_mutex_unlock(&mutex);
        return point2D;
    }

    static Point2D scaleVector(double constant, Point2D vector) {
        pthread_mutex_lock(&mutex);
        Point2D point2D = Point2D(vector.x * constant, vector.y * constant);
        pthread_mutex_unlock(&mutex);
        return point2D;
    }

    static double mod(Point2D vector) {
        pthread_mutex_lock(&mutex);
        double temp = sqrt(vector.x * vector.x + vector.y * vector.y);
        pthread_mutex_unlock(&mutex);
        return temp;
    }

    void copyRevers(Point2D point) {
        pthread_mutex_lock(&mutex);
        this->x = point.x * -1;
        this->y = point.y * -1;
        pthread_mutex_unlock(&mutex);
    }

private:
    static pthread_mutex_t mutex;
};

pthread_mutex_t Point2D::mutex = PTHREAD_MUTEX_INITIALIZER;

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
        Point2D firstPoint = this->getPoint();
        Point2D secondPoint = body.getPoint();
        double denominator = pow(Point2D::mod(Point2D::subtractVectors(firstPoint, secondPoint)), 3);

        if (denominator < E) {
            denominator = E;
        }

        forces[body.number] = Point2D::addVectors(forces[this->number],
                                                  Point2D::scaleVector(GravConstant * body.m / denominator,
                                                                       Point2D::subtractVectors(secondPoint,
                                                                                                firstPoint)));
        forces[body.number].isNotNull = true;
        pthread_rwlock_unlock(&rwlock);
    }

    void calculateForceSum() {
        pthread_rwlock_wrlock(&rwlock);
        for (int i = 0; i < bodiesCount; i++) {
            this->force = Point2D::addVectors(this->force, this->forces[i]);
            this->forces[i].isNotNull = false;
        }
        pthread_rwlock_unlock(&rwlock);
    }

    void calculatePosition() {
        pthread_rwlock_wrlock(&rwlock);
        this->point = Point2D::addVectors(this->point, Point2D::scaleVector(DT, this->speed));
        pthread_rwlock_unlock(&rwlock);
    }

    void calculateSpeed() {
        pthread_rwlock_wrlock(&rwlock);
        this->speed = Point2D::addVectors(this->speed, Point2D::scaleVector(DT, this->force));
        pthread_rwlock_unlock(&rwlock);
    }

    Point2D getPoint() {
        std::cout << "thread " << this->number << " read point " << std::endl;
        Point2D point2D = this->point;
        std::cout << "thread " << this->number << " stop reading point " << std::endl;
        return point2D;
    }


    Point2D getForce(int index) {
        pthread_rwlock_rdlock(&rwlock);
        Point2D point2D = forces[index];
        pthread_rwlock_unlock(&rwlock);
        return point2D;
    }


    void copyForceRevers(Body body) {
        pthread_rwlock_wrlock(&rwlock);
        this->forces[body.number].copyRevers(body.getForce(this->number));
        pthread_rwlock_unlock(&rwlock);
    }

    bool isNotNullForce(int index) {
        pthread_rwlock_rdlock(&rwlock);
        if (forces[index].isNotNull) {
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
        int temp = this->step;
        pthread_rwlock_unlock(&rwlock);
        return temp;
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
//    if(thread_num == 0){
//        usleep(1000);
//    }
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

            pthread_t threads[THREAD_COUNT];

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
                for (Body &body: bodies){
                    std::cout << body.number << std::endl;
                    for(Point2D point: body.forces){
                        std::cout << point.x << " " << point.y << " " << point.isNotNull << std::endl;
                    }
                }
                curStep++;
            }
        }
        return 0;
    }
}

#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>

void printMatrix(const std::vector<int>& matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << matrix[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(int* vector, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << vector[i] << " ";
    }
    std::cout << std::endl;
}

void readFile(const char* filename, int& rows, int& cols, std::vector<int>& matrix, int* vector) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error while oppening file " << filename << std::endl;
        exit(1);
    }

    file >> rows >> cols;

    matrix.resize(rows * cols);

    // Чтение матрицы
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file >> matrix[i * cols + j];
        }
    }

    // Чтение вектора
    for (int i = 0; i < cols; ++i) {
        file >> vector[i];
    }

    file.close();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 2) {
        if (rank == 0) {
            std::cerr << "Wrong file name: " << argv[0] << " <file name>\n";
        }
        MPI_Finalize();
        return 1;
    }

    const char *filename = argv[1];

    int rows, cols;
    std::vector<int> matrix;
    int* vector = new int[cols];

    if (rank == 0) {
        // Чтение данных только на корневом процессе
        readFile(filename, rows, cols, matrix, vector);

        if (cols % size != 0) {
            std::cerr << "Enter the number of processes that is a multiple of the number of matrix rows or columns";

            delete vector;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::cout << "Initial matrix:\n";
        printMatrix(matrix, rows, cols);

        std::cout << "Initial vector:\n";
        printVector(vector, cols);
    }

    // Рассылка данных о размерах матрицы и вектора всем процессам
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Выделение памяти на каждом процессе
    std::vector<int> local_matrix(rows * cols / size);
    std::vector<int> local_result(rows / size, 0);

    // Рассылка вектора всем процессам
    MPI_Bcast(vector, cols, MPI_INT, 0, MPI_COMM_WORLD);

    // Рассылка частей матрицы всем процессам
    MPI_Scatter(matrix.data(), rows * cols / size, MPI_INT,
                local_matrix.data(), rows * cols / size, MPI_INT, 0, MPI_COMM_WORLD);

    // Вывод данных после отправки
    std::cout << "Process " << rank << " after Scatter:\n";
    std::cout << "Local matrix:\n";
    printMatrix(local_matrix, rows / size, cols);

    std::cout << "Local vector:\n";
    printVector(vector, cols);

    // Локальное умножение части матрицы на вектор
    for (int i = 0; i < rows / size; ++i) {
        for (int j = 0; j < cols; ++j) {
            local_result[i] += local_matrix[i * cols + j] * vector[j];
        }
    }

    // Сбор результатов на корневом процессе
    int* result = new int[rows];
    MPI_Gather(local_result.data(), rows / size, MPI_INT,
               result, rows / size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Вывод результата на корневом процессе
        std::cout << "Result:\n";
        printVector(result, cols);
    }

    delete vector;
    delete result;
    MPI_Finalize();
    return 0;
}
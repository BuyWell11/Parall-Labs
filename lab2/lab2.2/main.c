#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

typedef struct {
    int       nb_proc;    // Количество процессов 
    MPI_Comm  grid_comm;  // Коммуникатор для сетки процессов
    MPI_Comm  row_comm;   // Коммуникатор для строк сетки процессов
    MPI_Comm  col_comm;   // Коммуникатор для столбцов сетки процессов
    int       order;      // Размерность сетки (m*n)
    int       pos_row;    // Позиция процесса в коммуникаторе для строк
    int       pos_col;    // Позиция процесса в коммуникаторе для столбцов
    int       grid_rank;  // Ранк сетки
} GRID;


void multiply_matrix(GRID* grid, float* block_A, float* block_B, float* block_C) {
	int i, j, k;
        int sqroot = (int)sqrt(grid->nb_proc);

    	for (i = 0; i < sqroot; i++)
        	for (j = 0; j < sqroot; j++)
            		for (k = 0; k < sqroot; k++)
                		block_C[i * sqroot + j] += block_A[i * sqroot + k] * block_B[k * sqroot + j];
}

void multiply_matrix_linear(int size, float* block_A, float* block_B, float* block_C) {
	int i, j, k;

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			for (k = 0; k < size; k++) {
					block_C[i * size + j] += block_A[i * size + k] * block_B[k * size + j];
			}
		}
	}
}

void cannon(GRID* grid, float* block_A, float* block_B, float* block_C) {
	int sqroot = sqrt(grid->nb_proc);
  	int shift_source, shift_dest;
	MPI_Status status;
	int up_rank, down_rank, left_rank, right_rank;
	int i;

	MPI_Cart_shift(grid->grid_comm, 1, -1, &right_rank, &left_rank); 
	MPI_Cart_shift(grid->grid_comm, 0, -1, &down_rank, &up_rank); 
	MPI_Cart_shift(grid->grid_comm, 1, -grid->pos_row, &shift_source, &shift_dest); 

	MPI_Sendrecv_replace(block_A, sqroot*sqroot, MPI_FLOAT, shift_dest, 1, shift_source, 1, grid->grid_comm, &status); 

	MPI_Cart_shift(grid->grid_comm, 0, -grid->pos_col, &shift_source, &shift_dest); 
	MPI_Sendrecv_replace(block_B, sqroot*sqroot, MPI_FLOAT, shift_dest, 1, shift_source, 1, grid->grid_comm, &status); 
   
    for (i=0; i<sqroot; i++) 
    { 
        multiply_matrix(grid, block_A, block_B, block_C); 
        MPI_Sendrecv_replace(block_A, grid->nb_proc, MPI_FLOAT, left_rank, 1, right_rank, 1, grid->grid_comm, &status); 
        MPI_Sendrecv_replace(block_B, grid->nb_proc, MPI_FLOAT, up_rank, 1, down_rank, 1, grid->grid_comm, &status); 
    } 
   
    MPI_Cart_shift(grid->grid_comm, 1, +grid->pos_row, &shift_source, &shift_dest); 
    MPI_Sendrecv_replace(block_B, grid->nb_proc, MPI_FLOAT, shift_dest, 1, shift_source, 1, grid->grid_comm, &status); 

    MPI_Cart_shift(grid->grid_comm, 0, +grid->pos_col, &shift_source, &shift_dest); 
    MPI_Sendrecv_replace(block_B, grid->nb_proc, MPI_FLOAT, shift_dest, 1, shift_source, 1, grid->grid_comm, &status); 	
}

void print_matrix(GRID* grid, float* mat) {
	if (grid->grid_rank == 0) {
		printf("Matrix: \n");

		int i, j;
		for (i = 0; i < sqrt(grid->order); i++) {
			for (j = 0; j < sqrt(grid->order); j++) {
				if(log10((int)mat[i * grid->order + j]) + 1 <= 2) {
					printf("%2d ", (int)mat[i * grid->order + j]);
				} else if (log10((int)mat[i * grid->order + j]) + 1 > 2) {
					printf("%4d ",(int)mat[i * grid->order + j]);
				}
			}
			printf("\n");                                                                                                                       
		}
		printf("\n");                                                                                                                                       
	}       
}


void init_grid(GRID* grid_info) {
    int rank;
    int dims[2];
    int period[2];
    int coords[2];
    int free_coords[2];

    MPI_Comm_size(MPI_COMM_WORLD, &(grid_info->nb_proc));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    grid_info->order = grid_info->nb_proc;
    dims[0] = dims[1] = (int) sqrt(grid_info->order);
    period[0] = period[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 1, &(grid_info->grid_comm));
    MPI_Comm_rank(grid_info->grid_comm, &(grid_info->grid_rank));
    MPI_Cart_coords(grid_info->grid_comm, grid_info->grid_rank, 2, coords);

    grid_info->pos_row = coords[0];
    grid_info->pos_col = coords[1];

    free_coords[0] = 0; 
    free_coords[1] = 1;
    MPI_Cart_sub(grid_info->grid_comm, free_coords, &(grid_info->row_comm));

    free_coords[0] = 1; 
    free_coords[1] = 0;
    MPI_Cart_sub(grid_info->grid_comm, free_coords, &(grid_info->col_comm));
} 


void read_matrix(GRID* info, char* namefile, float* mat, int order) {
	FILE* file = fopen(namefile, "r"); 

	if(file == NULL) {
		perror("invalid read");
		exit(1);  
	}

	int col; 
	int row;
	for(row = 0; row < sqrt(order); ++row)
	{
		for(col = 0; col < sqrt(order); ++col)
		{
			fscanf(file, "%f", &mat[row * order + col]);
		}
	}

	fclose(file);
}


int main(int argc, char** argv) {
	GRID grid;

	float* block_A; 
	float* block_B;
	float* block_C;
    
	float* mat_A;
    float* mat_B; 
    float* mat_C;
	
	double time1, time2, time3, time4, time5; 

	MPI_Init(&argc, &argv);
	init_grid(&grid);

	time1 = MPI_Wtime();

	if(grid.grid_rank == 0) {
		mat_A = (float*) malloc(grid.order*grid.order * sizeof(float));
		mat_B = (float*) malloc(grid.order*grid.order * sizeof(float));
		mat_C = (float*) malloc(grid.order*grid.order * sizeof(float));

		read_matrix(&grid, "A.txt", mat_A, grid.order);
		read_matrix(&grid, "B.txt", mat_B, grid.order);		
	
		int mat_col;
		int mat_row;

		for(mat_row = 0; mat_row < grid.order; ++mat_row) {
			for(mat_col = 0; mat_col < grid.order; ++mat_col) {
				mat_C[mat_row*grid.order+mat_col] = 0.0; 
			}
		}
	}

	if(argc > 0) {
        if(strcasecmp(argv[1], "linear") == 0) {
            if(grid.grid_rank==0) {
				printf("Linear computation\n \n");        
                multiply_matrix_linear(grid.nb_proc, mat_A, mat_B, mat_C);
                print_matrix(&grid, mat_C);
			}
			
			MPI_Finalize();
			return 0;  
		}
    }

	MPI_Datatype blocktype, type;	

	int array_size[2] = {grid.nb_proc, grid.nb_proc};
	int subarray_sizes[2] = {(int)sqrt(grid.nb_proc), (int) sqrt(grid.nb_proc)};
	int array_start[2] = {0,0};
	
	MPI_Type_create_subarray(2, array_size, subarray_sizes, array_start, MPI_ORDER_C, MPI_FLOAT, &blocktype); 
	MPI_Type_create_resized(blocktype, 0, (int)sqrt(grid.nb_proc)*sizeof(float), &type);
	MPI_Type_commit(&type);
	
	int i, j;
	int displs[grid.nb_proc];
	int send_counts[grid.nb_proc];

	block_A = (float*) malloc(grid.nb_proc*sizeof(float));
	block_B = (float*) malloc(grid.nb_proc*sizeof(float));
	block_C = (float*) malloc(grid.nb_proc*sizeof(float));

	for(i = 0; i < grid.nb_proc; ++i) {
		block_C[i] = 0.0;
	}

	if (grid.grid_rank == 0) {
		for(i=0; i < grid.nb_proc; i++) {
			send_counts[i] = 1;
		}

		int disp = 0;
		for (i = 0; i < (int)sqrt(grid.nb_proc); i++) {
			for (j = 0; j < (int)sqrt(grid.nb_proc); j++) {
				displs[i * (int)sqrt(grid.nb_proc) + j] = disp;
				disp += 1;
			}
			disp += ((grid.nb_proc / (int)sqrt(grid.nb_proc)-1)) * (int)sqrt(grid.nb_proc);
		}
	}

	time2 = MPI_Wtime();

	MPI_Scatterv(mat_A, send_counts, displs, type, block_A,grid.nb_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);	
	MPI_Scatterv(mat_B, send_counts, displs, type, block_B,grid.nb_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);                                              

	if(grid.grid_rank == 0) {
		printf("Cannon's algorithm\n \n");
	}

	time3 = MPI_Wtime();

	cannon(&grid, block_A, block_B, block_C);

	time4 = MPI_Wtime();

	MPI_Gatherv(block_C, grid.nb_proc,  MPI_FLOAT, mat_C, send_counts, displs, type, 0, MPI_COMM_WORLD);

	time5 = MPI_Wtime();

	print_matrix(&grid, mat_A);
	print_matrix(&grid, mat_B);
	print_matrix(&grid, mat_C);

	if(grid.grid_rank == 0) {
		printf("%f, %f, %f, %f, %f \n", time5-time4, time4-time3, time3-time2, time2-time1, time5-time1 );
	}

	MPI_Finalize();
	free(block_A);
    free(block_B);
    free(block_C);
    free(mat_A);
    free(mat_B);
    free(mat_C);
	
	return 0;
}
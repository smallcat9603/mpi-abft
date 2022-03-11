// modified from https://github.com/mperlet/matrix_multiplication
// Author: huyao 220310

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

typedef struct {
        unsigned int rows;
        unsigned int cols;
        double **mat_data;
} matrix_struct;

matrix_struct *get_matrix_struct(char matrix[]) {
    matrix_struct *m = malloc(sizeof(matrix_struct));
    m->rows = 0;
    m->cols = 0;
    FILE* myfile = fopen(matrix, "r");
    
    if(myfile == NULL) {
        printf("Error: The file you entered could not be found.\n");
        exit(EXIT_FAILURE);
    }
    // get the rows and columns
    int ch = 0;
    do {
        ch = fgetc(myfile);
        
        // count the columns at the first line (looking for "\t")
        if(m->rows == 0 && ch == '\t')
            m->cols++;
        
        // count the rows with "\n"
        if(ch == '\n')
            m->rows++;
            
    } while (ch != EOF);
    
    // write rows and cols to struct
    m->cols++;
    
    // allocate memory for matrix data
    m->mat_data = calloc(m->rows, sizeof(double*)); 
    int i;
    for(i=0; i < m->rows; ++i)
        m->mat_data[i]=calloc(m->cols, sizeof(double));
        
    
    rewind(myfile);
    int x,y;
    
    // fill matrix with data
    for(x = 0; x < m->rows; x++) {
        for(y = 0; y < m->cols; y++) {
            if (!fscanf(myfile, "%lf", &m->mat_data[x][y])) 
            break;
        }
    }
    
    fclose(myfile);

    return m;
}

void print_matrix(matrix_struct *matrix_to_print){
    int i,j;
    for(i = 0; i < matrix_to_print->rows; i++) {
        for(j = 0; j < matrix_to_print->cols; j++) {
            printf("%lf\t",matrix_to_print->mat_data[i][j]); //Use lf format specifier, \n is for new line
        }
        printf("\n");
    }
}

void free_matrix(matrix_struct *matrix_to_free) {
    for(int i = 0; i < matrix_to_free->rows; i++) {
        free(matrix_to_free->mat_data[i]);
    }
    free(matrix_to_free->mat_data);
    free(matrix_to_free);
}

double *mat_2D_to_1D(matrix_struct *m) {
    double *matrix = malloc( (m->rows * m->cols) * sizeof(double) );
    for (int i = 0; i < m->rows; i++) {
        memcpy( matrix + (i * m->cols), m->mat_data[i], m->cols * sizeof(double) );
    }
    return matrix;
}

uint32_t do_crc32(unsigned char *data_bits, int bytes)
{
    uint32_t crc = crc32(0L, Z_NULL, 0);

    for (int i = 0; i < bytes; i++)
    {
        crc = crc32(crc, data_bits + i, 1);
    }

    return crc;
}

//calculate block size (bytes)
int block_size(int data_bytes)
{
  double ber = 1e-5; //array length bug in hamming_encode() if 1e-7=0.0000001, hamming encode too slow if 1e-6 (because block size is too large)
  uint64_t bits = 1/ber; //10000000
  uint64_t bytes = bits/8; //1250000
  int bs = data_bytes; 
  if(bs > bytes)
  {
    bs = bytes;
  }
  return bs;
}

//1 bit --> char '0' or '1' (8 bits)
void cast_bits_to_char(unsigned char* bits, char* data, int bytes)
{
  for(int i = 0; i < bytes; i++)
  {
    for(int j = 0; j < 8; j++)
    {
      int pos = i * 8 + j;
      data[pos] = ((bits[i] >> (7-j)) & 1) + '0';
    }
  }
}

//calculate r
int hmLength(int k)
{
  int r = 0, flag = 1;
  while(flag)
  {
    int temp = pow(2, r);
    temp = temp - 1;
    flag = (temp-r-k<0);
    r++;
  }
  return r-1;
}

//c: hamming code, length(data) = k, length(c) = r+1 (secded)
void hamming_code(char* data, char* c, int k, int r)
{
  for(int i=0; i<r; i++)
  {
    int sum = 0, dnum = 0, cnum = 0;
    for(int j=1; j<r+k+1; j++)
    {
      if(j == (int)pow(2, cnum))
      {
        cnum++;
      }
      else
      {
        int x = pow(2, i);
        int y = j%(x*2);
        x = y/x;
        sum += (data[dnum]-'0')*x;
        dnum++;
      }
    }
    c[i] = sum%2 == 0?'0':'1';
  }

  //one additional parity bit for two bit error detection (secded)
  int sum = 0;
  for(int i=0; i<k; i++)
  {
    sum += data[i] - '0';
  }
  for(int i=0; i<r; i++)
  {
    sum += c[i] - '0';
  }
  c[r] = sum%2 + '0';
}

//hamming verify --> v
void hamming_verify_bit(unsigned char* bits, char* c, int bytes, int r, char* v)
{
  for(int i=0; i<r; i++)
  {
    int sum = 0, bnum = 0, cnum = 0;
    for(int j=1; j<r+bytes*8+1; j++)
    {
      if(j == (int)pow(2, cnum))
      {
        cnum++;
      }
      else
      {
        int x = pow(2, i);
        int y = j%(x*2);
        x = y/x;
        int bits_num = bnum/8;
        int bit_num = bnum%8;       
        sum += ((bits[bits_num] >> (7 - bit_num)) & 1) * x;
        bnum++;
      }
    }
    v[i] = sum%2 == (c[i]-'0')?'0':'1';
  }

  int sum = 0;
  for(int i=0; i<bytes; i++)
  {
    for(int j=0; j<8; j++)
    {
      sum += (bits[i] >> j) & 1;
    }
  }
  for(int i=0; i<r; i++)
  {
    sum += c[i] - '0';
  }
  v[r] = sum%2 == (c[r]-'0')?'0':'1';
}

//rectify error bit if one error (except parity bit)
void hamming_rectify_bit(unsigned char* bits, char* c, int bytes, int r, int error_bit_pos)
{
  int bnum = 0, cnum = 0;

  for(int j=1; j<r+bytes*8+1; j++)
  {
    if(j == (int)pow(2, cnum))
    {   
      if(j == error_bit_pos)
      {
        c[cnum] = c[cnum] == '0'?'1':'0';
        break;
      }
      else
      {
        cnum++;
      }  
    }
    else
    {  
      if(j == error_bit_pos)
      {
        int bits_num = bnum/8;
        int bit_num = bnum%8;
        bits[bits_num] ^= 1 << (7 - bit_num);
        break;
      }
      else
      {
        bnum++;
      } 
    }
  }
}

//calculate bit error position if possible one bit error
int error_info(char* v, int r, int* error_bit_pos)
{
  int error_type = 0;

  for(int i=0; i<r; i++)
  {
    *error_bit_pos += (v[i]-'0')*pow(2, i);
  }

  if(*error_bit_pos > 0 && v[r] == '0') //two bit error
  {
    error_type = 1;
  }
  else if(*error_bit_pos == 0 && v[r] == '1') //parity bit error
  {
    error_type = 2;
  }
  else if(*error_bit_pos > 0 && v[r] == '1') //one bit error
  {
    error_type = 3;
  }

  return error_type;
}

//data bits --> char[]
void hamming_encode(unsigned char* bits, char** c, int bytes, int* r)
{
  char data[bytes*8];
  cast_bits_to_char(bits, data, bytes);

  *r = hmLength(bytes*8);
  *c = (char*) malloc(sizeof(char)*(*r+1));
  hamming_code(data, *c, bytes*8, *r);
}

int hamming_decode(unsigned char* bits, char* c, int bytes, int r)
{
  char v[r+1];
  hamming_verify_bit(bits, c, bytes, r, v);

  int error_bit_pos = 0; // 1 for most left bit
  int error_type = error_info(v, r, &error_bit_pos);  
  switch(error_type)
  {
    case 0:
      // printf("no error\n");
      break;
    case 1:
      printf("two-bit error\n");
      break;
    case 2:
      printf("parity error\n");
      c[r] = c[r] == '0'?'1':'0';
      break;
    case 3:
      printf("one bit error: pos = %d\n", error_bit_pos);
      hamming_rectify_bit(bits, c, bytes, r, error_bit_pos);
      break;
    default:
      printf("ERROR");
  }

  return error_type;
}

//bit flip
void bit_flip(unsigned char* bits, int bytes)
{
  int num = rand() % (bytes*8);
  int byte = num/8;
  int bit = num%8;

  bits[byte] ^= 1 << (7 - bit);
}

//abft: crc --> hamming
//root = 0
void MPI_Bcast_abft(double *buffer, int count, int rank, int procs, int* resend)
{
  uint32_t crc = 0;
  uint32_t crc_check = 0;
  unsigned char crc_ok = 'y';
  unsigned char* crc_ok_recv = NULL;   

  int data_bytes = count*sizeof(double);
  unsigned char* data_bits = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes);

  if(rank == 0)
  {
      memcpy(data_bits, buffer, data_bytes);
      crc = do_crc32(data_bits, data_bytes);		
  }

  //hamming for blocks
  int bs = block_size(data_bytes);
  int bs_num = data_bytes/bs;
  int bs_last = data_bytes%bs;
  if (bs_last > 0) bs_num++;
  // printf("bs_num = %d, bs_last = %d \n", bs_num, bs_last);
  unsigned char* blocks[bs_num];
  char* c[bs_num];
  int r[bs_num];
  if(rank == 0)
  { 
    for(int i=0; i<bs_num; i++)
    {
      // c[i] = NULL;
      // r[i] = 0;
      int bytes = bs;
      if(bs_last > 0 && i == bs_num-1)
      {
        bytes = bs_last;
      } 

      // blocks[i] = (unsigned char*)malloc(bytes*sizeof(unsigned char));
      // memcpy(blocks[i], &data_bits[i*bs], bytes); 
      // printf("i = %d, bs = %d, bytes = %d\n", i, bs, bytes);
      hamming_encode(&data_bits[i*bs], &c[i], bytes, &r[i]);
    }
  }
  //printf("rank = %d, bs = %d, bs_last = %d, bs_num = %d, r[0] = %d\n", rank, bs, bs_last, bs_num, r[0]);

  MPI_Bcast(data_bits, data_bytes, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&crc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  //hamming
  MPI_Bcast(r, bs_num, MPI_INT, 0, MPI_COMM_WORLD);
  for(int i=0; i<bs_num; i++)
  {
    if(rank != 0)
    {
      c[i] = (char*) malloc(sizeof(char)*(r[i]+1));
    }
    MPI_Bcast(c[i], r[i]+1, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  // printf("r[0] = %d, c[0][1] = %c\n", r[0], c[0][1]);


  if(rank != 0)
  {
        // uint64_t to = 1e7;
        // int errors = data_bytes*8/to;
        // for (int n = 0; n < errors; n++)
        // {
        // bit_flip(data_bits, data_bytes);
        // }

      crc_check = do_crc32(data_bits, data_bytes);
      
      if (crc == crc_check)
      {
          // printf("CRC passed\n");
          crc_ok = 'y';
      }  
      else
      {
          // printf("CRC NOT passed\n");
          crc_ok = 'y';
          for(int i=0; i<bs_num; i++)
          {
            int bytes = bs;
            if(bs_last > 0 && i == bs_num-1)
            {
              bytes = bs_last;
            } 
            int error_type = hamming_decode(&data_bits[i*bs], c[i], bytes, r[i]);
            if(error_type == 1) //two-bit error
            {
              crc_ok = 'n';
              break;
            }
            // else if(error_type == 2 || error_type == 3)
            // {
            //   (*hamming_correct)++;
            // }
          }
      }            
  }
  else
  {
      crc_ok_recv = (unsigned char *)malloc(procs*1*sizeof(unsigned char));
  }

  MPI_Gather(&crc_ok, 1, MPI_UNSIGNED_CHAR, crc_ok_recv, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  
  if(rank == 0)
  {
      for(int i = 0; i < procs; i++)
      {
          if(i != 0 && crc_ok_recv[i] == 'n')
          {
              MPI_Send(data_bits, data_bytes, MPI_UNSIGNED_CHAR, i, i, MPI_COMM_WORLD);
              (*resend)++;
          }
      }
  }
  else if(crc_ok == 'n')
  {
      MPI_Recv(data_bits, data_bytes, MPI_UNSIGNED_CHAR, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if(rank != 0) memcpy(buffer, data_bits, data_bytes);

  free(data_bits);
}

int main(int argc, char *argv[]) {
    /** Matrix Properties
     * [0] = Rows of Matrix A
     * [1] = Cols of Matrix A
     * [2] = Rows of Matrix B
     * [3] = Cols of Matrix B
     **/
    int matrix_properties[4];
     
    double *m_a = NULL;
    double *m_b = NULL;
    double *final_matrix = NULL;
    
    int num_worker, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    /** the master initializes the data **/
    if (rank == 0) {
        
        if(argc != 3){
            printf("ERROR: Please specify only 2 files.\n");
            exit(EXIT_FAILURE);
        }
            
        matrix_struct *m_1 = get_matrix_struct(argv[1]);
        matrix_struct *m_2 = get_matrix_struct(argv[2]);

        if(m_1->cols != m_2->rows){
            printf("ERROR: The number of columns of matrix A must match the number of rows of matrix B.\n");
            exit(EXIT_FAILURE);
        }
        
        if (m_1->rows % num_worker != 0) {
            printf("ERROR: Matrix can not be calculated with this number of tasks.\n");
            exit(EXIT_FAILURE);
        }
        
        // fill the property-array for workers
        matrix_properties[0] = m_1->rows;
        matrix_properties[1] = m_1->cols;
        matrix_properties[2] = m_2->rows;
        matrix_properties[3] = m_2->cols;
        
        /* generate 1D matrices for workers 
         * m_a is the 1D Matrix of m_1 
         * m_a is the 1D Matrix of m_1 
        */
        m_a = mat_2D_to_1D(m_1);
        m_b = mat_2D_to_1D(m_2);

        free_matrix(m_1);
        free_matrix(m_2);
    }

    // send the matrix properties to the workers
    MPI_Bcast(&matrix_properties, 4, MPI_INT, 0, MPI_COMM_WORLD);

    // calculate the 1D-sizes of the matrices
    int size_a   = matrix_properties[0] * matrix_properties[1];
    int size_b   = matrix_properties[2] * matrix_properties[3];
    int size_res = matrix_properties[0] * matrix_properties[3];
    
    // allocate memory for 1D-matrices
    if(rank == 0) {
        final_matrix = malloc( size_res * sizeof(double) );
    } else {
        m_a = malloc( size_a * sizeof(double) );
        m_b = malloc( size_b * sizeof(double) );
    }
 
    int resent = 0;
    
    // send 1D matrices to workers
    // MPI_Bcast(m_a, size_a , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(m_b, size_b , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast_abft(m_a, size_a, rank, num_worker, &resent);
    MPI_Bcast_abft(m_b, size_b, rank, num_worker, &resent); 
    
    // calculate the start- and endrow for worker  
    int startrow = rank * ( matrix_properties[0] / num_worker);
    int endrow = ((rank + 1) * ( matrix_properties[0] / num_worker)) -1;
    
    /* calculate sub matrices */
    int number_of_rows = size_res / num_worker;
    double *result_matrix = calloc(number_of_rows, sizeof(double));

    int position = 0;

    for (int i = startrow; i <= endrow; i++) {
        for (int j = 0; j < matrix_properties[3]; j++) {
            for (int k = 0; k < matrix_properties[2]; k++) {
                result_matrix[position] +=
                    m_a[ (i * matrix_properties[1] + k) ] *
                    m_b[ (k * matrix_properties[3] + j) ];
            }
            position++;
        }
    }
    
    free(m_a);
    free(m_b);
    
    /* collect the results */
    MPI_Gather(result_matrix, number_of_rows, MPI_DOUBLE,
           final_matrix, number_of_rows,  MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /** The master presents the results on the console */
    if (rank == 0){
        FILE* fp = fopen("mat_result.txt", "w");

        int size = matrix_properties[0] * matrix_properties[3];
        int i = 0;
        while (i < size) {
            // printf("%lf\t", final_matrix[i]);
            fprintf(fp, "%lf\t", final_matrix[i]);
            i++;
        
            if (i % matrix_properties[3] == 0)
                // printf("\n");
                fprintf(fp, "\n");
        }

        fclose(fp); 

		// printf("--------------------------------------------------\n");
		// printf("FINAL RESULTS:\n");	  
    //     printf("resent = %d (percentage = %f)\n", resent, resent/(2.0*(num_worker-1)));    
        // printf("hamming_correct = %d (percentage = %f)\n", hamming_correct, 1.0*hamming_correct/(resent+hamming_correct));         
    }
    
    free(result_matrix);
    free(final_matrix);
    
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

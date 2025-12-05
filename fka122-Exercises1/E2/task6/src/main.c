#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tools.h"
//#include "attempt.h"


/*
 * The following functions are to be used:
*/
/*
double
average(
        double *v,
        unsigned int len
       );
double
variance(
         double *v,
         unsigned int len
        );
*/


/* ***************************************
*
* Autocorrelation
*
* Parameters
* ----------
*  data - the raw data from which the
*         autocorrelation should be
*         calculated
*  data_len     - the len of data
*  time_lag_ind - the lag (index) at 
*                 which the autocorrelation
*                 should be calculated
*
* Returns
* -------
*  The autocorrelation at a specific time
*  lag
*
* ***************************************/
double autocorrelation(
			   double *data,
			   int data_len,
			   int time_lag_ind)
{
    double ret = 0;
    double fi_average = average(data, data_len);
    double new_data[data_len];
    for(unsigned int i = 0; i < data_len - time_lag_ind; i++){
        new_data[i]= data[i] * data[i + time_lag_ind];
    }
    double new_data_average = average(new_data, data_len - time_lag_ind);
    double data_variance = standard_deviation(data, data_len);
    ret = (new_data_average - fi_average * fi_average) / (data_variance*data_variance);


    return ret;
}

/* ***************************************
*
* Block average
*
* Parameters
* ----------
*  data       - the raw data from which the
*               autocorrelation should be
                calculated
*  data_len   - the length of data
*  block_size - the size of the block
*
* Returns
* -------
* The statistical inefficency for a given
* block size
*
* ***************************************/
double block_average(double *data,
	                int data_len,
	                int block_size
	               )
{
    double ret = 0;
    int n_blocks = data_len / block_size;
    double block_averages[n_blocks];
    for(unsigned int i = 0; i < n_blocks; i++){
        double sum = 0.;
        for(unsigned int j = 0; j < block_size; j++){
            sum += data[i*block_size + j];
        }
        block_averages[i] = sum / block_size;
    }
    ret = standard_deviation(block_averages, n_blocks)*standard_deviation(block_averages , n_blocks) / (standard_deviation(data, data_len)*standard_deviation(data, data_len));

    return ret;
}


int save_2d_csv(const char *filename, double **data,
                unsigned int rows, unsigned int cols)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen");
        return -1;
    }

    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < cols; ++j) {
            fprintf(fp, "%.10e", data[i][j]);
            if (j + 1 < cols) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}


int main()
{

    return 0;
}
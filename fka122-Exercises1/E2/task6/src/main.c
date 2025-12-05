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
double variance(
         double *v1,
         unsigned int len
        )
{
    double avg = average(v1, len);
    double sum = 0.;
    for(unsigned int i = 0; i < len; i++){
        sum += pow(v1[i]-avg, 2);
    }
    return sum/len;
}


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
    double block_variance = standard_deviation(block_averages, n_blocks)*standard_deviation(block_averages , n_blocks);
    double data_variance = standard_deviation(data, data_len)*standard_deviation(data, data_len);
    ret = block_size*block_variance / data_variance;
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
    //read MC.txt
    FILE *fp = fopen("MC.txt", "r");
    if (!fp) {
        perror("fopen");
        return -1;
    }  
    int n_lines = 0;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), fp)) {
        n_lines++;
    }
    rewind(fp); 
    double *data = malloc(sizeof(double) * n_lines);
    for (int i = 0; i < n_lines; i++) {
        fscanf(fp, "%lf", &data[i]);
        //printf("the line %d: %f\n", i, data[i]);//test print

    }
    printf("the data: %f\n", data[0]);//test print
    fclose(fp);
    //calculate autocorrelation for different time lags
    double *autocorr_values = malloc(sizeof(double) * (n_lines/2));
    int auto_correlation_s=0;
    for (int lag = 0; lag < 2500; lag++) {
        autocorr_values[lag] = autocorrelation(data, n_lines, lag);
        printf("Lag %d: Autocorrelation = %.5f\n", lag, autocorr_values[lag]);
        if(fabs(autocorr_values[lag]-0.135) < 0.001){
            printf("The integrated autocorrelation time is approximately: %d\n", lag);
            auto_correlation_s = lag;
            break;
        }
    }
    double block_s = block_average(data, n_lines, 1000);

    printf("block_s: %f, auto_correlation_s: %d\n", block_s, auto_correlation_s);
    //save autocorrelation to csv
    double **autocorr_2d = create_2D_array(n_lines/2, 2);
    for (int i = 0; i < 15000; i++) {
        autocorr_2d[i][0] = i;
        autocorr_2d[i][1] = autocorr_values[i];
    }
    save_2d_csv("autocorrelation.csv", autocorr_2d, n_lines/2, 2);
    //free memory
    free(data);
    free(autocorr_values);
    destroy_2D_array(autocorr_2d, n_lines/2);

    return 0;
}
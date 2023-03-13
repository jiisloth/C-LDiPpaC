#include "hmatrix.h"


#define HARDMAXITER 16



void hardDecisionDecoding(int max_iter, int message[BLOCKSIZE], int output[BLOCKSIZE], int *iterations){

#pragma HLS INTERFACE mode=s_axilite port=iterations
#pragma HLS INTERFACE mode=s_axilite port=output
#pragma HLS INTERFACE mode=s_axilite port=message
#pragma HLS INTERFACE mode=s_axilite port=max_iter
#pragma HLS INTERFACE mode=s_axilite port=return


    int c;
    int i;
    int j;

    int maxi = max_iter;
    int itercount = 0;

    int v_nodes[BLOCKSIZE];
    int decode[BLOCKSIZE];

    for (j = 0; j < BLOCKSIZE; ++j) {
    	decode[j] = message[j];
        v_nodes[j] = decode[j];
    }



    for (int iter = 0; iter < HARDMAXITER; ++iter) {
    	if (iter < maxi){
    	    int satisfied = 0;
            for (c = 0; c < CNODES; ++c) {
                int c_node = 0;
                for (i = 0; i < PARITYCHECKS; ++i) {
                	for (j = 0; j < BLOCKSIZE; ++j) {
                        if (h_matrix[j][i] == c) {
                            c_node ^= decode[j];
                        }
                    }
                }
                satisfied |= c_node;
                for (i = 0; i < PARITYCHECKS; ++i) {
                	for (j = 0; j < BLOCKSIZE; ++j) {
                        if (h_matrix[j][i] == c) {
                            v_nodes[j] += decode[j]^c_node;
                        }
                    }
                }
            }

            if (satisfied == 0){
                itercount = iter +1;
                break;
            }

            for (j = 0; j < BLOCKSIZE; ++j) {
                if (v_nodes[j] > (PARITYCHECKS+1)/2){
                	decode[j] = 1;
                } else if (v_nodes[j] < (PARITYCHECKS+1)/2){
                	decode[j] = 0;
                }
                v_nodes[j] = decode[j];
            }
    	}
    }
    *iterations = itercount;
    for (j = 0; j < BLOCKSIZE; ++j) {
    	output[j] = decode[j];
    }
}

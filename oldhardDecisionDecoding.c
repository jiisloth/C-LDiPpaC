#include "hmatrix.h"


#define true 1
#define false 0

void hardDecisionDecoding(unsigned int max_iter, unsigned int message[MSGLEN], unsigned int output[MSGLEN], unsigned int *iterations){

	#pragma HLS INTERFACE s_axilite port=return bundle=CRTLS
	#pragma HLS INTERFACE s_axilite port=max_iter bundle=CRTLS
	#pragma HLS INTERFACE s_axilite port=message bundle=CRTLS
	#pragma HLS INTERFACE s_axilite port=output bundle=CRTLS
	#pragma HLS INTERFACE s_axilite port=iterations bundle=CRTLS

	//returns -1 if max_iter is reached
    *iterations = -1;

    unsigned int i = 0;
    unsigned int j = 0;

    for (unsigned int iter = 0; iter < max_iter; ++iter) {
        for (i = 0; i < MSGLEN; ++i) {
            output[i] = message[i];
        }
        short int satisfied = true;
        unsigned int c_nodes[HROWS] = {0};
        for (i = 0; i < HROWS; ++i) {
            for (j = 0; j < MSGLEN; ++j) {
                if (h_matrix[j][i] == 1){
                    c_nodes[i] ^= output[j];
                }
            }
            for (j = 0; j < MSGLEN; ++j) {
                if (h_matrix[i][j] == 1) {
                	message[j] += output[j] ^ c_nodes[i];
                }
            }
            if (c_nodes[i] == 1){
                satisfied = false;
            }
        }
        if (satisfied){
            *iterations = iter +1;
            break;
        }
        for (j = 0; j < MSGLEN; ++j) {
            unsigned int div = 1;
            for (i = 0; i < HROWS; ++i) {
                div += h_matrix[i][j];
            }
            float v = (float)message[j] /(float) div;
            if (v > 0.5){
            	message[j] = 1;
            } else {
            	message[j] = 0;
            }
        }
        for (i = 0; i < MSGLEN; ++i) {
        	output[i] = message[i];
        }
    }
}

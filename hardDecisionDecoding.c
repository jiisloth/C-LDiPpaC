#define MSGLEN 50
#define HROWS 20


#define true 1
#define false 0


void hardDecisionDecoding(int max_iter, int h_matrix[HROWS][MSGLEN], int message[MSGLEN], int output[MSGLEN], int *iterations){
    //returns -1 if max_iter is reached
    *iterations = -1;

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < MSGLEN; ++i) {
            output[i] = message[i];
        }
        int satisfied = true;
        int c_nodes[HROWS] = {0};
        for (int i = 0; i < HROWS; ++i) {
            for (int j = 0; j < MSGLEN; ++j) {
                if (h_matrix[i][j] == 1){
                    c_nodes[i] ^= message[j];
                }
            }
            for (int j = 0; j < MSGLEN; ++j) {
                if (h_matrix[i][j] == 1) {
                    output[j] += message[j] ^ c_nodes[i];
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
        for (int j = 0; j < MSGLEN; ++j) {
            int div = 1;
            for (int i = 0; i < HROWS; ++i) {
                div += h_matrix[i][j];
            }
            float v = (float)output[j] /(float) div;
            if (v > 0.5){
                output[j] = 1;
            } else {
                output[j] = 0;
            }
        }
        for (int i = 0; i < MSGLEN; ++i) {
            message[i] = output[i];
        }
    }
}

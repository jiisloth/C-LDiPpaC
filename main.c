#include <stdio.h>
#include <math.h>

#define SIGMA 1
#define MSGLEN 15
#define HROWS 6
#define MAXITER 5

#define H2ROWS 4
#define MSGLEN2 8

#define true 1
#define false 0

double message[MSGLEN] = {-1.3436883, 0.92061442, 1.00983062, -0.06886619, 1.71192052, -1.15801087, -1.06911068, 1.96699053, 1.35679392, 1.79001563, -1.09911672, 1.45552054, -0.87099213, 0.10570174, -0.58468938};
int H[HROWS][MSGLEN] = {{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0}, {0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1}, {0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1}, {0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0}, {1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0}};

int H2[H2ROWS][MSGLEN2] = {{0, 1, 0, 1, 1, 0, 0, 1},{1, 1, 1, 0, 0, 1, 0, 0},{0, 0, 1, 0, 0, 1, 1, 1},{1, 0, 0, 1, 1, 0, 1, 0}};
int message2[MSGLEN2] = {1,1,0,1,0,1,0,1};
//{{1, 1, 0, 0, 0, 0, 0},{0, 1, 1, 1, 0, 0, 0},{0, 0, 1, 1, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0}, {0, 0, 0, 0, 1, 1, 0}, {0, 0, 0, 1, 0, 0, 1}};


void calculateBits(double p0[], double p1[], int *bits){
    for (int i = 0; i < MSGLEN; ++i) {
        if (p0[i] > p1[i]){
            bits[i] = 0;
        } else {
            bits[i] = 1;
        }
    }
}

void print_bits(int* bits){
    for(int loop = 0; loop < MSGLEN; loop++)
        printf("%d", bits[loop]);
    printf("\n");
}

void sumProductDecoding() {
    double p0[MSGLEN];
    double p1[MSGLEN];
    double P0[HROWS][MSGLEN];
    double P1[HROWS][MSGLEN];
    int bits[MSGLEN];
    // init arrays
    for (int i = 0; i < MSGLEN; ++i) {
        p1[i] = 1 / (1 + exp((2 * message[i]) / SIGMA));
        p0[i] = 1 - p1[i];
    }
    // get original message
    calculateBits(p0, p1, bits);
    print_bits(bits);
    // init P matrix's
    for (int i = 0; i < HROWS; ++i) {
        for (int j = 0; j < MSGLEN; ++j) {
            P0[i][j] = H[i][j] * p0[j];
            P1[i][j] = H[i][j] * p1[j];
        }
    }
    // actual algorithm
    for (int iter = 0; iter < MAXITER; ++iter) {
        // DeltaP
        double deltaP[HROWS][MSGLEN];
        for (int i = 0; i < HROWS; ++i) {
            for (int j = 0; j < MSGLEN; ++j) {
                deltaP[i][j] = P0[i][j] - P1[i][j];
            }
        }
        // Q:s
        double deltaQ[HROWS][MSGLEN];
        double Q0[HROWS][MSGLEN];
        double Q1[HROWS][MSGLEN];
        for (int i = 0; i < HROWS; ++i) {
            for (int j = 0; j < MSGLEN; ++j) {
                double product = 1;
                if (deltaP[i][j] != 0){
                    for (int j2 = 0; j2 < MSGLEN; ++j2) {
                        if (j != j2 && deltaP[i][j2] != 0){
                            product = product * deltaP[i][j2];
                        }
                        deltaQ[i][j] = product;
                        Q0[i][j] = (0.5 * (1+product));
                        Q1[i][j] = (0.5 * (1-product));
                    }
                } else {
                    deltaQ[i][j] = 0;
                    Q0[i][j] = 0;
                    Q1[i][j] = 0;
                }
            }
        }
        //new P's
        double np0[MSGLEN] = {0};
        double np1[MSGLEN] = {0};
        for (int i = 0; i < HROWS; ++i) {
            for (int j = 0; j < MSGLEN; ++j) {
                double product = 1;
                double product2 = 1;
                if (Q0[i][j] != 0) {
                    for (int i2 = 0; i2 < HROWS; ++i2) {
                        if (Q0[i2][j] != 0) {
                            product2 = product2 * Q0[i2][j];
                            if (i != i2) {
                                product = product * Q0[i2][j];
                            }
                        }
                    }
                    np0[j] = p0[j] * product2;
                    P0[i][j] = p0[j] * product;
                } else {
                    P0[i][j] = 0;
                }
                product = 1;
                product2 = 1;
                if (Q1[i][j] != 0) {
                    for (int i2 = 0; i2 < HROWS; ++i2) {
                        if (Q1[i2][j] != 0) {
                            product2 = product2 * Q1[i2][j];
                            if (i != i2) {
                                product = product * Q1[i2][j];
                            }
                        }
                    }
                    np1[j] = p1[j] * product2;
                    P1[i][j] = p1[j] * product;
                } else {
                    P1[i][j] = 0;
                }
            }
        }
        for (int i = 0; i < MSGLEN; ++i) {
            p0[i] = np0[i];
            p1[i] = np1[i];
        }
        //scale
        for (int j = 0; j < MSGLEN; ++j) {
            for (int i = 0; i < HROWS; ++i) {
                if (P0[i][j] + P1[i][j] != 0) {
                    double factor = 1 / (P0[i][j] + P1[i][j]);
                    P0[i][j] = P0[i][j] * factor;
                    P1[i][j] = P1[i][j] * factor;
                }
            }
            if (p0[j] + p1[j] != 0) {
                double factor = 1 / (p0[j] + p1[j]);
                p0[j] = p0[j]*factor;
                p1[j] = p1[j]*factor;
            }

        }
        calculateBits(p0, p1, bits);
        print_bits(bits);


    }
}


void hardDecisionDecoding(){
    for (int iter = 0; iter < MAXITER; ++iter) {
        int out[MSGLEN2];
        for (int i = 0; i < MSGLEN2; ++i) {
            out[i] = message2[i];
        }
        int satisfied = true;
        int c_nodes[H2ROWS] = {0};
        for (int i = 0; i < H2ROWS; ++i) {
            for (int j = 0; j < MSGLEN2; ++j) {
                if (H2[i][j] == 1){
                    // TODO: fix this with or statement / use single bits instead of int..
                    c_nodes[i] ^= message2[j];
                    if (c_nodes[i] == 2){
                        c_nodes[i] = 0;
                    }
                }
            }
            for (int j = 0; j < MSGLEN2; ++j) {
                if (H2[i][j] == 1) {
                    int test = message2[j] + c_nodes[i];
                    if (test == 2){
                        test = 0;
                    }
                    out[j] += test;

                }
            }
            if (c_nodes[i] == 1){
                satisfied = false;
            }
        }
        if (satisfied){
            break;
        }
        for (int i = 0; i < H2ROWS; ++i) {
            printf("%d ", c_nodes[i]);
        }
        printf("\n");
        for (int j = 0; j < MSGLEN2; ++j) {
            int div = 1;
            for (int i = 0; i < H2ROWS; ++i) {
                div += H2[i][j];
            }
            double v = (double)out[j] /(double) div;
            if (v > 0.5){
                out[j] = 1;
            } else {
                out[j] = 0;
            }
        }
        for (int i = 0; i < MSGLEN2; ++i) {
            message2[i] = out[i];
        }
        for (int i = 0; i < MSGLEN2; ++i) {
            printf("%d ", message2[i]);
        }
        printf("\n");
    }
    printf("Final result:\n");
    for (int i = 0; i < MSGLEN2; ++i) {
        printf("%d ", message2[i]);
    }
    printf("\n");
}


int main() {
    printf("Sum product:\n");
    sumProductDecoding();
    //printf("Hard decision:\n");
    //hardDecisionDecoding();
    return 0;
}
#include <stdio.h>
#include "testcases.h"
#include "Testbench.h"

#define MAXITER 5

int main(){
    printf("Running Hard Decision Decoding..\n");
    int errors = 0;
    for (int t = 0; t < TESTCASECOUNT; t++){

        printf("\nTestcase: %d", t);

        errors += run_case(*cases[t]);
    }

    printf("\nFailed decoding: %d/%d\n", errors, TESTCASECOUNT);

    return 0;
}


int run_case(Testcase custard){
    // custard = case but cant use case as variable...

    //get binary message and initial errors..
    int errorsinmessage = 0;
    int message[MSGLEN];
    int i;

    for(i = 0; i < MSGLEN; i++) {
        if (custard.message[i] < 0){
            message[i] = 1;
        } else {
            message[i] = 0;
        }
        if (message[i] != custard.reference[i]){
            errorsinmessage = 1;
        }
    }
    printf("Error bits: %d.\n", errorsinmessage);

    int output[MSGLEN];
    int iterations = -1;


    hardDecisionDecoding(MAXITER, custard.h_matrix, message, output, &iterations);

    if (iterations < 0){
        printf("- Max iterations reached.\n");
    } else {
        printf("- Finished on iteration # %d.\n", iterations);
    }

    int errors = 0;
    int errorline[MSGLEN] = {0};

    printf("- Output:\n");

    for(i = 0; i < MSGLEN; i++){
        if(output[i] != custard.reference[i]){
            errors += 1;
            errorline[i] = 1;
        }
        printf("%d", output[i]);
    }
    if (errors > 0) {
        for (i = 0; i < MSGLEN; i++) {
            if (errorline[i] == 1) {
                printf("*");
            } else {
                printf(" ");
            }
        }
        printf("- Decoding failed with %d errors!\n", errors);
        return 1;
    } else {
        printf("- No errors found!\n");
        return 0;
    }
}
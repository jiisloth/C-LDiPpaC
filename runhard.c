#include <stdio.h>
#include "platform.h"
#include "xil_printf.h"
#include <xparameters.h>
#include "testcases.h"
#include <xsumproductdecoding.h>
#include <xharddecisiondecoding.h>
#include "xtime_l.h"

#define MAXITER 5
#define ENDSAME 2

#define SUMPRODUCT 0
#define HARDDECISION 1

int run_case(Testcase custard, int decoder, int *iterations);
void run_sumProductDecoding(int max_iter, int end, int h_matrix[HROWS][MSGLEN], float message[MSGLEN], int output[MSGLEN], int *iterations);
void run_hardDecisionDecoding(int max_iter, int h_matrix[HROWS][MSGLEN], int message[MSGLEN], int output[MSGLEN], int *iterations);




int main(){

    int status;
    int spderrors = 0;
    int spditerations = 0;
    int hdderrors = 0;
    int hdditerations = 0;

    float spdtime;
    float hddtime;

    XTime gbl_time_before_test;
    XTime gbl_time_after_test;

    init_platform();

    XHarddecisiondecoding doHarddecisiondecoding;
    XHarddecisiondecoding_Config *doHarddecisiondecoding_cfg;

    XSumproductdecoding doSumproductdecoding;
    XSumproductdecoding_Config *doSumproductdecoding_cfg;

    doSumproductdecoding_cfg = XSumproductdecoding_LookupConfig(XPAR_SUMPRODUCTDECODING_0_DEVICE_ID);
    if (!doSumproductdecoding_cfg) {
        printf("Error loading conf for spd\n");
    }
    status = XSumproductdecoding_CfgInitialize(&doSumproductdecoding, doSumproductdecoding_cfg);
    if (status != XST_SUCCESS){
        printf("Error initializing spd\n");
    }

    doHarddecisiondecoding_cfg = XHarddecisiondecoding_LookupConfig(XPAR_HARDDECISIONDECODING_0_DEVICE_ID);
    if (!doHarddecisiondecoding_cfg) {
        printf("Error loading conf for hdd\n");
    }
    status = XHarddecisiondecoding_CfgInitialize(&doHarddecisiondecoding, doHarddecisiondecoding_cfg);
    if (status != XST_SUCCESS){
        printf("Error initializing hdd\n");
    }


    printf("\n######################################\n\n");


    printf("Running Sum Product Decoding..\n");


    XTime_GetTime(&gbl_time_before_test);

    for (int t = 0; t < TESTCASECOUNT; t++){

        printf("\nTestcase: %d, ", t+1);

        spderrors += run_case(*cases[t], SUMPRODUCT, &spditerations);
    }

    XTime_GetTime(&gbl_time_after_test);

    spdtime = (float) (gbl_time_after_test - gbl_time_before_test)/(COUNTS_PER_SECOND);

    printf("\nDecoded %d/%d messages correctly!\n", TESTCASECOUNT-spderrors, TESTCASECOUNT);
    printf("Took in average %.2f iterations\n", (float)spditerations/(float)TESTCASECOUNT);

    printf("\n######################################\n\n");

    printf("Running Hard Decision Decoding..\n");

    XTime_GetTime(&gbl_time_before_test);

    for (int t = 0; t < TESTCASECOUNT; t++){

        printf("\nTestcase: %d, ", t+1);

        hdderrors += run_case(*cases[t], HARDDECISION, &hdditerations);
    }

    XTime_GetTime(&gbl_time_after_test);

    printf("\nDecoded %d/%d messages correctly!\n", TESTCASECOUNT-hdderrors, TESTCASECOUNT);
    printf("Took in average %.2f iterations\n", (float)hdditerations/(float)TESTCASECOUNT);

    printf("\n######################################\n\n");

    XTime_GetTime(&gbl_time_after_test);

    hddtime = (float) (gbl_time_after_test - gbl_time_before_test)/(COUNTS_PER_SECOND);

    printf("SPD: %.1f%% correct with avg: %.2f iterations in %.3f seconds.\n", (float)(TESTCASECOUNT-spderrors)/(float)TESTCASECOUNT*100.0, (float)spditerations/(float)TESTCASECOUNT, spdtime);
    printf("HDD: %.1f%% correct with avg: %.2f iterations in %.3f seconds.\n", (float)(TESTCASECOUNT-hdderrors)/(float)TESTCASECOUNT*100.0, (float)hdditerations/(float)TESTCASECOUNT, hddtime);

    cleanup_platform();
    return 0;
}

int run_case(Testcase custard, int decoder, int *totaliter){
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

    if (decoder == SUMPRODUCT){
        run_sumProductDecoding(MAXITER, ENDSAME, custard.h_matrix, custard.message, output, &iterations);
    }
    else if (decoder == HARDDECISION){
        run_hardDecisionDecoding(MAXITER, custard.h_matrix, message, output, &iterations);
    }
    else {
        printf("\n WRONG DECODER PARAMETER.\n");
        return 1;

    }

    if (iterations < 0){
        printf("- Max iterations reached.\n");
        *totaliter += MAXITER;
    } else {
        printf("- Finished on iteration # %d.\n", iterations);
        *totaliter += iterations;
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
        printf("\n");
        for (i = 0; i < MSGLEN; i++) {
            if (errorline[i] == 1) {
                printf("*");
            } else {
                printf(" ");
            }
        }
        printf("\n- Decoding failed with %d errors!\n", errors);
        return 1;
    } else {
        printf("\n- Decoding successful!\n");
        return 0;
    }
}


void run_sumProductDecoding(int max_iter, int end, int h_matrix[HROWS][MSGLEN], float message[MSGLEN], int output[MSGLEN], int *iterations) {

    XSumproductdecoding_Set_max_iter(&doSumproductdecoding, max_iter);
    XSumproductdecoding_Set_end_r(&doSumproductdecoding, end);
    XSumproductdecoding_Write_h_matrix_Words(&doSumproductdecoding, 0, h_matrix, HROWS*MSGLEN);
    XSumproductdecoding_Write_message_Words(&doSumproductdecoding, 0, message, MSGLEN);
    XSumproductdecoding_Write_output_r_Words(&doSumproductdecoding, 0, output, MSGLEN);


    XSumproductdecoding_Start(&doSumproductdecoding);
    while(!XSumproductdecoding_IsDone(&doSumproductdecoding));


    XSumproductdecoding_Read_output_r_Words(&doSumproductdecoding, 0, output, MSGLEN);
    *iterations = XSumproductdecoding_Get_iterations(&doSumproductdecoding);
}

void run_hardDecisionDecoding(int max_iter, int h_matrix[HROWS][MSGLEN], int message[MSGLEN], int output[MSGLEN], int *iterations){
    XHarddecisiondecoding_Set_max_iter(&doHarddecisiondecoding, max_iter);
    XHarddecisiondecoding_Write_h_matrix_Words(&doHarddecisiondecoding, 0, h_matrix, HROWS*MSGLEN);
    XHarddecisiondecoding_Write_message_Words(&doHarddecisiondecoding, 0, message, MSGLEN);
    XHarddecisiondecoding_Write_output_r_Words(&doHarddecisiondecoding, 0, output, MSGLEN);

    XHarddecisiondecoding_Start(&doHarddecisiondecoding);
    while(!XHarddecisiondecoding_IsDone(&doHarddecisiondecoding));

    XHarddecisiondecoding_Read_output_r_Words(&doHarddecisiondecoding, 0, output, MSGLEN);
    *iterations = XHarddecisiondecoding_Get_iterations(&doHarddecisiondecoding);
}


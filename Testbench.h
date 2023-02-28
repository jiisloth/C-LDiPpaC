void sumProductDecoding(int max_iter, int end, int h_matrix[HROWS][MSGLEN], float message[MSGLEN], int output[MSGLEN], int *iterations);
void hardDecisionDecoding(int max_iter, int h_matrix[HROWS][MSGLEN], int message[MSGLEN], int output[MSGLEN], int *iterations);
int run_case(Testcase custard);
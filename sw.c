#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_REF "AATTGGGGCCCC"
#define TEST_QUERY "AATTGGGGCCCC"
#define REF_LEN strlen(TEST_REF)
#define QUERY_LEN strlen(TEST_QUERY)
#define MATCH 2.0
#define MISMATCH (-2.0)
#define GAP_OPEN (-3.0)
#define GAP_EXTENSION_DECAY (0.3)


void print_matrix(float matrix[REF_LEN+1][QUERY_LEN+1]);
void align(void);
float max(float list[], int len);
float sw_matrix[REF_LEN +1][QUERY_LEN +1] = {{0}};
/* main: ***********************************/
//static float decision_matrix[REF_LEN +1][QUERY_LEN +1];
int main(){
    align();
    print_matrix(sw_matrix);
    return 0;
}

/*******************************************/



void score(int, int);
/* align: fill in score for all sw_matrix entries */
void align(){
    int row, col;
    for (row = 1; row < REF_LEN+1; row++)
        for (col = 1; col < QUERY_LEN+1; col++)
            score(row, col);
}


float score_mm(int, int);
float score_insert(int, int);
float score_delete(int, int);
float max(float list[], int len);
/* score: compute the score for the [row][col] entry for sw_matrix */
void score(int row, int col){
    float best_score;
    float choices[4] = {0.0,
                        score_mm(row, col), 
                        score_insert(row, col), 
                        score_delete(row, col)};
    best_score = max(choices, 4);
    sw_matrix[row][col] = best_score;
    print_matrix(sw_matrix);
    printf("\n");
}


float penalise_gap(int gap_len);
float score_mm(int row, int col){
    float prefix_score = sw_matrix[row -1][col-1];
    if (TEST_REF[row -1] == TEST_QUERY[col -1])
        return prefix_score + MATCH;
    else
        return prefix_score + MISMATCH;
}


float score_insert(int row, int col){
    int gap_len;
    float score = 0;
    float new_score;
    for (gap_len = col; gap_len > 0; gap_len--){
        new_score= sw_matrix[row][col-gap_len] + penalise_gap(gap_len);
        if (new_score > score)
            score = new_score;
    }
    return score;
}

float score_delete(int row, int col){
    int gap_len;
    float score = 0; // min score on sw_matrix is 0
    float new_score;
    for (gap_len = row; gap_len > 0; gap_len--){
        new_score = sw_matrix[row - gap_len][col] + penalise_gap(gap_len); 
        if (new_score > score)
           score =  new_score;
    }
    return score;
}

float penalise_gap(int gap_len){
    float penalty;
    penalty = GAP_OPEN + GAP_EXTENSION_DECAY * gap_len;
    return (penalty < 0)? penalty : 0.0;
} 
    


//*********************************************//
/* max: return the largest float in list[] of
 * length len */
float max(float list[], int len){
    float maximum = list[0]; // assume len > 0
    register int i;
    for (i = 1; i < len; i++)
        if (list[i] > maximum)
            maximum = list[i];
    return maximum;
}

float min(float list[], int len){
    float minimum = list[0];
    int i;
    for (i = 1; i < len; i++)
        if (list[i] < minimum)
            minimum = list[i];
    return minimum;
}

/* print_matrix: print all entries of a given matrix of floats */
void print_matrix(float matrix[REF_LEN+1][QUERY_LEN+1]){
    int row, col;
    for (row = 0; row < REF_LEN +1; row++){
        for (col = 0; col < QUERY_LEN +1; col++)
            printf("%.2f ", matrix[row][col]);
        printf("\n");
    }
}


#include <string.h>
#include <stdio.h>

#define TEST_REF ("AATTGGGGACAGGGGCTATATATCGATCGATGGCTAGCGGCGGCGCGCGGGGG"\
                  "GGGCTCATACAGTGACGTACGTAGCATGACTGCATGTACGTAGTGCTGCGGGG"\
                  "GCGCGATCGATGCATGTCACCCACGTGACGATGCATGGCGGCGTATCGATCGA"\
                  "CGATCGACGATCGATCGTACGTGACTGAC")
#define TEST_QUERY TEST_REF
#define REF_LEN strlen(TEST_REF)
#define QUERY_LEN strlen(TEST_QUERY)

#define MATCH 2.0
#define MISMATCH (-1.0)
#define GAP_OPEN (-1.0)
#define GAP_EXTENSION_DECAY (0.0)

#define MAX_SEQ_LEN 60

void print_matrix(float matrix[MAX_SEQ_LEN+1][MAX_SEQ_LEN+1], int nrow, int ncol);
void align(char *ref, char *query);
static float sw_matrix[MAX_SEQ_LEN +1][MAX_SEQ_LEN +1] = {{0}};


int main(int argc, char **argv){
    char *ref, *query;
    ref = argv[1];
    query = argv[2];
    printf("%s\n%s\n", ref, query);
    align(ref, query);

    int nrow = strlen(ref);
    int ncol = strlen(query);
    print_matrix(sw_matrix, nrow, ncol);
    printf("length of ref sequence = %lu\n", strlen(ref));
    printf("length of query sequence = %lu\n", strlen(query));
    return 0;
}

/*****************************************************************************/



void score(int, int, char *, char *);
/* align: fill in score for all sw_matrix entries */
void align(char *ref, char *query){
    int row, col;
    int ref_len = strlen(ref);
    int query_len = strlen(query);

    for (row = 1; row < ref_len +1; row++)
        for (col = 1; col < query_len +1; col++)
            score(row, col, ref, query);
}


float score_mm(int, int, char *, char *);
float score_insert(int, int);
float score_delete(int, int);
float max(float list[], int len);
/* score: compute the score for the [row][col] entry for sw_matrix */
void score(int row, int col, char *ref, char *query){
    float best_score;
    float choices[4] = {0.0,
                        score_mm(row, col, ref, query), 
                        score_insert(row, col), 
                        score_delete(row, col)};
    best_score = max(choices, 4);
    sw_matrix[row][col] = best_score;
}


float penalise_gap(int gap_len);
float score_mm(int row, int col, char *ref, char *query){
    float prefix_score = sw_matrix[row -1][col-1];
    if (ref[row -1] == query[col -1])
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
void print_matrix(float matrix[MAX_SEQ_LEN+1][MAX_SEQ_LEN+1], int nrow, int ncol){
    int row, col;
    for (row = 0; row < nrow +1; row++){
        for (col = 0; col < ncol +1; col++)
            printf("%.2f ", matrix[row][col]);
        printf("\n");
    }
}


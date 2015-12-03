/* This module will be the main module in which
 * the SWinC aligner (Smith-Waterman alignmer) will
 * be implemented. 
 * In addition to the vanilla Smith-Waterman 
 * (bottom-up Dynamic Programming with
 * gap-scoring scheme), the algorithm will be modified
 * to suit its use in the context of DNA nucleotide 
 * sequence alignment. In particular, thermodynamic
 * consideration would be incorporated into the
 * scoring scheme. 
 *
 * The algorithm is described in the following paper:
 * Kaderali, L, "Primer Design for Multiplexed Genotyping",
 * Methods In Molecular Biology, Humana Press, Totowa,
 * New Jersy, pg. 269-285
 */

#include <stdio.h>
#include <string.h>

/* define maximum primer + heel length acceptable */
#define MAX_SEQ_LEN 60

/* produce record of nearest neighbour thermodynamic
 * data here. It should be in "hash table" like 
 * structure
 */


/* Record score and decision at every sw_matrix entry
 * the decision are 
 *   'M' for match
 *   'm' for mismatch
 *   'D' for deletion
 *   'I' for insertion
 *   '\0' for termimate alignment */
struct sw_entry {
    float score;
    char decision;
};

float swalign(char *ref, char *query);
void fill_matrix(int ref_len, int query_len);
void score(int row, int col);
struct sw_entry score_mm(int row, int col);
struct sw_entry score_insert(int row, int col);
struct sw_entry score_delete(int row, int col);
float penalise_gap(int gap_len);
struct sw_entry max(struct sw_entry list[], int list_len);
void print_matrix(int nrow, int ncol);

float match_score = 2.0;
float mismatch_penalty = -1.0;
float gap_open = -1.0;
float gap_extension = 0.0;

static struct sw_entry null_entry = {0.0, '\0'};
/* initialise sw_matrix as (MAX_SEQ_LEN+1) x (MAX_SEQ_LEN+1)
 * all entries will be initialised with struct sw_entry = {0.0,'\0'} */
static struct sw_entry sw_matrix[MAX_SEQ_LEN +1][MAX_SEQ_LEN +1];

/* external variable recording the reference and query sequence
 * visible to all the alignment/scoring related functions */
char *g_ref;
char *g_query;


/******************************************************************************/

/* main() */

int main(int argc, char **argv){
    char *ref = argv[1];
    char *query = argv[2];
    int ref_len = strlen(ref);
    int query_len = strlen(query);
    printf("ref=%s\nquery=%s\nref_len=%i query_len=%i\n", ref, query, ref_len, query_len);
    swalign(ref, query);
    print_matrix(query_len +1, ref_len +1);
    return 0;
}




/* swalign: produce the best alignment score for the 2 input string */
float swalign(char *ref, char *query){
    g_ref = ref;
    g_query = query;
    int ref_len = strlen(ref);
    int query_len = strlen(query);

    fill_matrix(ref_len, query_len);
    return 0.0;
}


void fill_matrix(int ref_len, int query_len){
    int row, col;
    for (row = 1; row < query_len +1; row++)
        for (col = 1; col < ref_len +1; col++)
            score(row, col);
}


/* score:
 * for the specifed row and col in sw_matrix
 * assign the best  */
void score(int row, int col){
    struct sw_entry choices[4] = 
                       {null_entry,
                        score_mm(row, col),
                        score_insert(row, col),
                        score_delete(row, col)};
    sw_matrix[row][col] = max(choices, 4);
}

struct sw_entry score_mm(int row, int col){
    float prefix_score = sw_matrix[row-1][col-1].score;
    struct sw_entry mm_entry;
    if (g_ref[col -1] == g_query[row -1]){
        mm_entry.score = prefix_score + match_score;
        mm_entry.decision = 'M';
    } else{
        mm_entry.score = prefix_score + mismatch_penalty;
        mm_entry.decision = 'm';
    };
    return mm_entry;
}

struct sw_entry score_insert(int row, int col){
    int gap_len;
    struct sw_entry insert_entry = {0.0, 'I'};
    float new_score;
    for (gap_len = col; gap_len > 0; gap_len--){
        new_score = sw_matrix[row][col-gap_len].score + penalise_gap(gap_len);
        if (new_score > insert_entry.score)
            insert_entry.score = new_score;
    };
    return insert_entry;
}

struct sw_entry score_delete(int row, int col){
    int gap_len;
    struct sw_entry delete_entry = {0.0, 'D'};
    float new_score;
    for (gap_len = row; gap_len > 0; gap_len--){
        new_score = sw_matrix[row-gap_len][col].score + penalise_gap(gap_len);
        if (new_score > delete_entry.score)
            delete_entry.score = new_score;
    };
    return delete_entry;
}

float penalise_gap(int gap_len){
    float penalty = gap_open + gap_extension * gap_len;
    return penalty;
}




/***** Utilities *********/

/* max:
 * return the entry in list[] with maximal score */
struct sw_entry max(struct sw_entry list[], int list_len){
    struct sw_entry maximum = list[0]; // assume list_len > 0
    register int i;
    for (i = 1; i < list_len; i++)
        if (list[i].score > maximum.score)
            maximum = list[i];
    return maximum;
}


void print_matrix(int nrow, int ncol){
    int row, col;
    struct sw_entry entry;
    for (row = 0; row < nrow; row++){
        for (col = 0; col < ncol; col++){
            entry = sw_matrix[row][col];
            printf("%.2f%c ", entry.score, entry.decision);
        };
        putchar('\n');
    }
}

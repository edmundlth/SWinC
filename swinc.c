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
#include <ctype.h>
#include "swinc.h"

/* main() */
int main(int argc, char **argv){
    //char *ref = argv[1];
    //char *query = argv[2];
    //int ref_len = strlen(ref);
    //int query_len = strlen(query);
    //float align_score;
    //printf("ref=%s\nquery=%s\nref_len=%i query_len=%i\n", ref, query, ref_len, query_len);
    //align_score = swalign(ref, query);
    //print_sw_matrix(query_len +1, ref_len +1);
    //printf("Alignment score = %.2f\n\n", align_score);

    char *filename = argv[3];
    int num_primer;
    float max_interaction[num_primer], mean_interaction[num_primer];
    num_primer = get_primers(filename);
    printf("num_primer = %i\n", num_primer);

    align_pool(num_primer);
    print_interaction_matrix(num_primer, num_primer);
    return 0;
}

/* align_pool:
 * given an array of strings, pool, return a matrix of size
 * |pool|x|pool| with matrix[i][j] representing the alignment score
 * between string_i and string_j (a complete graph).
 * The matrix will be symmetric and the diagonal element will be 0.0, i.e.
 * by default we don't want information about self-alignment.*/
void align_pool(int pool_size){
    int i, j;
    for (i = 0; i < pool_size; i++){
        for (j = 0; j < pool_size; j++){
            char *this_query;
            this_query = rev_complement(pool[j]);
            interaction_matrix[i][j] = (i==j)? 
                                       0.0 : 
                                       swalign(pool[i], this_query);
        }
    }
}


int get_primers(char *filename){
    FILE *file_handle = fopen(filename, "r");
    int i = 0;
    char buff[MAX_SEQ_LEN];
    char *temp;
    while ((temp = fgets(buff, MAX_SEQ_LEN, file_handle)) != NULL){
        temp = trim_whitespace(temp); //temp[strlen(temp)-1] = '\0';//remove trailing newline
        strcpy(pool[i], temp);
        i++;
    }
    fclose(file_handle);
    return i;
}

/******************** SW alignment algorithm ****************************/

/* swalign: produce the best alignment score for the 2 input string */
float swalign(char *ref, char *query){
    g_ref = ref;
    g_query = query;
    float best_score;
    int ref_len = strlen(ref);
    int query_len = strlen(query);

    best_score = fill_matrix(ref_len, query_len);
    return best_score;
}

/* fill_matrix:
 * fill the matrix with the the best scores and decisions that lead
 * to those scores. It return the best scored found in the matrix */
float fill_matrix(int ref_len, int query_len){
    register int row, col;
    float best_score = 0.0;
    float new_score;
    for (row = 1; row < query_len +1; row++)
        for (col = 1; col < ref_len +1; col++){
            new_score = score(row, col);
            best_score = (new_score > best_score)? new_score : best_score;
        }
    return best_score;
}


/* score:
 * for the specifed row and col in sw_matrix
 * assign the best  */
float score(int row, int col){
    SW_entry best_choice;
    SW_entry choices[4] = 
                       {null_entry,
                        score_mm(row, col),
                        score_insert(row, col),
                        score_delete(row, col)};
    best_choice = max(choices, 4);
    sw_matrix[row][col] = best_choice;
    return best_choice.score;
}

SW_entry score_mm(int row, int col){
    float prefix_score = sw_matrix[row-1][col-1].score;
    SW_entry mm_entry;
    if (g_ref[col -1] == g_query[row -1]){
        mm_entry.score = prefix_score + match_score;
        mm_entry.decision = 'M';
    } else{
        mm_entry.score = prefix_score + mismatch_penalty;
        mm_entry.decision = 'm';
    };
    return mm_entry;
}

SW_entry score_insert(int row, int col){
    int gap_len;
    SW_entry insert_entry = {0.0, 'I'};
    float new_score;
    for (gap_len = col; gap_len > 0; gap_len--){
        new_score = sw_matrix[row][col-gap_len].score + penalise_gap(gap_len);
        if (new_score > insert_entry.score)
            insert_entry.score = new_score;
    };
    return insert_entry;
}

SW_entry score_delete(int row, int col){
    int gap_len;
    SW_entry delete_entry = {0.0, 'D'};
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
SW_entry max(SW_entry list[], int list_len){
    SW_entry maximum = list[0]; // assume list_len > 0
    register int i;
    for (i = 1; i < list_len; i++)
        if (list[i].score > maximum.score)
            maximum = list[i];
    return maximum;
}


void print_sw_matrix(int nrow, int ncol){
    int row, col;
    SW_entry entry;
    for (row = 0; row < nrow; row++){
        for (col = 0; col < ncol; col++){
            entry = sw_matrix[row][col];
            printf("%.2f%c ", entry.score, entry.decision);
        };
        putchar('\n');
    }
}

void print_interaction_matrix(int nrow, int ncol){
    register int row, col;
    float max_interaction, temp;
    float average_interaction;
    for (row = 0; row < nrow; row++){
        printf("%60s ", pool[row]);
        max_interaction = 0.0;
        for (col = 0; col < ncol; col++){
            temp = interaction_matrix[row][col];
            printf("%.2f ", temp);
            if (temp > max_interaction)
                max_interaction = temp;
        }
        average_interaction = mean(interaction_matrix[row], ncol);
        printf("\t max= %.2f \tmean= %.2f\n", max_interaction, average_interaction);
    }
}

char complement(char base){
    base = toupper(base);
    switch (base){
        case 'A':
            return 'T';
            break;
        case 'T':
            return 'A';
            break;
        case 'G':
            return 'C';
            break;
        case 'C':
            return 'G';
            break;
        default:
            return 'N';
            break;
    }
}


char *rev_complement(char *seq){
    int len = strlen(seq);
    static char result[MAX_SEQ_LEN];
    int i;
    for (i = 0; i < len; i++)
        result[i] = complement(seq[len-i-1]);
    result[len] = '\0';
    return result;
}

float mean(float num_list[], int list_len){
    float result;
    register int i;
    for (i = 0; i < list_len; i++)
        result += num_list[i];
    return result / (float) list_len;
}

char *trim_whitespace(char *input){
    while (isspace(*input)) 
        input++;
    int len = strlen(input);
    char *endpointer = input + len -1;
    while (isspace(*endpointer) && endpointer != input)
        endpointer--;
    if (input + len -1 != endpointer) // if trailing space exist
        *(endpointer +1) = '\0'; // terminate string before space
    return input;
}

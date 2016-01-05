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
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include "swinc.h"

static SW_entry *init_matrix(int nrow, int ncol);
static void add_matrix_entry_repr(char *matrix_repr, Decision_Record record);

/* main():
 * read commandline parameters, obtain the primers in the primer pool
 * recorded in the specified file and then print out the interaction matrix
 * together with max_interaction and mean_interaction information */
int main(int argc, char **argv){
    User_Inputs user_inputs = parse_args(argc, argv);
    swalign(user_inputs.ref, user_inputs.query, user_inputs.scoring_param);
    int ref_len = strlen(user_inputs.ref);
    int query_len = strlen(user_inputs.query);
    printf("the alignment matrix is:\n");
    print_sw_matrix(query_len +1, ref_len+1);
    print_alignment(query_len +1, ref_len+1);
    return 0;
}


User_Inputs parse_args(int argc, char **argv)
{
    User_Inputs user_inputs;
    user_inputs.ref = "";
    user_inputs.query = "";
    user_inputs.primer_filename = "";
    user_inputs.verbose_flag = 0;
    user_inputs.scoring_param.match_score = DEFAULT_MATCH_SCORE;
    user_inputs.scoring_param.mismatch_penalty = DEFAULT_MISMATCH_PENALTY;
    user_inputs.scoring_param.gap_open_penalty = DEFAULT_GAP_OPEN_PENALTY;
    user_inputs.scoring_param.gap_extension_penalty = DEFAULT_GAP_EXTENSION_PENALTY;

    int opt;
    while ((opt = getopt(argc, argv, "rqfmxpe:")) != -1)
    {
        switch (opt)
        {
            case 'r':
                user_inputs.ref = optarg;
                break;
            case 'q':
                user_inputs.query = optarg;
                break;
            case 'f':
                user_inputs.primer_filename = optarg;
                break;
            case 'm':
                user_inputs.scoring_param.match_score = atof(optarg);
                break;
            case 'x':
                user_inputs.scoring_param.mismatch_penalty = atof(optarg);
                break;
            case 'p':
                user_inputs.scoring_param.gap_open_penalty = atof(optarg);
                break;
            case 'e':
                user_inputs.scoring_param.gap_extension_penalty = atof(optarg);
                break;
            case 'v':
                user_inputs.verbose_flag = 1;
                break;
            case '?':
                fprintf(stderr, 
                        "option %c isn't defined or missing its argument",
                        optopt);
                break;
            default:
                exit(EXIT_FAILURE);
                break;
        }
    }
    return user_inputs;
}







/******************* pool alignment routines *******************************/
/* align_pool:
 * given an array of strings, pool, return a matrix of size
 * |pool|x|pool| with matrix[i][j] representing the alignment score
 * between string_i and string_j (a complete graph).
 * The matrix will be symmetric and the diagonal element will be 0.0, i.e.
 * by default we don't want information about self-alignment.*/

void align_pool(int pool_size)
{
    int i, j;
    for (i = 0; i < pool_size; i++)
    {
        for (j = 0; j < pool_size; j++)
        {
            char *this_query = rev_complement(pool[j], KMER_SIZE);
            interaction_matrix[i][j] = (i==j)? 
                                       0.0 : 
                                       swalign(pool[i], this_query);
        }
    }
}


int get_primers(char *filename)
{
    FILE *file_handle = fopen(filename, "r");
    int count = 0;
    char buff[MAX_SEQ_LEN];
    char *temp;
    while ((temp = fgets(buff, MAX_SEQ_LEN, file_handle)) != NULL)
    {
        temp = trim_whitespace(temp); // remove left and right whitespaces
        if (count <= MAX_POOL_SIZE)
        {
            strcpy(pool[count], temp);
        }else
        {
            printf("Something's wrong");
        }
        count++;
    }
    fclose(file_handle);
    return count;
}

/******************** SW alignment algorithm ****************************/

/* swalign: produce the best alignment score for the 2 input string */
float verbose_swalign(User_Inputs user_inputs)
{
    int nrow = strlen(query) +1;
    int ncol = strlen(ref) +1;
    SW_entry sw_matrix[nrow][ncol];
    float best_score = fill_matrix(sw_matrix, 
                                   user_inputs.ref, 
                                   user_inputs.query, 
                                   user_inputs.score_param);
    // !! print things here
    return best_score;
}

/* swalign:
 * a non-verbose sw aligner that return only the score
 * of the best alignment.
 */
float swalign(char *ref, char *query, Score_Param score_param)
{
    int nrow = strlen(query) +1;
    int ncol = strlen(ref) +1;
    sw_entry sw_matrix[nrow][ncol];
    return fill_matrix(sw_matrix, ref, query, score_param);
}

/* fill_matrix:
 * fill the matrix with the the best scores and decisions that lead
 * to those scores. It return the best scored found in the matrix */
float fill_matrix(SW_entry **sw_matrix, 
                  char *ref, char *query, 
                  Score_Param score_param)
{
    register int row, col;
    float best_score = 0.0;
    float new_score;
    ref_len = strlen(ref);
    query_len = strlen(query);
    int nrow = query_len +1;
    int ncol = ref_len +1;
    //initiate the first row and first col to null entries.
    init_matrix(sw_matrix, nrow, ncol);
    for (row = 1; row < nrow; row++)
    {
        for (col = 1; col < ncol; col++)
        {
            new_score = score(sw_matrix, ref, query, row, col, score_param);
            best_score = (new_score > best_score)? new_score : best_score;
        }
    }
    return best_score;
}


/* init_matrix:
 * return a sw_matrix with the desired number of rows and cols with
 * all of first row, first column initialised to null entries
 */
static SW_entry *init_matrix(SW_entry **sw_matrix, int nrow, int ncol)
{
    // initialise the matrix first row and first col to null.
    Decision_Record null_decision_record = {0, "TT"};
    SW_entry null_entry = {
                           null_decision_record,
                           null_decision_record,
                           null_decision_record
                          };
    for (row = 0, col = 0; col < ncol; col++)
    {
        initiated_matrix[row][col] = null_entry;
    }
    for (col = 0, row = 0; row < nrow; row++)
    {
        initiated_matrix[row][col] = null_entry;
    }
    
    return initiated_matrix;
}



/* score:
 * for the specifed row and col in sw_matrix
 * assign the best  */
float score(SW_entry **sw_matrix, char *ref, char *query, 
            int row, int col, Score_Param score_param)
{
    SW_entry new_entry = 
                       {
                           score_match_mismatch(sw_matrix, 
                                                ref, query, 
                                                row, col, 
                                                score_param),
                           score_insert(sw_matrix, row, col, score_param),
                           score_delete(sw_matrix, row, col, score_param)
                       };
    sw_matrix[row][col] = new_entry;
    return max_record(new_entry, 3).score;
    /* there're 3 decision records in any sw_matrix entry
     * the best score is return and kept track of so that
     * we dont have to walk through 3mn decisions in mn entries to 
     * track down the alignment score for a[1,m], b[1,n].*/
}

/* score_mm: score mismatch. Given the position of 
 * an entry determine the score of the entry if the
 * alignment arrive at the entry through a 
 * match or mismatch (ie a diagonal movment
 * in the matrix */
Decision_Record score_match_mismatch(SW_entry **sw_matrix, 
                                     char *ref, char *query, 
                                     int row, int col, 
                                     Score_Param score_param)
{
    SW_entry previous_entry = sw_matrix[row-1][col-1];
    Decision_Record max_prev_record = max_record(previous_entry);
    char continued_from = max_prev_record.decision[1];
    Decision_Record match_mismatch_record;
    if (g_ref[col -1] == g_query[row -1])
    {
        match_mismatch_record.score = max_prev_record.score + \
                                      score_param.match_score;
        match_mismatch_record.decision = {continued_from, 'M', '\0'};
    } else
    {
        match_mismatch_record.score = max_prev_record.score + \
                                      score_param.mismatch_penalty;
        match_mismatch_record.decision = {continued_from, 'X', '\0'};
    }
    return mm_entry;
}

/* score_insert: compute the score at the
 * given position in the sw_matrix if 
 * the arrival at the pos is by an insertion
 * (a horizontal movement in the matrix) */

/* current insert can be a continuation of previous insert or
 * from a previous match/mismatch situation. 
 * If continue from insert, the penalty is for extending gap.
 * If it is from match/mismatch, the penalty is for opening the gap.
 * 
 * The same is for deletion.
 */
Decision_Record score_insert(SW_entry **sw_matrix, 
                             int row, int col, 
                             Score_Param score_param)
{
    SW_entry prev_entry = sw_matrix[row -1][col];
    // inspect the continuation from previous insertion first
    Decision_Record insert_record = {prev_entry.insert_record.score + \
                                       score_param.gap_extension_penalty,
                                       "II"};
    // now check if coming from match mismatch is better
    from_match_mismatch_score = prev_entry.match_record.score + \
                                score_param.gap_open_penalty;
    if (from_match_mismatch_score > insert_record.score)
    {
        insert_record.score = from_match_mismatch_score;
        insert_record.decision = {prev_entry.match_record.decision[1], 'I', '\0'};
    }
    return insert_record;
}


/* score_delete: compute the score at the 
 * given position in the sw_matrix if the
 * one arrive at the position by deletion
 * (a vertical movement in the matrix) */
Decision_Record score_delete(SW_entry **sw_matrix, 
                             int row, int col,
                             Score_Param score_param)
{
    SW_entry prev_entry = sw_entry[row][col-1];
    // inspect continuation from previous deletion first.
    Decision_Record delete_record = {prev_entyr.delete_record.score + \
                                     score_param.gap_extension_penalty,
                                     "DD"};
    from_match_mismatch_score = prev_entry.match_record.score + \
                                score_param.gap_open_penalty;
    if (from_match_mismatch_score > delete_record.score)
    {
        delete_record.score = from_match_mismatch_score;
        delete_record.decision = {prev_entry.match_record.decision[1], 'D', '\0'};
    }
    return delete_record;
}






/***** Utilities *********/

/* max_record:
 * given a list of Decision_Record type, return the one with the
 * maximum score. */
Decision_Record max_record(Decision_Record records[], int num_record)
{
    Decision_Record maximum = records[0];
    register int i;
    for (i = 1; i < num_record; i++)
    {
        if (records[i].score > maximum.score)
        {
            maximum = list[i]
        }
    }
    return maximum;
}


/* print_sw_matrix:
 * output a representation of sw_matrix.
 * Since each entry consist of 3 separate decision
 * we represent the whole matrix as 3 separate matrices
 * each for insert, delete and match separately.
 */
void print_sw_matrix(SW_entry **sw_matrix, int nrow, int ncol)
{
    int row, col;
    char *match_matrix = "";
    char *insert_matrrix = "";
    char *delete_matrix = "";
    char *entry_holder;
    SW_entry entry;
    for (row = 0; row < nrow; row++)
    {
        for (col = 0; col < ncol; col++)
        {
            entry = sw_matrix[row][col];
            add_matrix_entry_repr(match_matrix, entry.match_record);
            add_matrix_entry_repr(insert_matrix, entry.insert_record);
            add_matrix_entry_repr(delete_matrix, entry.delete_record);
        }
        strcat(match_matrix, "\n");
        strcat(insert_matrix, "\n");
        strcat(delete_matrix, "\n");
    }
    printf("Matrix for match/mismatch records\n%s", match_matrix);
    printf("Matrix for insert records\n%s", insert_matrix);
    printf("Matrix for delete records\n%s", delete_matrix);

}

static void add_matrix_entry_repr(char *matrix_repr, Decision_Record record)
{
    char *entry_holder;
    sprintf(entry_holder, "%.2f%c ", record.score, record.decision);
    strcat(matrix_repr, entry_holder);
}


void print_alignment(SW_entry **sw_matrix, 
                     char *ref, char *query)
{
    int nrow = strlen(query) +1;
    int ncol = strlen(ref) +1;
    int *best_entry_pos = best_entry(nrow, ncol);
    int row = best_entry_pos[0];
    int col = best_entry_pos[1];
    int num_insertion = 0;
    int num_deletion = 0;
    char decision;
    char *ref_string = "";
    char *query_string = "";
    SW_entry entry = sw_matrix[row][col];
    float best_score = entry.score;

    printf("Reference sequence = %s\n", ref);
    printf("Query sequence = %s\n", query);
    while ((decision = entry.decision) != "TT")
    { // not terminated yet
        switch (decision)
        {
            case ('M'): case ('m'):
                ref_string = prepend_char(ref_string, ref[col-1]);
                query_string = prepend_char(query_string, query[row-1]);
                row--;
                col--;
                break;
            case ('I'):
                ref_string = prepend_char(ref_string, '-');
                query_string = prepend_char(query_string, query[row-1]);
                row--;
                num_deletion++;
                break;
            case ('D'):
                ref_string = prepend_char(ref_string, ref[col-1]);
                query_string = prepend_char(query_string, '-');
                col--;
                num_insertion++;
                break;
        }
        entry = sw_matrix[row][col];
    }
    printf("Alignment score: %.2f\n", best_score);
    printf("Reference position: [%d, %lu)\n", 
            col, col + strlen(ref_string) - num_deletion);
    printf("Query position: [%d, %lu)\n", 
            row, row + strlen(query_string) - num_insertion);
    printf("%s\n",ref_string);
    printf("%s\n",query_string);
}

/* prepend_char: insert character c at the begining of string */
char *prepend_char(char *string, char c)
{
    char *result = malloc(strlen(string) +1 +1);
    assert(result != NULL);
    result[0] = c;
    strcpy(result+1, string);
    return result;
}

/* best_entry: search through sw_matrix and 
 * return the {row, col} of the best score entry
 */
int *best_entry(SW_entry **sw_entry, int nrow, int ncol)
{
    int *coord;
    coord = (int *) malloc(2 * sizeof(int));
    assert(coord != NULL);
    coord[0] = 0;
    coord[1] = 0;
    int row, col;
    int best_score = 0;
    for (row = 0; row < nrow; row++)
    {
        for (col = 0; col < ncol; col++)
        {// bottom rightmost maximum will be chosen among equal maxima
            if (sw_matrix[row][col].score >= best_score)
            {
                best_score = sw_matrix[row][col].score;
                coord[0] = row;
                coord[1] = col;
            }
        }
    }
    return coord;
}


void print_interaction_matrix(int nrow, int ncol)
{
    register int row, col;
    float max_interaction, temp;
    float average_interaction;
    for (row = 0; row < nrow; row++)
    {
        printf("%60s ", pool[row]);
        max_interaction = 0.0;
        for (col = 0; col < ncol; col++)
        {
            temp = interaction_matrix[row][col];
            printf("%.2f ", temp);
            if (temp > max_interaction)
            {
                max_interaction = temp;
            }
        }
        average_interaction = mean(interaction_matrix[row], ncol);
        printf("\t max= %.2f \tmean= %.2f\n", max_interaction, average_interaction);
        fprintf(outfile, "%.4f\t%.4f\n", max_interaction, average_interaction);
    }
}

/* complement:
 * given a nucleotide base letter,
 * return its complement in upper case.
 * return '\0' otherwise. */
char complement(char base){
    base = toupper(base);
    switch (base)
    {
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

/* rev_complement:
 * return the reverse complement of the inpput seq but we
 * only take the first result_len bases of the reverse
 * complement sequence */
char *rev_complement(char *seq, int result_len)
{
    char *result = malloc(sizeof(char) * (result_len+1));
    if (result == NULL)
    {
        error_handle(ERROR_MEM_ALLOC);
        exit(ERROR_MEM_ALLOC);
    }
    // if result_len should be smaller than seq_len
    int seq_len = strlen(seq);
    result_len = (seq_len < result_len)? seq_len : result_len;
    int i;
    for (i = 0; i < result_len; i++)
    {
        result[i] = complement(seq[seq_len-i-1]);
    }
    result[result_len] = '\0';
    return result;
}

/* mean:
 * given a list of float of length list_len,
 * compute and return their mean.*/
float mean(float num_list[], int list_len)
{
    float result;
    register int i;
    for (i = 0; i < list_len; i++)
    {
        result += num_list[i];
    }
    return result / (float) list_len;
}

/*trim_whitespace:
 * remove the whitespaces from the beginning and
 * end of input */
char *trim_whitespace(char *input)
{
    while (isspace(*input)) input++;
    int len = strlen(input);
    char *endpointer = input + len -1;
    while (isspace(*endpointer) && endpointer != input) 
        endpointer--;
    if (input + len -1 != endpointer)
    { // if trailing space exist
        *(endpointer +1) = '\0'; // terminate string before space
    }
    return input;
}





/********* Unused or depreciated ***************/
/* max_sw_entry:
 * return the entry in list[] with maximal score */
SW_entry max_sw_entry(SW_entry list[], int list_len)
{
    SW_entry maximum = list[0]; // assume list_len > 0
    register int i;
    for (i = 1; i < list_len; i++)
    {
        if (list[i].match_pair.score > maximum.score
            || list[i].delete_pair.score > maximum.
                )
        {
            maximum = list[i];
        }
    }
    return maximum;
}

/* penalise_gap: return the score given the gap length
 * according to affine gap scoring scheme */
float penalise_gap(int gap_len)
{
    float penalty = gap_open + gap_extension * gap_len;
    return penalty;
}

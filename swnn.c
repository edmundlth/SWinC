#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "swnn.h"

static SW_Entry **initialise_matrix(int nrow, int ncol);
static Decision_Record best_record(Decision_Record records[], int nrecord);

int main()
{
    float score = swnnalign("AAA",
                            "TTT");
    printf("%.2f\n", score);
    return 0;
}


/************************** ALIGNMENT ROUTINES ******************************/

/* swnnalign:
 * align ref and query string using swnn and return the 
 * score of the best alignment. */
float swnnalign(char *ref, char *query)
{
    float best_score = 0.0; // return value
    // Imagine the matrix with reference layout horizontally on top
    // and query running vertically downwards.
    // We also need to add extra 1 row and col to the matrix.
    int nrow = strlen(query) +1;
    int ncol = strlen(ref) +1;
    // reserve memory for sw_matrix (uninitialised)
    // This matrix will be pass around by all the scoring
    // routines.
    // The matrix will be initialised with all of first
    // row and first column being the null entry
    // { {0, STOP, STOP, 0}, ... }
    SW_Entry **sw_matrix = initialise_matrix(nrow, ncol);
    // Now fill up the matrix.
    register int row, col;
    for (row = 1; row < nrow; row++)
    {
        for (col = 1; col < ncol; col++)
        {
            sw_matrix[row][col] = compute_entry(sw_matrix, 
                                                row, col, 
                                                ref, query);
        }
    }
    best_score = find_best_score(sw_matrix, nrow, ncol);
    return best_score;
}

SW_Entry compute_entry(SW_Entry **sw_matrix,
                       int row, int col,
                       char *ref, char *query)
{
    SW_Entry entry;
    entry.match = score_match(sw_matrix,
                              row, col,
                              ref, query);
    entry.insertion = score_insertion(sw_matrix,
                                      row, col,
                                      ref, query);
    entry.deletion = score_deletion(sw_matrix,
                                    row, col,
                                    ref, query);
    return entry;
}



float find_best_score(SW_Entry **sw_matrix, int nrow, int ncol)
{
    register int row, col;
    float lowest_delG = 0.0;
    float new_delG;
    SW_Entry current_entry;
    for (row = 0; row < nrow; row++)
    {
        for (col = 0; col < ncol; col++)
        {
            current_entry = sw_matrix[row][col];
            Decision_Record three_options[] = {current_entry.match,
                                               current_entry.insertion,
                                               current_entry.deletion};
            new_delG = best_record(three_options, 3).delG;
            lowest_delG = (new_delG < lowest_delG) ? 
                           new_delG : lowest_delG;
        }
    }
    return lowest_delG;
}







/************************** SCORING ROUTINES ******************************/

Decision_Record score_match(SW_Entry **sw_matrix,
                            int row, int col,
                            char *ref, char *query)
{
    SW_Entry prev_entry = sw_matrix[row -1][col -1];
    Decision_Record match_match; // read match then match
    Decision_Record insertion_match; 
    // read insertion then match, 
    // ie current match comes is continued from previous insertion
    Decision_Record  deletion_match; // read deletion then match

    /* Handles match_match, cases:
     * double match   : add on nearest neighbour data
     * double mismatch: !! handle bulge loop
     * match mismatch : add on nearest neighbour data */
    Neighbour nn_config = {ref[col-1], ref[col], 
                           query[row-1], query[row]};
    int is_current_complement = is_complement(nn_config.top3, 
                                              nn_config.bottom5);
    int is_previous_complement = is_complement(nn_config.top5, 
                                               nn_config.bottom3);

    match_match.previous_decision = (is_previous_complement) ? MATCH : MISMATCH;
    match_match.current_state = (is_current_complement) ? MATCH : MISMATCH;
    if (is_current_complement && is_previous_complement)
    {// the nicest case where we just zip the nearest neighbour delG value
        match_match.delG = prev_entry.match.delG + get_delG_internal(nn_config);
    } else if (is_current_complement && prev_entry.match.loop_len > 1)
    {//!! previous mismatch is part of internal loop that assume current complement
        match_match.delG = prev_entry.match.delG;
        match_match.loop_len = 0;
    }else if (is_current_complement && prev_entry.match.loop_len == 1)
    {// this is just a MXM situation, we have nearest neighbour data
        match_match.delG = prev_entry.match.delG + get_delG_internal(nn_config);
        match_match.loop_len = 0;
    } else if (is_previous_complement)
    {/* maybe just a mismatch or the begining of internal loop,
        we assume mismatch. If the next one is mismatch, it should break
        and extend */
        match_match.delG = prev_entry.match.delG + get_delG_internal(nn_config);
        match_match.loop_len = 1;
    } else
    {/* double mismatch must be internal loop of part of a larger loop,
        break previous loop and extend loop length by one */
        match_match.delG = prev_entry.match.delG + \
                           extend_internal_loop(prev_entry.match.loop_len);
        match_match.loop_len = prev_entry.match.loop_len +1;
    }

    // !! closing of a loop cost nothing?? Since the previous loop assume this match?
    insertion_match.delG = prev_entry.insertion.delG; 
    insertion_match.previous_decision = INSERTION;
    if (is_current_complement)
    { // if current base match up the loop len is 0
        insertion_match.current_state = MATCH;
        insertion_match.loop_len = 0;
    } else
    { // in a loop of len 1, possibly just mismatch or might be extended later
        insertion_match.current_state = MISMATCH;
        insertion_match.loop_len = 1;
    }

    // !! see prev 
    // !! This could have been merge into previous conditional statements
    // !! but decided against it for code clarity
    deletion_match.delG = prev_entry.deletion.delG;
    if (is_current_complement)
    {
        deletion_match.current_state = MATCH;
        deletion_match.loop_len = 0;
    } else
    {
        deletion_match.current_state = MISMATCH;
        deletion_match.loop_len = 1;
    }
    
    Decision_Record three_options[] = {match_match, 
                                       insertion_match, 
                                       deletion_match};
    return best_record(three_options, 3);
}




Decision_Record score_insertion(SW_Entry **sw_matrix, 
                                int row, int col,
                                char *ref, char *query)
{
    /* insertion is gap in reference, ie the col coord doesn't vary
     * (don't add more reference bases) while row coord vary 
     * (adding more query bases). */
    SW_Entry prev_entry = sw_matrix[row-1][col];
    //only 2 previous decision will lead to an insertion
    //since we don't allow immediate succesion of insertion and deletion
    Decision_Record match_insertion; // bulge loop of size 1
    Decision_Record insertion_insertion; // bulge loop of size prev_size+1

    // size 1 bulge is special. Need consider the 2 bases flanking the bulge
    Neighbour intervening_bases_config = {ref[col -1], ref[col +1],
                                          query[row -1], query[row]};
    match_insertion.delG = prev_entry.match.delG + \
                           size_1_bulge(intervening_bases_config);
    match_insertion.current_state = INSERTION;
    match_insertion.previous_decision = (is_complement(ref[col], query[row-1])) ?
                                        MATCH : MISMATCH;
    match_insertion.loop_len = 1;

    insertion_insertion.delG = prev_entry.insertion.delG + \
                               extend_bulge_loop(prev_entry.insertion.loop_len);
    insertion_insertion.current_state = INSERTION;
    insertion_insertion.previous_decision = INSERTION;
    insertion_insertion.loop_len = prev_entry.insertion.loop_len +1;

    Decision_Record two_options[] = {match_insertion, insertion_insertion};
    return best_record(two_options, 2);
}








/* !!! This is too much like score_insertion, basically repeated code
 */
Decision_Record score_deletion(SW_Entry **sw_matrix, 
                                int row, int col,
                                char *ref, char *query)
{
    /* deletion is gap in query, ie the row coord doesn't vary
     * (don't add more query base) while col coord vary 
     * (adding more reference base). */
    SW_Entry prev_entry = sw_matrix[row][col -1];
    //only 2 previous decision will lead to an deletion
    //since we don't allow immediate succesion of deletion and deletion
    Decision_Record match_deletion; // bulge loop of size 1
    Decision_Record deletion_deletion; // bulge loop of size prev_size+1

    // size 1 bulge is special. Need consider the 2 bases flanking the bulge
    Neighbour intervening_bases_config = {ref[col -1], ref[col],
                                          query[row -1], query[row +1]};
    match_deletion.delG = prev_entry.match.delG + \
                           size_1_bulge(intervening_bases_config);
    match_deletion.current_state = INSERTION;
    match_deletion.previous_decision = (is_complement(ref[row -1], query[col])) ? 
                                        MATCH : MISMATCH;
    match_deletion.loop_len = 1;

    deletion_deletion.delG = prev_entry.deletion.delG + \
                               extend_bulge_loop(prev_entry.deletion.loop_len);
    deletion_deletion.current_state = INSERTION;
    deletion_deletion.previous_decision = INSERTION;
    deletion_deletion.loop_len = prev_entry.deletion.loop_len +1;

    Decision_Record two_options[] = {match_deletion, deletion_deletion};
    return best_record(two_options, 2);
}




/********************** THERMODYNAMICS ROUTINES ****************************/


#define INTERNAL_A 0
#define INTERNAL_C 1
#define INTERNAL_G 2
#define INTERNAL_T 3
#define TERMINAL_DOT 0
#define TERMINAL_A 1
#define TERMINAL_C 2
#define TERMINAL_G 3
#define TERMINAL_T 4

#define NUM_SYS_BASE_INTERNAL 4
#define NUM_SYS_BASE_TERMINAL 5

int _digit_internal(char base)
{
    int digit;
    switch (base)
    {
        case 'A':
             digit = INTERNAL_A;
             break;
        case 'C':
             digit = INTERNAL_C;
             break;
        case 'G':
             digit = INTERNAL_G;
             break;
        case 'T':
             digit = INTERNAL_T;
             break;
        default:
        digit = -1;
    }
    return digit;
}

int _digit_terminal(char base)
{
    int digit;
    switch (base)
    {
        case 'A':
             digit = TERMINAL_A;
             break;
        case 'C':
             digit = TERMINAL_C;
             break;
        case 'G':
             digit = TERMINAL_G;
             break;
        case 'T':
             digit = TERMINAL_T;
             break;
        default:
        digit = -1;
    }
    return digit;
}


    
int _get_index_internal(Neighbour nn_config)
{
    int index = 0;
    printf("%c%c/%c%c\n", nn_config.top5, nn_config.top3, nn_config.bottom3, nn_config.bottom5);
    index += _digit_internal(nn_config.top5) * pow(NUM_SYS_BASE_INTERNAL, 3);
    index += _digit_internal(nn_config.top3) * pow(NUM_SYS_BASE_INTERNAL, 2);
    index += _digit_internal(nn_config.bottom3) * pow(NUM_SYS_BASE_INTERNAL, 1);
    index += _digit_internal(nn_config.bottom5) * pow(NUM_SYS_BASE_INTERNAL, 0);
    return index;
}

int _get_index_terminal(Neighbour nn_config)
{
    int index = 0;
    index += _digit_terminal(nn_config.top5) * pow(NUM_SYS_BASE_TERMINAL, 3);
    index += _digit_terminal(nn_config.top3) * pow(NUM_SYS_BASE_TERMINAL, 2);
    index += _digit_terminal(nn_config.bottom3) * pow(NUM_SYS_BASE_TERMINAL, 1);
    index += _digit_terminal(nn_config.bottom5) * pow(NUM_SYS_BASE_TERMINAL, 0);
    return index;
}


float get_delG_internal(Neighbour nn_config)
{
    int index = _get_index_internal(nn_config);
    Therm_Param record = nn_data_internal[index];
    return record.delH * 1000.0 - (Reaction_Temperature + ABSOLUTE_ZERO_OFFSET) * record.delS;
}


/* !! The idea for extension is add on gradient(prev_len). */
float extend_internal_loop(int previous_loop_len)
{

    return 0.0;
}

float extend_bulge_loop(int previous_loop_len)
{
    return 0.0;
}

/* size_1_bulge: */
float size_1_bulge(Neighbour intervening_bases_config)
{
    float bulge_penalty = 4.0;
    float AT_penalty = 1.0;
    return get_delG_internal(intervening_bases_config) + bulge_penalty + AT_penalty;
}



/************************** UTILITIES ROUTINES ******************************/

/* initialise_matrix:
 * given a sw_matrix of dimension nrow x nrow, populate all of its first
 * row and first col with null entry */
static SW_Entry **initialise_matrix(int nrow, int ncol)
{
    register int i;
    SW_Entry **sw_matrix = malloc(sizeof(SW_Entry *) * nrow);
    for (i = 0; i < nrow; i++)
    {
        sw_matrix[i] = malloc(sizeof(SW_Entry) * ncol);
    }
    if (sw_matrix == NULL)
    {
        fprintf(stderr, "swnn: memory allocation error");
        exit(EXIT_FAILURE);
    }
    Decision_Record null_decision = {0.0, STOP, STOP, 0};
    SW_Entry null_entry = {
                           null_decision,
                           null_decision,
                           null_decision
                          };
    register int row, col;
    for (row = 0; row < nrow; row++)
    {
        sw_matrix[row][0] = null_entry;
    }
    for (col = 0; col < ncol; col++)
    {
        sw_matrix[0][col] = null_entry;
    }
    return sw_matrix;
}



/* complement: return the complement of given base
 * character in upper case */
char complement(char base)
{
    switch (toupper(base))
    {
        case 'A':
            return 'T';
            break;
        case 'T':
            return 'A';
            break;
        case 'C':
            return 'G';
            break;
        case 'G':
            return 'C';
            break;
        default:
            return '\0';
            break;
    }
}

/* is_complement: True if base1 is the complement of base2.
 * Not case sensitive */
int is_complement(char base1, char base2)
{
    if (complement(toupper(base1)) == toupper(base2))
    {
        return TRUE;
    } else
    {
        return FALSE;
    }
}


/* best_record: private routine that select the best decision
 * among the list of decision record.
 * !! best is currently defined as lowest delG value. */
static Decision_Record best_record(Decision_Record records[], int nrecord)
{
    Decision_Record best_record = records[0]; // assume at least 1 in list
    register int i;
    for (i = 0; i < nrecord; i++)
    {
        if (records[i].delG < best_record.delG)
        {//!! decide that in swnn, minimising delG is more meaningful
            best_record = records[i];
        }
    }
    return best_record;
}




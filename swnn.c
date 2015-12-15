#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "swnn.h"

static void initialise_matrix(SW_Entry **sw_matrix, int nrow, int ncol);
static Decision_Record best_record(Decision_Record records[], int nrecord);

int main()
{
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
    register int row, col;
    // reserve memory for sw_matrix (uninitialised)
    // This matrix will be pass around by all the scoring
    // routines.
    SW_Entry sw_matrix[nrow][ncol];
    initialise_matrix(sw_matrix, nrow, ncol);
    for (row = 1; row < nrow; row++)
    {
        for (col = 1; col < ncol; col++)
        {
            sw_matrix[row][col] = compute_entry(sw_matrix, 
                                                row, col, 
                                                ref, query);
        }
    }
    best_score = find_best_score(sw_matrix);
    return best_score;
}

SW_Entry compute_entry(SW_entry **sw_matrix,
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
        match_match.delG = prev_entry.match.delG + get_delG(nn_config);
    } else if (is_current_complement && prev_entry.match.loop_len > 1)
    {//!! previous mismatch is part of internal loop that assume current complement
        match_match.delG = prev_entry.match.delG;
        match_match.loop_len = 0;
    }else if (is_current_complement && prev_entry.match.loop_len == 1)
    {// this is just a MXM situation, we have nearest neighbour data
        match_match.delG = prev_entry.match.delG + get_delG(nn_config);
        match_match.loop_len = 0;
    } else if (is_previous_complement)
    {/* maybe just a mismatch or the begining of internal loop,
        we assume mismatch. If the next one is mismatch, it should break
        and extend */
        match_match.delG = prev_entry.match.delG + get_delG(nn_config);
        match_match.loop_len = 1
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

    int is_previous_complement = is_complement(match_insertion_config.top5,
                                               match_insertion_config.bottom3);


    // size 1 bulge is special. Need consider the 2 bases flanking the bulge
    Neighbour intervening_bases_config = {ref[col -1], ref[col +1],
                                          query[row -1], query[row]};
    match_insertion.delG = prev_entry.match.delG + \
                           size_1_bulge(intervening_bases_config);
    match_insertion.current_state = INSERTION;
    match_insertion.previous_decision = (is_previous_complement) ? MATCH : MISMATCH;
    match_insertion.loop_len = 1;

    insertion_insertion.delG = prev_entry.insertion.delG + \
                               extend_bulge_loop(prev_entry.insertion.loop_len);
    insertion_insertion.current_state = INSERTION;
    insertion_insertion.previous_state = INSERTION;
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

    int is_previous_complement = is_complement(match_deletion_config.top5,
                                               match_deletion_config.bottom3);


    // size 1 bulge is special. Need consider the 2 bases flanking the bulge
    Neighbour intervening_bases_config = {ref[col -1], ref[col],
                                          query[row -1], query[row +1]};
    match_deletion.delG = prev_entry.match.delG + \
                           size_1_bulge(intervening_bases_config);
    match_deletion.current_state = INSERTION;
    match_deletion.previous_decision = (is_previous_complement) ? MATCH : MISMATCH;
    match_deletion.loop_len = 1;

    deletion_deletion.delG = prev_entry.deletion.delG + \
                               extend_bulge_loop(prev_entry.deletion.loop_len);
    deletion_deletion.current_state = INSERTION;
    deletion_deletion.previous_state = INSERTION;
    deletion_deletion.loop_len = prev_entry.deletion.loop_len +1;

    Decision_Record two_options[] = {match_deletion, deletion_deletion};
    return best_record(two_options, 2);
}




/********************** THERMODYNAMICS ROUTINES ****************************/

float get_delG(Neighbour nn_config);
float size_1_bulge(Neighbour intervening_bases_config);

/* !! The idea for extension is add on gradient(prev_len). */
float extend_internal_loop(int previous_loop_len);
float extend_bulge_loop(int previous_loop_len);

/* size_1_bulge: */
float size_1_bulge(Neighbour intervening_bases_config)
{
    return get_delG(intervening_bases_config) + bulge_penalty + AT_penalty;
}



/************************** UTILITIES ROUTINES ******************************/

/* initialise_matrix:
 * given a sw_matrix of dimension nrow x nrow, populate all of its first
 * row and first col with null entry */
static void initialise_matrix(SW_Entry **sw_matrix, int nrow, int ncol)
{
    Decision_Record null_decision = {0.0, STOP, STOP};
    SW_Entry null_entry = {
                           null_decision,
                           null_decision,
                           null_decision
                          };
    register int row, col;
    for (row = 0, col = 0; row < nrow; row++)
    {
        sw_matrix[row][col] = null_entry;
    }
    for (col = 0, row = 0; col < ncol; col++)
    {
        sw_matrix[row][col] = null_entry;
    }
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
            return NULL;
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




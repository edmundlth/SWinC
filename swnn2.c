/* Implementation of secondary structure prediction algorithm using
 * idea from Smith-Waterman Alignement algorithm and Nearest neighbour
 * thermodynamics. 
 */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swnn2.h"

static SW_Entry **initialise_matrix(int nrow, int ncol);
static Decision_Record best_record(Decision_Record records[], int nrecord);

int main()
{
    Neighbour nn_config = {'A', 'C', 'T','A'};
    Therm_Param from_record = nn_data_terminal[_get_index_terminal(nn_config)];
    printf("%s %f %f\n", from_record.neighbour, from_record.delH, from_record.delS);
    printf("delG = %f\n", from_record.delH * 1000.0 - (Reaction_Temperature + ABSOLUTE_ZERO_OFFSET) * from_record.delS);
    printf("delG_function = %f\n", get_delG_terminal(nn_config));
    return 0;
}


/************************** ALIGNMENT ROUTINES ******************************/

/* duplex_matrix:
 * Given a reference sequence (sense sequence)
 * and a query sequence (antisense sequence)
 * return a DP matrix recording all the decisions and scores. 
 */
SW_Entry **duplex_matrix(char *ref, char *query)
{
    int nrow = strlen(query) +1;
    // the matrix layout:
    // reference is on the horizontal while
    // query is at the bottom. 
    // We add an extra row for initialisation.
    int ncol = strlen(ref) +1;
    // The matrix will be initialise with all
    // of first row and first col being null entry
    // eg. the element [0][0] will be 
    // {{0, STOP, STOP, 0}, 
    //  {0, STOP, STOP, 0}, 
    //  {0, STOP, STOP, 0}, 
    //  {0, STOP, STOP,0}}
    //  !!! THIS NEED TO BE CHANGED!!! INITIATION IS DEPENDEND ON SEQUENCE
    //  !!! CONTEXT, WE NEED TO ACCOUNT FOR INTIATION_ENTROPY AND 
    //  !!! TERMINAL EFFECTS (DANGLING ENDS AND MISMATCHES).
    SW_Entry **sw_matrix = initialise_matrix(nrow, ncol);
    // now we fill up the matrix
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
    return sw_matrix;
}





SW_Entry compute_entry(SW_Entry **sw_matrix,
                       int row, int col,
                       char *ref, char *query)
{
    SW_Entry entry;
    entry.bind = score_bind(sw_matrix,
                              row, col,
                              ref, query);
    entry.top_bulge = score_top_bulge(sw_matrix,
                                      row, col,
                                      ref, query);
    entry.bottom_bulge = score_bottom_bulge(sw_matrix,
                                            row, col,
                                            ref, query);
    entry.stop = score_stop(sw_matrix,
                            row, col,
                            ref, query);
    return entry;
}



Coord find_best_entry_coord(SW_Entry **sw_matrix, int nrow, int ncol)
{
    register int row, col;
    float lowest_delG = 0.0;
    float new_delG;
    Coord best_coord = { 0, 0, STOP};
    SW_Entry current_entry;
    Decision_Record record;
    for (row = 0; row < nrow; row++)
    {
        for (col = 0; col < ncol; col++)
        {
            current_entry = sw_matrix[row][col];
            Decision_Record four_options[] = {current_entry.bind,
                                               current_entry.top_bulge,
                                               current_entry.bottom_bulge,
                                               current_entry.stop};
            record = best_record(four_options, 4);
            new_delG = record.delG;
            if (new_delG < lowest_delG)
            {
                lowest_delG = new_delG;
                best_coord.row = row;
                best_coord.col = col;
                best_coord.current_decision = record.current_decision;
            }
        }
    }
    return best_coord;
}

/************************** SCORING ROUTINES *****************************
 * In sw_matrix, each entry contain 4 "Decision Records" we call 
 * "current decision": They are: bind, top_bulge, bottom_bulge, stop. 
 * 
 * These current decisions are continuation from a "previous decision",
 * or decision from previous entry's "current decision".
 * But there're limit to which prev decision a particular type of current
 * decision can continue from:
 * A bind can continue from [bind, top_bulge, bottom_bulge]
 * A top_bulge can continue from [bind, top_bulge]
 * A bottom_bulge can continue from [bind, bottom_bulge]
 * A stop can continue from [bind, top_bulge, bottom_bulge]
 * Notice that nothing continue from a "stop" (as the name implied). 
 * 
 * A separate scoring function is written for each current decision, 
 * i.e. it takes 4 routines to "fill" the entry. 
 * Each function evaluate score of all legal continuation and write to 
 * the entry only the best scored one (best score == lowest delG).
 * Each function depend only on the value in the appropriate "previous entry",
 * the reference and the query sequences. 
 * !!! (There is implicit dependency on the reaction temperature and salt condition
 * which are set as global variables. This should be changed to a single 
 * structure "Reaction Condition" and pass it around)
 *
****************************************************************************/

/* score_bind:
 * evaluate the score (delG) of deciding on matching (or mismatching) current
 * bases as continuation from [bind, bottom_bulge, top_bulge]. 
 * Return the record of best continuation. 
 */
Decision_Record score_bind(SW_Entry **sw_matrix,
                           int row, int col,
                           char *ref, char *query)
{
    SW_Entry prev_entry = sw_matrix[row-1][col -1];
    Decision_Record prev_decision_record;
    char current_decision = (is_complement(query[row], ref[col])) ? 'M' : 'X';
    Decision_Record continue_from_bind = {0, MATCH, current_decision, 0, 0};
    Decision_Record continue_from_top_bulge = {0, TOP_BULGE, current_decision, 0, 0};
    Decision_Record continue_from_bottom_bulge = {0, BOTTOM_BULGE, current_decision, 0, 0};

    // handle continue from previous binding:
    // 2 cases,
    // - previous binding is a match: simply add on delG from match table
    // - previous binding is a mismatch: 2 subcases
    //     -> current binding is a match: another 2 subcases
    //         => previous mismatch is a single mismatch: delG from mismatch table
    //         => previous mismatch part of loop: loop penalty and close of the loop
    //     -> current binding is a mismatch: another 2 subcases
    //         => previous mismatch is single: backtrack the delG addition and add on loop penalty instead
    //         => previous mismatch part of loop: loop penalty
    prev_decision_record = prev_entry.bind;
    switch (prev_decision_record.current_decision)
    {
        case (MATCH): // continue from previous match
             {// do further zipping, internal delG handle both match and mismatch
                 Neighbour nn_config = {ref[col -1], ref[col], query[row -1], query[row]};
                 continue_from_bind.delG = prev_decision_record.delG + get_delG_internal(nn_config);
                 continue_from_bind.top_loop_len = (current_decision == 'M') ? 0 : 1;
                 continue_from_bind.bottom_loop_len = continue_from_bind.top_loop_len;
             }
             break;
        case (MISMATCH):
            {// continue from previous mismatch
               switch (current_decision)
               {
                   case (MATCH):
                       { // mismatch then match
                          if (prev_decision_record.top_loop_len == 1 && prev_decision_record.bottom_loop_len == 1)
                          {// do further zipping
                              Neighbour nn_config = {ref[col -1], ref[col], query[row -1], query[row]};
                              continue_from_bind.delG = prev_decision_record.delG + get_delG_internal(nn_config);
                              continue_from_bind.previous_decision = MISMATCH;
                          } else if (prev_decision_record.top_loop_len > 1 || prev_decision_record.bottom_loop_len > 1)
                          {// do nothing since previous internal loop calculation already assume this is match
                              continue_from_bind.delG = prev_decision_record.delG;
                              continue_from_bind.previous_decision = MISMATCH;
                          }
                       }
                       break;
                   case (MISMATCH): // mismatch and mismatch
                          {
                              continue_from_bind.previous_decision = MISMATCH;
                              if (prev_decision_record.top_loop_len == 1 && prev_decision_record.bottom_loop_len == 1)
                              {// backtrack delG addition and score loop
                                  Neighbour nn_config = {ref[col -1], ref[col], query[row -1], query[row]};
                                  continue_from_bind.top_loop_len = 2; // adding 1 to previous loop of 1
                                  continue_from_bind.bottom_loop_len = 2; // doing the same to bottom, since if a mismatch
                                  continue_from_bind.delG = prev_decision_record.delG \
                                                            - get_delG_internal(nn_config) \
                                                            + internal_loop_score(continue_from_bind.top_loop_len, 
                                                                                    continue_from_bind.bottom_loop_len);
                              } else if (prev_decision_record.top_loop_len > 1 || prev_decision_record.bottom_loop_len > 1)
                              {// no backtracking required, already in loop
                                  continue_from_bind.top_loop_len = prev_decision_record.top_loop_len +1;
                                  continue_from_bind.bottom_loop_len = prev_decision_record.bottom_loop_len +1;
                                  continue_from_bind.delG = prev_decision_record.delG \
                                                            + internal_loop_score(continue_from_bind.top_loop_len,
                                                                                    continue_from_bind.bottom_loop_len);
                              }
                          }
                          break;
                   default:
                   fprintf(stderr, "Neither MATCH nor MISMATCH at bind record");
                   exit(EXIT_FAILURE);
               }
            }
    }
    // now handle continuing from previous top_bulge:
    // previously, we add 1 to top_loop_len while none to bottom_loop_len,
    // so bottom_loop_len might potentially be zero. 
    // Now if current_decision is a MISMATCH, regardless of the size of
    // top and bottom loop, an addition of 1 to both send us into 
    // an internal loop situation. 
    // On the other hand, if current decision is a MATCH, 
    // we might remain in a bulge or continue a previous internal loop,
    // but since both case assume current base is a match, we do nothing
    // except carry the record over. 
    prev_decision_record = prev_entry.top_bulge;
    switch (current_decision)
    {
        case (MATCH):
             continue_from_top_bulge.delG = prev_decision_record.delG;
             break;
        case (MISMATCH):
             {
                 continue_from_top_bulge.top_loop_len = prev_decision_record.top_loop_len +1;
                 continue_from_top_bulge.bottom_loop_len = prev_decision_record.bottom_loop_len +1;
                 continue_from_top_bulge.delG = prev_decision_record.delG \
                                                + internal_loop_score(continue_from_top_bulge.top_loop_len,
                                                                      continue_from_top_bulge.bottom_loop_len);
             }
             break;
    }
    // now handle continuing from previous bottom_bulge,
    // the resoning is the same as continuing from top_bulge. 
    prev_decision_record = prev_entry.bottom_bulge;
    switch (current_decision)
    {
        case (MATCH):
             continue_from_bottom_bulge.delG = prev_decision_record.delG;
             break;
        case (MISMATCH):
             {
                 continue_from_bottom_bulge.top_loop_len = prev_decision_record.top_loop_len +1;
                 continue_from_bottom_bulge.bottom_loop_len = prev_decision_record.bottom_loop_len +1;
                 continue_from_bottom_bulge.delG = prev_decision_record.delG \
                                                + internal_loop_score(continue_from_bottom_bulge.top_loop_len,
                                                                      continue_from_bottom_bulge.bottom_loop_len);
             }
             break;
    }
    // we have evaluated and recorded all 3 continuation, 
    // now return the best one.
    Decision_Record three_continuation_records[] = {continue_from_bind,
                                                    continue_from_top_bulge,
                                                    continue_from_bottom_bulge};
    return best_record(three_continuation_records,3);
}

/* score_top_bulge: 
 * compute the 
 */
Decision_Record score_top_bulge(SW_Entry **sw_matrix, 
                                int row, int col,
                                char *ref, char *query)
{
    SW_Entry prev_entry = sw_matrix[row][col -1];
    Decision_Record previous_decision_record;
    Decision_Record continue_from_bind = {0, MATCH, TOP_BULGE, 0, 0};
    Decision_Record continue_from_top_bulge = {0, TOP_BULGE, TOP_BULGE, 0, 0};
    // handle continue from preivous match or mismatch: 2 cases
    // previous binding is a match: introduce a size1 bulge which is to be score
    //                              differently from bulge size > 1. Flanking
    //                              base contribution has to be included. 
    //                              Thus bases from 
    //                                  ref[col -1] ref[col +1]
    //                                  query[row] query[row +1]
    //                              to be included
    // previous binding is a mismatch: this is an internal loop, but there's 2 cases
    //     -> previous mismatch is single: need to backtrack previous delG addition
    //                                     then extend to form internal loop
    //     -> previous mismatch part of loop: continue the internal loop
    previous_decision_record = prev_entry.bind;
    if (previous_decision_record.current_decision == MATCH)
    {
        continue_from_bind.top_loop_len = 1; // increment from previous 0
        /* neighbour configuration:
         *     MBM
         *     m M    where B is at current position (col)
         *            m is at current row
         *     current ref base doesn't bind,
         *     previous ref base bind with current query base
         *     while next ref base bind with next query base */
        Neighbour nn_config = {ref[col -1], ref[col +1],
                              query[row], query[row +1]};
        continue_from_bind.delG = previous_decision_record.delG \
                                  + bulge_score(continue_from_bind.top_loop_len)
                                  + get_delG_internal(nn_config);
    } else if (previous_decision_record.current_decision == MISMATCH)
    {
        if (previous_decision_record.top_loop_len == 1 && previous_decision_record.bottom_loop_len == 1)
        {// need to backtrack previous delG addition
            continue_from_bind.top_loop_len = 2; // add one from previous single mismatch
            continue_from_bind.bottom_loop_len = 1; // from previous single mismatch

            /* neighbour:
             *     MXB
             *     Mx   where B is the current bulge, which is current col
             *          while x is current row.
             *     previous binding was treated as a mismatch, 
             *     but now it's part of a loop, the mismatch contribution
             *     needs to be backtracked. */
            Neighbour nn_config = {ref[col -2], ref[col -1],
                                  query[row -1], query[row]};
            continue_from_bind.delG = previous_decision_record.delG \
                                      - get_delG_internal(nn_config) \
                                      + internal_loop_score(continue_from_bind.top_loop_len,
                                                            continue_from_bind.bottom_loop_len);
        } else if (previous_decision_record.top_loop_len > 1 || previous_decision_record.bottom_loop_len > 1)
        {
            continue_from_bind.top_loop_len = previous_decision_record.top_loop_len +1;
            continue_from_bind.bottom_loop_len = previous_decision_record.bottom_loop_len;
            continue_from_bind.delG = previous_decision_record.delG \
                                      + internal_loop_score(continue_from_bind.top_loop_len,
                                                            continue_from_bind.bottom_loop_len);
        }
    }
    // now handle continue from previous top_bulge: 2 cases
    // previous bulge has size 1: need to backtrack the special size one intervening delG addition
    // previous bulge size > 1: simply extend bulge size.
    previous_decision_record = prev_entry.bottom_bulge;
    if (previous_decision_record.bottom_loop_len == 1 && previous_decision_record.top_loop_len == 0)
    {
        /* neighbour configuration:
         *     Mbm
         *     m M  where b is previous bulge.
         *          top "m" is current col which needs to be changed to bulge.
         *          bottom m is current row */
        Neighbour nn_config = {ref[col -2], ref[col],
                               query[row], query[row +1]};
        continue_from_top_bulge.top_loop_len = 2; // added 1 to prev len
        continue_from_top_bulge.delG = previous_decision_record.delG \
                                          - get_delG_internal(nn_config) \
                                          + internal_loop_score(continue_from_top_bulge.top_loop_len,
                                                                continue_from_top_bulge.bottom_loop_len);
    } else if (previous_decision_record.top_loop_len > 1 || previous_decision_record.bottom_loop_len > 0)
    {
        continue_from_top_bulge.top_loop_len = previous_decision_record.top_loop_len +1;
        continue_from_top_bulge.bottom_loop_len = previous_decision_record.bottom_loop_len;
        continue_from_top_bulge.delG = previous_decision_record.delG \
                                          + internal_loop_score(continue_from_top_bulge.top_loop_len,
                                                                continue_from_top_bulge.bottom_loop_len);
    }

    Decision_Record two_continuation_records[] = {continue_from_bind,
                                                continue_from_top_bulge};
    return best_record(two_continuation_records, 2);

}

// !!!!!!!!!! all the scoring function should be a function of "prev entry" only
// should just pass in "prev_entry" not the whole matrix. 
// this way, one should be able to handle both score_bulge in a single function.

/* score_bottom_bulge:
 * 
 */
Decision_Record score_bottom_bulge(SW_Entry **sw_matrix, 
                                int row, int col,
                                char *ref, char *query)
{
    SW_Entry prev_entry = sw_matrix[row -1][col];
    Decision_Record previous_decision_record;
    Decision_Record continue_from_bind = {0, MATCH, BOTTOM_BULGE, 0, 0};
    Decision_Record continue_from_bottom_bulge = {0, BOTTOM_BULGE, BOTTOM_BULGE, 0, 0};
    // handle continue from preivous match or mismatch: 2 cases
    // previous binding is a match: introduce a size1 bulge which is to be score
    //                              differently from bulge size > 1. Flanking
    //                              base contribution has to be included. 
    //                              Thus bases from 
    //                                  ref[col] ref[col +1]
    //                                  query[row-1] query[row +1]
    //                              to be included
    // previous binding is a mismatch: this is an internal loop, but there's 2 cases
    //     -> previous mismatch is single: need to backtrack previous delG addition
    //                                     then extend to form internal loop
    //     -> previous mismatch part of loop: continue the internal loop
    previous_decision_record = prev_entry.bind;
    if (previous_decision_record.current_decision == MATCH)
    {
        continue_from_bind.bottom_loop_len = 1;
        /* neighbour configuration:
         *     m M
         *     MBM    where B is at current position (row), m is at current row
         *     current query base doesn't bind, 
         *     previous query base bind with current ref base
         *     while next query base bind with next ref base */
        Neighbour nn_config = {ref[col], ref[col +1],
                              query[row -1], query[row +1]};
        continue_from_bind.delG = previous_decision_record.delG \
                                  + bulge_score(continue_from_bind.bottom_loop_len)
                                  + get_delG_internal(nn_config);
    } else if (previous_decision_record.current_decision == MISMATCH)
    {
        if (previous_decision_record.top_loop_len == 1 && previous_decision_record.bottom_loop_len == 1)
        {// need to backtrack previous delG addition
            continue_from_bind.top_loop_len = 1; // the one from previous mismatch
            continue_from_bind.bottom_loop_len = 2; // add one for current bulge

            /* neighbour:
             *     MX
             *     MXB   where B is the current bulge 
             *     previous binding was treated as a mismatch, 
             *     but now it's part of a loop, the mismatch contribution
             *     needs to be backtracked. */
            Neighbour nn_config = {ref[col -1], ref[col],
                                  query[row -2], query[row -1]};
            continue_from_bind.delG = previous_decision_record.delG \
                                      - get_delG_internal(nn_config) \
                                      + internal_loop_score(continue_from_bind.top_loop_len,
                                                            continue_from_bind.bottom_loop_len);
        } else if (previous_decision_record.top_loop_len > 1 || previous_decision_record.bottom_loop_len > 1)
        {
            continue_from_bind.bottom_loop_len = previous_decision_record.bottom_loop_len +1;
            continue_from_bind.top_loop_len = previous_decision_record.top_loop_len;
            continue_from_bind.delG = previous_decision_record.delG \
                                      + internal_loop_score(continue_from_bind.top_loop_len,
                                                            continue_from_bind.bottom_loop_len);
        }
    }
    // now handle continue from previous bottom_bulge: 2 cases
    // previous bulge has size 1: need to backtrack the special size one intervening delG addition
    // previous bulge size > 1: simply extend bulge size.
    previous_decision_record = prev_entry.bottom_bulge;
    if (previous_decision_record.bottom_loop_len == 1 && previous_decision_record.top_loop_len == 0)
    {
        /* neighbour configuration:
         *     m M  the m at the top is current col
         *     Mbm  where b is previous bulge.
         *     "m" is current position (row) which needs to be changed to bulge. 
         *     the addition of mM/Mm delG needs to be backtracked. */
        Neighbour nn_config = {ref[col], ref[col +1],
                               query[row -2], query[row]};
        continue_from_bottom_bulge.bottom_loop_len = 2; // added 1 to prev len
        continue_from_bottom_bulge.delG = previous_decision_record.delG \
                                          - get_delG_internal(nn_config) \
                                          + internal_loop_score(continue_from_bottom_bulge.top_loop_len,
                                                                continue_from_bottom_bulge.bottom_loop_len);
    } else if (previous_decision_record.bottom_loop_len > 1 || previous_decision_record.top_loop_len > 0)
    {
        continue_from_bottom_bulge.bottom_loop_len = previous_decision_record.bottom_loop_len +1;
        continue_from_bottom_bulge.top_loop_len = previous_decision_record.top_loop_len;
        continue_from_bottom_bulge.delG = previous_decision_record.delG \
                                          + internal_loop_score(continue_from_bottom_bulge.top_loop_len,
                                                                continue_from_bottom_bulge.bottom_loop_len);
    }

    // now both records are ready, return the better one.
    Decision_Record two_continuation_records[] = {continue_from_bind,
                                                continue_from_bottom_bulge};
    return best_record(two_continuation_records, 2);

}

/* score_stop:
 * 
 */
Decision_Record score_stop(SW_Entry **sw_matrix, 
                           int row, int col,
                           char *ref, char *query)
{
    SW_Entry prev_entry = sw_matrix[row][col];
    Decision_Record continue_from_bind = {0, MATCH, STOP, 0, 0};
    Decision_Record continue_from_top_bulge = {0, TOP_BULGE, STOP, 0, 0};
    Decision_Record continue_from_bottom_bulge = {0, BOTTOM_BULGE, STOP, 0, 0};
    
    // handle continue from previous match or mismatch,
    // 


    Decision_Record two_continuation_records[] = {continue_from_bind,
                                                continue_from_bottom_bulge};
    return best_record(two_continuation_records, 2);

}


/********************** THERMODYNAMICS ROUTINES ****************************/

float internal_loop_score(int top_loop_len, int bottom_loop_len)
{
    return 0.0;
}

float bulge_score(int loop_len)
{
    return 0.0;
}



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
        case '.':
             digit = TERMINAL_DOT;
             break;
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

float get_delG_terminal(Neighbour nn_config)
{
    int index = _get_index_terminal(nn_config);
    Therm_Param record = nn_data_terminal[index];
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
    Decision_Record null_decision = {0.0, STOP, STOP, 0, 0};
    SW_Entry null_entry = {
                           null_decision,
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




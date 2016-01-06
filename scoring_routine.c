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

#include <stdlib.h>
#include <stdio.h>
#include "swnn2.h"




/* score_bind:
 * evaluate the score (delG) of deciding on matching (or mismatching) current
 * bases as continuation from [bind, bottom_bulge, top_bulge]. 
 * Return the record of best continuation. 
 *
 * As mentioned, a binding condition can be a continuation from previous
 * bind, top_bulge or bottom_bulge.
 * There're multiple subcases in all 3 possibilities.
 * The function check all possibilities and return the best (lowest delG)
 * continuation.
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
 * compute the best (lowest delG) continuation if we are to choose
 * "top_bulge" as "current_decicision". As mentioned, a top_bulge
 * can be continuation from previous bind or top_bulge. Both possibilities
 * are evaluated and the best one is returned.
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
 * compute the best (lowest delG) continuation if we are to choose
 * "bottom_bulge" as "current_decicision". As mentioned, a bottom_bulge
 * can be continuation from previous bind or bottom_bulge. Both possibilities
 * are evaluated and the best one is returned.
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
 * compute the score if we are to stop the hybridisation at this point.
 * The stop can be continuation from previous bind, top_bulge or bottom_bulge
 * The best continuation is returned.
 * The appropriate "previous_entry" is the diagonal (towards top left).
 * Since the each decision made there can admit a stop continuation.
 *     a match / mismatch can be backtracked and 
 *     replace by terminal match/mismatch.
 *     a bulge (top or bottom) already assume that the next
 *     pair (a pair is a diagonal movement, where both index are incremented)
 *     of base is binding thus allowing a stop.
 *
 *
 * !!! see "!!!" comments in the function. If this function is intented for
 * !!! internal bases and without exhausting either one of the sequence,
 * !!! the situation doesn't count as terminal or dangling end, then
 * !!! this function is entirely useless. The only relevant stop option
 * !!! will be at the last row and last column.
 */
Decision_Record score_stop(SW_Entry **sw_matrix, 
                           int row, int col,
                           char *ref, char *query)
{
    SW_Entry prev_entry = sw_matrix[row -1][col -1];
    Decision_Record previous_decision_record;
    Decision_Record continue_from_bind = {0, MATCH, STOP, 0, 0};
    Decision_Record continue_from_top_bulge = {0, TOP_BULGE, STOP, 0, 0};
    Decision_Record continue_from_bottom_bulge = {0, BOTTOM_BULGE, STOP, 0, 0};
    
    // handle continue from previous match or mismatch: 2 cases
    // !!! REMAIN TO BE ANSWERED: which of the below to choose?
    // !!! Don't know if terminate before neither sequence is
    // !!! exhausted count as terminal or dangling end
    // Scheme 1: count as terminal
    // previous binding is a match: add in terminal AT penalty
    // previous binding is a mismatch: 2 subcases
    //     -> loop len 1: backtrack the mismatch and replace it 
    //                    with terminal mismatch data
    //     -> loop len > 1: it's an internal loop, already assume current
    //                      base is a match, thus do nothing
    //                      (!!! add one single binding energy????)
    // then, whichever the case, we add in both 3' and 5' dangling end
    //
    // Scheme 2: not counted as terminal.
    // In all cases, do nothing. (we will use this first ... )
    continue_from_bind.previous_decision = prev_entry.bind.current_decision;
    continue_from_bind.delG = prev_entry.bind.delG;

    // handle continue from previous top_bulge.
    // !!! SAME SITUATION AS ABOVE.
    continue_from_top_bulge.previous_decision = TOP_BULGE;
    continue_from_top_bulge.delG = prev_entry.top_bulge.delG;
    
    // handle continue from previus bottom_bugle.
    // !!! SAME AS ABOVE
    continue_from_bottom_bulge.previous_decision = BOTTOM_BULGE;
    continue_from_bottom_bulge.delG = prev_entry.bottom_bulge.delG;

    Decision_Record three_continuation_records[] = {continue_from_bind,
                                                    continue_from_top_bulge,
                                                    continue_from_bottom_bulge};
    return best_record(three_continuation_records, 3);
}


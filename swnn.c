/* Implementation of secondary structure prediction algorithm using
 * idea from Smith-Waterman Alignement algorithm and Nearest neighbour
 * thermodynamics. 
 */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swnn.h"


void _test_get_delG();
void _test();

int main(int argc, char **argv)
{
    char *ref = argv[1];
    char *query = argv[2];
    _test(ref, query);
    return EXIT_SUCCESS;
}

void _test(char *ref, char *query)
{
    SW_Entry **sw_matrix = complete_duplex_matrix(ref, query);
    Coord best_coord = find_best_decision_coord(sw_matrix, strlen(query), strlen(ref));
    Decision_Record best_decision;
    SW_Entry best_entry = sw_matrix[best_coord.row][best_coord.col];
    Decision_Record options[] = {best_entry.bind, best_entry.top_bulge, best_entry.bottom_bulge};
    best_decision = best_record(options , 3);

    printf("delG in best decision = %f\n decision = %c at row=%i, col=%i\n", 
           best_decision.delG, best_decision.current_decision, best_coord.row, best_coord.col);

    print_duplex(sw_matrix, best_coord, ref, query);
}




/****************************************************************************
 * DP ROUTINES 
 * *************************************************************************/

/* complete_duplex_matrix:
 * Given a reference sequence (sense sequence)
 * and a query sequence (antisense sequence)
 * return a DP matrix recording all the decisions and scores. 
 */
SW_Entry **complete_duplex_matrix(char *ref, char *query)
{
    // initialise matrix base on 
    // the matrix layout:
    // reference is on the horizontal while
    // query is at the bottom. 
    // Initiation is considerate of dangling ends and
    // init_AT or init_GC scenarios.
    int nrow = strlen(query);
    int ncol = strlen(ref);
    SW_Entry **sw_matrix = initialise_duplex_matrix(ref, query);
    // now we fill up the matrix
    register int row, col;
    // start from 1, since the 0th row and col 
    // had been filled during initiation
    for (row = 1; row < nrow -1; row++)
    {
        for (col = 1; col < ncol -1; col++)
        {
            sw_matrix[row][col] = compute_internal_entry(sw_matrix,
                                                         row, col,
                                                         ref, query);
        }
    }
    sw_matrix = process_last_row_col(sw_matrix, ref, query);
    return sw_matrix;
}


SW_Entry compute_internal_entry(SW_Entry **sw_matrix,
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
    // entry.stop = score_stop(sw_matrix,
    //                        row, col,
    //                        ref, query);
    return entry;
}


/* find_best_decision: find the best recorded decision in the whole
 * sw_matrix. We look through only the last row and column.
 */
Coord find_best_decision_coord(SW_Entry **sw_matrix, int nrow, int ncol)
{
    register int row, col;
    Coord best_coord;
    SW_Entry this_entry;
    Decision_Record this_record;
    float lowest_delG = 0.0;
    const int num_choice = 3;
    int choice;
    // look through the last column then the last row
    // we favour the bottom right entries.
    // first, the last column
    for (col = ncol -1, row = 0; row < nrow -1; row++)
    {
        this_entry = sw_matrix[row][col];
        Decision_Record all_records[] = {this_entry.bind,
                                       this_entry.top_bulge,
                                       this_entry.bottom_bulge};
        for (choice = 0; choice < num_choice; choice++)
        {
            this_record = all_records[choice];
            if (this_record.delG < lowest_delG)
            {
                lowest_delG = this_record.delG;
                best_coord = (Coord) {row, col, this_record.current_decision};
            }
        }
    }
    // then the last row
    for (col = 0, row = nrow -1; col < ncol; col++) 
    {
        this_entry = sw_matrix[row][col];
        Decision_Record all_records[] = {this_entry.bind,
                                       this_entry.top_bulge,
                                       this_entry.bottom_bulge};
        for (choice = 0; choice < num_choice; choice++)
        {
            this_record = all_records[choice];
            if (this_record.delG < lowest_delG)
            {
                lowest_delG = this_record.delG;
                best_coord = (Coord) {row, col, this_record.current_decision};
            }
        }
    }
    return best_coord;
}


    
    
Coord find_best_entry_coord(SW_Entry **sw_matrix, int nrow, int ncol)
{
    register int row, col;
    float lowest_delG = 0.0;
    float new_delG;
    Coord best_coord = {0, 0, MATCH};
    SW_Entry current_entry;
    Decision_Record record;
    for (row = 0; row < nrow; row++)
    {
        for (col = 0; col < ncol; col++)
        {
            current_entry = sw_matrix[row][col];
            Decision_Record three_options[] = {current_entry.bind,
                                               current_entry.top_bulge,
                                               current_entry.bottom_bulge};
            record = best_record(three_options, 3);
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




/************************** UTILITIES ROUTINES ******************************/

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
Decision_Record best_record(Decision_Record records[], int nrecord)
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

void print_duplex(SW_Entry **sw_matrix, Coord coord, char *ref, char *query)
{
    int ref_len = strlen(ref);
    int query_len = strlen(query);
    int repr_len = ref_len + query_len; //fmax(ref_len, query_len);
    char ref_string[repr_len +1];
    char bindings_string[repr_len +1];
    char query_string[repr_len +1];
    ref_string[repr_len] = bindings_string[repr_len] = query_string[repr_len] = '\0';

    Decision_Record current_decision_record;
    // back tracing
    int row = coord.row;
    int col = coord.col;
    char decision = coord.current_decision;
    current_decision_record = get_decision_from_entry(sw_matrix[row][col],
                                                      decision);
    register int i;
    for (i = repr_len -1; i >= 0; i--)
    {
        if (i <= row)
        {
            if (i <= col)
            {
                decision = current_decision_record.current_decision;
                if (decision == MATCH || 
                    decision == MISMATCH)
                {
                    ref_string[i] = ref[col];
                    bindings_string[i] = (decision == MATCH) ? BOND : XBOND;
                    query_string[i] = query[row];
                    row -= 1;
                    col -= 1;
                } else if (decision == TOP_BULGE)
                {
                    ref_string[i] = ref[col];
                    bindings_string[i] = EMPTY_SPACE;
                    query_string[i] = BULGE_GAP;
                    col -= 1;
                } else if (decision == BOTTOM_BULGE)
                {
                    ref_string[i] = BULGE_GAP;
                    bindings_string[i] = EMPTY_SPACE;
                    query_string[i] = query[row];
                    row -= 1;
                }
                if (row > 0 && col > 0)
                {
                    current_decision_record = \
                        get_decision_from_entry(sw_matrix[row][col],
                                                current_decision_record.previous_decision);
                }
            } else if (i > col)
            {
                ref_string[i] = EMPTY_SPACE;
                bindings_string[i] = EMPTY_SPACE;
                query_string[i] = query[i];
            }
        } else if (i > row)
        {
            bindings_string[i] = EMPTY_SPACE;
            if (i <= col)
            {
                ref_string[i] = ref[i];
                query_string[i] = EMPTY_SPACE;
            } else if (i > col)
            {
                ref_string[i] = EMPTY_SPACE;
                query_string[i] = EMPTY_SPACE;
            }
        }
    }
    printf("%s\n%s\n%s\n", ref_string, bindings_string, query_string);
}




Decision_Record get_decision_from_entry(SW_Entry entry, char decision)
{
    Decision_Record current_decision;
    switch (decision)
    {
        case (MATCH): case (MISMATCH):
            current_decision = entry.bind;
            break;
        case (TOP_BULGE):
            current_decision = entry.top_bulge;
            break;
        case (BOTTOM_BULGE):
            current_decision = entry.bottom_bulge;
            break;
        default:
            {
                fprintf(stderr, 
                    "swnn get_decision_from_entry: "
                    "illegal decision character: %c\n", decision);
                exit(EXIT_FAILURE);
            }
            break;
    }
    return current_decision;
}



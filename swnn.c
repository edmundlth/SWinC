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
Decision_Record best_record(Decision_Record records[], int nrecord);

int main()
{
    Neighbour nn_config = {'A', 'G', 'T','C'};
    extern const Therm_Param GLOBAL_nn_data_internal[];
    extern float GLOBAL_Reaction_Temperature;
    Therm_Param from_record = GLOBAL_nn_data_internal[_get_index_internal(nn_config)];
    printf("%s %f %f\n", from_record.neighbour, from_record.delH, from_record.delS);
    printf("delG = %f\n", from_record.delH * 1000.0 - (GLOBAL_Reaction_Temperature + ABSOLUTE_ZERO_OFFSET) * from_record.delS);
    printf("delG_function = %f\n", get_delG_internal(nn_config));
    return 0;
}


/************************** ALIGNMENT ROUTINES ******************************/

/* duplex_matrix:
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


/* find_best_decision: find the best recorded decision in the whole
 * sw_matrix. We look through only the last row and column.
 */
Coord find_best_decision(SW_Entry **sw_matrix, int nrow, int ncol)
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




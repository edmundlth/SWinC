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




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

#define MAX_SEQ_LEN 60
#define MAX_POOL_SIZE 3000

void _test(char **pool, int pool_size, char *outfilename);
void _test_pool(char **pool, int pool_size);
void _test_pool_complement(char **pool, int pool_size);
int extract_pool(char *filename, char *buffer[]);
char *reverse(char *string, char buffer[]);
int main(int argc, char **argv)
{
    char **pool = malloc(sizeof(char *) * MAX_POOL_SIZE);
    int i;
    for (i = 0; i < MAX_POOL_SIZE; i++)
    {
        pool[i] = malloc(sizeof(char) * (MAX_SEQ_LEN +1));
    }
    int num_seq = extract_pool(argv[1], pool);
    _test(pool, num_seq, argv[2]);
    return EXIT_SUCCESS;
}

int extract_pool(char *filename, char *buffer[])
{
    FILE *inputfile = fopen(filename, "r");
    int i;
    for (i = 0; i < MAX_POOL_SIZE; i++)
    {
        buffer[i] = fgets(buffer[i], (MAX_SEQ_LEN +1), inputfile);
        if (buffer[i] != NULL)
        {
            buffer[i][strlen(buffer[i]) -2] = '\0'; // remove newline
        }
        else
        {
            break;
        }
    }
    return i;
}


void _test(char **pool, int pool_size, char *outfilename)
{
    FILE *outfile = fopen(outfilename, "r");
    float best_delG, new_delG;
    char *ref;
    char query[MAX_SEQ_LEN];
    char *best_partner;
    Coord best_coord;
    int ref_len, query_len;
    int i, j;

    for (i = 0; i < pool_size; i++)
    {
        best_delG = 0.0;
        ref = pool[i];
        ref_len = strlen(ref);
        best_partner = "-";
        for (j = 0; j < pool_size; j++)
        {
            if (i != j)
            {
                reverse(pool[j], query);
                query_len = strlen(query);
                SW_Entry **sw_matrix = complete_duplex_matrix(ref, query);
                best_coord = find_best_decision_coord(sw_matrix, ref_len, query_len);
                new_delG = get_decision_from_entry(
                              sw_matrix[best_coord.row][best_coord.col],
                              best_coord.current_decision).delG;
                if (new_delG < best_delG)
                {
                    best_delG = new_delG;
                    best_partner = query;
                }
            }
        }
        fprintf(outfile, "%s\t%s\t%f", ref, best_partner, best_delG);
    }
    fclose(outfile);
}

char *reverse(char *string, char *buffer)
{
    int i;
    int length = strlen(string);
    for (i = length; i >= 1; i--)
    {
        buffer[i -1] = string[length -i];
    }
    buffer[length] = '\0';
    return buffer;
}

char *complement_seq(char *string, char buffer[])
{
    int length = strlen(string);
    int i;
    buffer[length] = '\0';
    for (i = 0; i < length; i++)
    {
        buffer[i] = complement(string[i]);
    }
    return buffer;
}

void _test_pool(char **pool, int pool_size)
{
    char ref[61];
    char query[61];
    SW_Entry **sw_matrix;
    Coord best_coord;
    int i, j;
    for (i = 0; i < pool_size; i++)
    {
        for (j = i; j < pool_size; j++)
        {
            strcpy(ref, pool[i]);
            reverse(pool[j], query);
            printf("(i, j) = (%i, %i)\n", i, j);
            printf("ref  = %s\nquery= %s\n", ref, query);
            sw_matrix = complete_duplex_matrix(ref, query);
            best_coord = find_best_decision_coord(sw_matrix, strlen(query), strlen(ref));
            print_duplex(sw_matrix, best_coord, ref, query);
        }
    }
}

void _test_pool_complement(char **pool, int pool_size)
{
    char *ref;
    char query[61];
    SW_Entry **sw_matrix;
    Coord best_coord;
    int i;
    for (i = 0; i < pool_size; i++)
    {
        ref = pool[i];
        complement_seq(pool[i], query);
        printf("i= %i\nref  = %s\nquery= %s\n", i, ref, query);
        sw_matrix = complete_duplex_matrix(ref, query);
        best_coord = find_best_decision_coord(sw_matrix, strlen(query), strlen(ref));
        print_duplex(sw_matrix, best_coord, ref, query);
    }
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
    Coord best_coord = {0, 0, MATCH};
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
    int start_row = coord.row;
    int start_col = coord.col;
    int ref_len = strlen(ref);
    int query_len = strlen(query);
    int repr_len = ref_len + query_len;
    char *ref_string = calloc(repr_len +1, sizeof(char));
    char *bindings_string = calloc(repr_len +1, sizeof(char));
    char *query_string = calloc(repr_len +1, sizeof(char));
    int j;
    for (j = 0; j < repr_len; j++)
    {
        ref_string[j] = bindings_string[j] = query_string[j] = ' ';
    }


    printf("delG(duplex) = %f\n", get_decision_from_entry(sw_matrix[start_row][start_col],
                                                        coord.current_decision).delG);
    // back tracing
    int row = start_row;
    int col = start_col;
    int i;
    char decision = coord.current_decision;
    Decision_Record current_decision_record;
    repr_len = ref_len + query_len;
    for (i = repr_len -1; i >= 0; i--)
    {
        if (row >= 0 && col >= 0)
        {
            current_decision_record = \
                    get_decision_from_entry(sw_matrix[row][col],
                                            decision);
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
            decision = current_decision_record.previous_decision;
        } else if (row >= 0 && col < 0)
        {
            ref_string[i] = EMPTY_SPACE;
            bindings_string[i] = EMPTY_SPACE;
            query_string[i] = query[row];
            row--;
        } else if (row < 0 && col >= 0)
        {
            ref_string[i] = ref[col];
            bindings_string[i] = EMPTY_SPACE;
            query_string[i] = EMPTY_SPACE;
            col--;
        } else
        {
            ref_string[i] = bindings_string[i] = query_string[i] = EMPTY_SPACE;
        }
//        printf("i=%i\nref  =%s\nbind =%s\nquery=%s\n\n", i, ref_string, bindings_string, query_string);
    }
    printf("5'-%s-3'\n   %s   \n3'-%s-5'\n\n", ref_string, bindings_string, query_string);
    free(ref_string);
    free(bindings_string);
    free(query_string);
}



Decision_Record get_decision_from_entry(SW_Entry entry, char decision)
{
    Decision_Record current_decision;
    switch (decision)
    {
        case (MATCH): 
            current_decision = entry.bind;
            break;
        case (MISMATCH):
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



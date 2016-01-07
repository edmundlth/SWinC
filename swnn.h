/* Header file for SWNN algorithm 
 * !! Put summary of the various 
 * !! names, variables and routine defined here.
 */
#define TRUE 1
#define FALSE 0
#define MATCH 'M'
#define MISMATCH 'X'
#define TOP_BULGE 'T'
#define BOTTOM_BULGE 'B'
#define STOP 'S'

#define ABSOLUTE_ZERO_OFFSET 273.15

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
/* type that record the score
 * for each possible decision 
 * Recording current state is redundant since
 * given an initial state, the path can be traced
 * by only information about the next state. 
 * previous decision of the node where the current
 * node is pointed from should always equal to the
 * current state.
 * But it improve clarity of the code. */
typedef struct 
{
    float delG;
    char previous_decision; // how did we get here?
    char current_decision; // where are we now?
    int top_loop_len; // for match, loop_len should be 0
    int bottom_loop_len;
} Decision_Record;

/* The structure that holds information about all 3 different
 * possible decisions at each  sw_matrix entry */
typedef struct
{
    Decision_Record bind;
    Decision_Record top_bulge;
    Decision_Record bottom_bulge;
    Decision_Record stop;
} SW_Entry;

/* The structure that represent the coordinate
 * of a particular decision in the sw_matrix.
 * It records the entry coordinate (row, col)
 * and the decision made in in the entry, which 
 * will be one of {'M', 'X', 'T', 'B', 'S'}
 * 'M' -> Match
 * 'X' -> Mismatch
 * 'T' -> Top_Bulge
 * 'B' -> Bottom_Bulge
 * 'S' -> Stop / Termination (terminate alignment here) */
typedef struct
{
    int row;
    int col;
    char current_decision;
} Coord;


/* The structure that represent the a NN pairing. 
 * e.g. Given:
 * top5 = 'A', top3 = 'G', bottom3 = 'T', bottom5 = 'C'
 * We have the following binding:
 * 5'-AG-3'
 *    ||
 * 3'-TC-5'
 */
typedef struct
{
    char top5;
    char top3;
    char bottom3;
    char bottom5;
} Neighbour;

/* The structure that record NN_data.
 * e.g. Therm_Param record = {"AG/TC", -7.8, -21.0}
 * is the record of enthalpy and entropy value of
 * the neighbouring base paring 
 *     5'-AG-3'
 *     3'-TC-5'
 */
typedef struct {
    char *neighbour;
    float delH;
    float delS;
} Therm_Param;

/* This structure record the difference of entropy
 * going from loop_size -1 to loop_size for both
 * internal loop and bulge loop.
 * For any loop, enthalpy values are always assumed
 * to be 0, hence the change in enthalpy is also 0 */
typedef struct {
    int loop_size;
    float deldelS;
} Loop_Entropy_Diff;

/************************** ALIGNMENT ROUTINES ******************************/
SW_Entry **complete_duplex_matrix(char *ref, char *query);
SW_Entry **initialise_duplex_matrix(char *ref, char *query);
SW_Entry compute_entry(SW_Entry **sw_matrix, 
                       int row, int col, 
                       char *ref, char *query);
Coord find_best_entry_coord(SW_Entry **sw_matrix, int nrow, int ncol);

/************************** SCORING ROUTINES ******************************/
Decision_Record score_bind(SW_Entry **sw_matrix,
                            int row, int col,
                            char *ref, char *query);
Decision_Record score_top_bulge(SW_Entry **sw_matrix,
                                int row, int col,
                                char *ref, char *query);
Decision_Record score_bottom_bulge(SW_Entry **sw_matrix,
                                   int row, int col,
                                   char *ref, char *query);
Decision_Record score_stop(SW_Entry **sw_matrix,
                           int row, int col,
                           char *ref, char *query);

/********************** THERMODYNAMICS ROUTINES ****************************/
float internal_loop_score(int top_loop_len, int bottom_loop_len);
float bulge_score(int loop_len);
int _get_index_internal(Neighbour nn_config);
int _get_index_terminal(Neighbour nn_config);
float get_delG_internal(Neighbour nn_config);
float get_delG_terminal(Neighbour nn_config);
float init_delG(char base);

SW_Entry **_allocate_matrix(int nrow, int ncol);
SW_Entry **initialise_duple_matrix(char *ref, char *qeury);
SW_Entry _handle_first_entry(char first_ref, char first_query);
SW_Entry _handle_init_row_col(Neighbour nn_config);
SW_Entry **process_last_row_col(SW_Entry **sw_matrix, 
                                char *ref, char *query);
SW_Entry compute_last_entry(SW_Entry **sw_matrix,
                            int row, int col,
                            char *ref, char *query);
Decision_Record score_bind_terminal(SW_Entry **sw_matrix,
                                    int row, int col,
                                    char *ref, char *query);
Decision_Record score_top_bulge_terminal(SW_Entry **sw_matrix,
                                         int row, int col,
                                         char *ref, char *query);
Decision_Record score_bottom_bulge_terminal(SW_Entry **sw_matrix,
                                            int row, int col,
                                            char *ref, char *query);
float _get_dangling_end_delG(char *ref, char *query, int row, int col);
/************************** UTILITIES ROUTINES ******************************/
char complement(char base);
int is_complement(char base1, char base2);
Decision_Record best_record(Decision_Record records[], int nrecord);

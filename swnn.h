/* Header file for SWNN algorithm 
 * !! Put summary of the various 
 * !! names, variables and routine defined here.
 */
#define TRUE 1
#define FALSE 0
#define MATCH 'M'
#define MISMATCH 'X'
#define INSERTION 'I'
#define DELETION 'D'
#define STOP 'S'

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
    char current_state; // where are we now?
    int loop_len; // for match, loop_len should be 0
} Decision_Record

/* The structure that holds information about all 3 different
 * possible decisions at each  sw_matrix entry */
typedef struct
{
    Decision_Record match;
    Decision_Record insertion
    Decision_Record deletion;
} SW_Entry

/* The structure that represend the coordinate
 * of a particular decision in the sw_matrix.
 * It records the entry coordinate (row, col)
 * and the decision made in in the entry, which 
 * will be one of {'M', 'X', 'D', 'I', 'S'}
 * 'M' -> Match
 * 'X' -> Mismatch
 * 'D' -> Deletion
 * 'I' -> Insertion
 * 'S' -> Stop / Termination (terminate alignment here) */
typedef struct
{
    int row;
    int col;
    char current_decision;
} Coord


typedef struct
{
    char top5;
    char top3;
    char bottom3;
    char bottom5;
} Neighbour


/************************** ALIGNMENT ROUTINES ******************************/
float swnnalign(char *ref, char *query);
SW_Entry compute_entry(SW_Entry *sw_matrix, 
                       int row, int col, 
                       char *ref, char *query);
float find_best_score(SW_Entry **sw_matrix);

/************************** SCORING ROUTINES ******************************/
Decision_Record score_match(SW_Entry **sw_matrix,
                            int row, int col,
                            char *ref, char *query);
Decision_Record score_insertion(SW_Entry **sw_matrix,
                                int row, int col,
                                char *ref, char *query);
Decision_Record score_deletion(SW_Entry **sw_matrix,
                               int row, int col,
                               char *ref, char *query);

/********************** THERMODYNAMICS ROUTINES ****************************/
float get_delG(Neighbour nn_config);
float size_1_bulge(Neighbour intervening_bases_config);
float extend_internal_loop(int previous_loop_len);
float extend_bulge_loop(int previous_loop_len);


/************************** UTILITIES ROUTINES ******************************/
char complement(char base);
int is_complement(char base1, char base2);


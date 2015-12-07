
/* define maximum primer + heel length acceptable */
#define MAX_SEQ_LEN 60
#define MAX_POOL_SIZE 100

/* produce record of nearest neighbour thermodynamic
 * data here. It should be in "hash table" like 
 * structure
 */



/***** Parameters and variables for SW alignment ************/
float match_score = 1.0;
float mismatch_penalty = -1.5;
float gap_open = -1.5;
float gap_extension = 0.5;

/* Record score and decision at every sw_matrix entry
 * the decision are 
 *   'M' for match
 *   'm' for mismatch
 *   'D' for deletion
 *   'I' for insertion
 *   '\0' for termimate alignment */
typedef struct {
    float score;
    char decision;
} SW_entry;
static SW_entry null_entry = {0.0, '\0'};
/* initialise sw_matrix as (MAX_SEQ_LEN+1) x (MAX_SEQ_LEN+1)
 * all entries will be initialised with SW_entry = {0.0,'\0'} */
static SW_entry sw_matrix[MAX_SEQ_LEN +1][MAX_SEQ_LEN +1];
/* external variable recording the reference and query sequence
 * visible to all the alignment/scoring related functions */
char *g_ref;
char *g_query;


/***** Routines for sw alignment *************/
float swalign(char *ref, char *query);
float fill_matrix(int ref_len, int query_len);
float score(int row, int col);
SW_entry score_mm(int row, int col);
SW_entry score_insert(int row, int col);
SW_entry score_delete(int row, int col);
float penalise_gap(int gap_len);
SW_entry max(SW_entry list[], int list_len);



/******* Variables for pool alignment *****/
float interaction_matrix[MAX_POOL_SIZE][MAX_POOL_SIZE];
char pool[MAX_POOL_SIZE][MAX_SEQ_LEN];
/******* Routines for pool alignment ******/
void align_pool(int pool_size);
int get_primers(char *filename);





/**** Utilities Routines *****/
void print_sw_matrix(int nrow, int ncol);
void print_interaction_matrix(int nrow, int ncol);
char *rev_complement(char *seq);
float mean(float num_list[], int list_len);
char *trim_whitespace(char *input);


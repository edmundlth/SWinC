
/* maximum primer + heel length acceptable 
 * kmer size is the number of bases from the 
 * 3' end of query to be align to a reference*/
#define MAX_SEQ_LEN 60
#define MAX_POOL_SIZE 10000
#define KMER_SIZE 20
#define TEST_REF AAAATTTTGGGGCCCC
#define TEST_QUERY AAAATTAAAAGGGGCCCC

/***** Parameters and variables for SW alignment ************/
#define DEFAULT_MATCH_SCORE 1.0
#define DEFAULT_MISMATCH_PENALTY -1.0
#define DEFAULT_GAP_OPEN_PENALTY -2.0
#define DEFAULT_GAP_EXTENSION_PENALTY DEFAULT_GAP_OPEN_PENALTY
// if not specified, extending a gap cost as much as opening a gap

/* a record of the parameters that will be used in scoring */
typedef struct
{
    float match_score;
    float mismatch_penalty;
    float gap_open_penalty;
    float gap_extension_penalty;
} Score_Param;

/* a record of user commandline input */
typedef struct {
    char *ref;
    char *query;
    char *primer_filename;
    Score_Param score_param;
    int verbose_flag;
} User_Inputs;


/* Record score and decision at every sw_matrix entry
 * the decisions are 
 *   'M' for match
 *   'X' for mismatch
 *   'D' for deletion
 *   'I' for insertion
 *   'T' for termimate alignment 
 *   There will be 2 characters in each decision entry
 *   the first records the immediate previous decision 
 *   that leads to the current decision and the second
 *   record the current state.
 *   eg: {2.3, "DM"} says: previous alignment was a deletion
 *   and we follow on that with a match the score will
 *   be 2.3. */
typedef struct {
    float score;
    char *decision;
} Decision_Record;


/* Each entry in sw_matrix has 3 "points". 
 * Each of them represent a possible current state, which
 * are "in a match-mismatch", "in an insert state",
 * "in a deletion state". 
 * The later 2 states are instances of gapped state which
 * insert being gaps in the reference and deletion being
 * gaps in the query */
typedef struct {
    Decision_Record match_record;
    Decision_Record insert_record;
    Decision_Record delete_record;
} SW_entry;



/***** Routines for sw alignment *************/
float swalign(char *ref, char *query, Score_Param score_param);
float fill_matrix(SW_entry **sw_matrix, 
                  char *ref, char *query, 
                  Score_Param score_param);
float score(SW_entry **sw_matrix, 
            char *ref, char *query,
            Score_Param score_param);
Decision_Record score_match_mismatch(SW_entry **sw_matrix,
                              char *ref, char *query,
                              int row, int col,
                              Score_Param score_param);
Decision_Record score_insert(SW_entry **sw_matrix,
                      int row, int col,
                      Score_Param score_param);
Decision_Record score_delete(SW_entry **sw_matrix,
                      int row, int col,
                      Score_Param score_param);
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
void print_alignment(int nrow, int ncol);
char *prepend_char(char *string, char c);
int *best_entry(int nrow, int ncol);
void print_interaction_matrix(int nrow, int ncol);
char *rev_complement(char *seq, int result_len);
float mean(float num_list[], int list_len);
char *trim_whitespace(char *input);
User_Inputs parse_args(int argc, char **argv);


#define PROGRAM_NAME "swinc"
enum
{
    ERROR_MEM_ALLOC = 1
};

void error_handle(int error_code)
{
    fprintf(stderr, "%s error:\n", PROGRAM_NAME);
    switch (error_code)
    {
        case ERROR_MEM_ALLOC:
            fprintf(stderr, "Memory allocation error");
            break;
    }
}






/*************************************************
 * Nearest Neighbour Thermodynamics Parameters *
 * **********************************************/

#define ABSOLUTE_ZERO_OFFSET 273.15
float Reaction_Temperature;

typedef struct {
    char *neighbour;
    float dH_dS[2];
} Therm_Param;

Therm_Param Initialisation[] = {
    "init", {0, 0}, 
    "init_A/T", {2.3, 4.1}, 
    "init_G/C", {0.1, -2.8},
    "init_oneG/C", {0, 0}, 
    "init_allA/T", {0, 0}, 
    "init_5T/A", {0, 0},
    "sym", {0, -1.4},
};

Therm_Param Match[] = {
    "AA/TT", {-7.9, -22.2}, 
    "AT/TA", {-7.2, -20.4}, 
    "TA/AT", {-7.2, -21.3},
    "CA/GT", {-8.5, -22.7}, 
    "GT/CA", {-8.4, -22.4}, 
    "CT/GA", {-7.8, -21.0},
    "GA/CT", {-8.2, -22.2}, 
    "CG/GC", {-10.6, -27.2}, 
    "GC/CG", {-9.8, -24.4},
    "GG/CC", {-8.0, -19.9}
};

// Internal mismatch and inosine table {DNA}
// Allawi & SantaLucia {1997}, Biochemistry 36, 10581-10594
// Allawi & SantaLucia {1998}, Biochemistry 37, 9435-9444
// Allawi & SantaLucia {1998}, Biochemistry 37, 2170-2179
// Allawi & SantaLucia {1998}, Nucl Acids Res 26, 2694-2701
// Peyret et al. {1999}, Biochemistry 38, 3468-3477
// Watkins & SantaLucia {2005}, Nucl Acids Res 33, 6258-6267
Therm_Param Internal_Mismatch[] = {
    "AG/TT", {1.0, 0.9}, 
    "AT/TG", {-2.5, -8.3}, 
    "CG/GT", {-4.1, -11.7},
    "CT/GG", {-2.8, -8.0}, 
    "GG/CT", {3.3, 10.4}, 
    "GG/TT", {5.8, 16.3},
    "GT/CG", {-4.4, -12.3}, 
    "GT/TG", {4.1, 9.5}, 
    "TG/AT", {-0.1, -1.7},
    "TG/GT", {-1.4, -6.2},
    "TT/AG", {-1.3, -5.3},
    "AA/TG", {-0.6, -2.3},
    "AG/TA", {-0.7, -2.3},
    "CA/GG", {-0.7, -2.3},
    "CG/GA", {-4.0, -13.2},
    "GA/CG", {-0.6, -1.0},
    "GG/CA", {0.5, 3.2},
    "TA/AG", {0.7, 0.7},
    "TG/AA", {3.0, 7.4},
    "AC/TT", {0.7, 0.2},
    "AT/TC", {-1.2, -6.2},
    "CC/GT", {-0.8, -4.5},
    "CT/GC", {-1.5, -6.1},
    "GC/CT", {2.3, 5.4},
    "GT/CC", {5.2, 13.5},
    "TC/AT", {1.2, 0.7},
    "TT/AC", {1.0, 0.7},
    "AA/TC", {2.3, 4.6}, 
    "AC/TA", {5.3, 14.6}, 
    "CA/GC", {1.9, 3.7},
    "CC/GA", {0.6, -0.6}, 
    "GA/CC", {5.2, 14.2}, 
    "GC/CA", {-0.7, -3.8},
    "TA/AC", {3.4, 8.0}, 
    "TC/AA", {7.6, 20.2},
    "AA/TA", {1.2, 1.7}, 
    "CA/GA", {-0.9, -4.2}, 
    "GA/CA", {-2.9, -9.8},
    "TA/AA", {4.7, 12.9}, 
    "AC/TC", {0.0, -4.4}, 
    "CC/GC", {-1.5, -7.2},
    "GC/CC", {3.6, 8.9}, 
    "TC/AC", {6.1, 16.4}, 
    "AG/TG", {-3.1, -9.5},
    "CG/GG", {-4.9, -15.3}, 
    "GG/CG", {-6.0, -15.8}, 
    "TG/AG", {1.6, 3.6},
    "AT/TT", {-2.7, -10.8}, 
    "CT/GT", {-5.0, -15.8}, 
    "GT/CT", {-2.2, -8.4},
    "TT/AT", {0.2, -1.5},
    "AI/TC", {-8.9, -25.5}, 
    "TI/AC", {-5.9, -17.4}, 
    "AC/TI", {-8.8, -25.4},
    "TC/AI", {-4.9, -13.9}, 
    "CI/GC", {-5.4, -13.7}, 
    "GI/CC", {-6.8, -19.1},
    "CC/GI", {-8.3, -23.8}, 
    "GC/CI", {-5.0, -12.6},
    "AI/TA", {-8.3, -25.0}, 
    "TI/AA", {-3.4, -11.2}, 
    "AA/TI", {-0.7, -2.6},
    "TA/AI", {-1.3, -4.6}, 
    "CI/GA", {2.6, 8.9}, 
    "GI/CA", {-7.8, -21.1},
    "CA/GI", {-7.0, -20.0}, 
    "GA/CI", {-7.6, -20.2},
    "AI/TT", {0.49, -0.7}, 
    "TI/AT", {-6.5, -22.0}, 
    "AT/TI", {-5.6, -18.7},
    "TT/AI", {-0.8, -4.3}, 
    "CI/GT", {-1.0, -2.4}, 
    "GI/CT", {-3.5, -10.6},
    "CT/GI", {0.1, -1.0}, 
    "GT/CI", {-4.3, -12.1},
    "AI/TG", {-4.9, -15.8}, 
    "TI/AG", {-1.9, -8.5}, 
    "AG/TI", {0.1, -1.8},
    "TG/AI", {1.0, 1.0}, 
    "CI/GG", {7.1, 21.3}, 
    "GI/CG", {-1.1, -3.2},
    "CG/GI", {5.8, 16.9}, 
    "GG/CI", {-7.6, -22.0},
    "AI/TI", {-3.3, -11.9}, 
    "TI/AI", {0.1, -2.3}, 
    "CI/GI", {1.3, 3.0},
    "GI/CI", {-0.5, -1.3}
};

// Terminal mismatch table (DNA)
// SantaLucia & Peyret (2001) Patent Application WO 01/94611
Therm_Param Terminal_Mismatch[] = {
    "AA/TA", {-3.1, -7.8}, 
    "TA/AA", {-2.5, -6.3}, 
    "CA/GA", {-4.3, -10.7},
    "GA/CA", {-8.0, -22.5},
    "AC/TC", {-0.1, 0.5}, 
    "TC/AC", {-0.7, -1.3}, 
    "CC/GC", {-2.1, -5.1},
    "GC/CC", {-3.9, -10.6},
    "AG/TG", {-1.1, -2.1}, 
    "TG/AG", {-1.1, -2.7}, 
    "CG/GG", {-3.8, -9.5},
    "GG/CG", {-0.7, -19.2},
    "AT/TT", {-2.4, -6.5}, 
    "TT/AT", {-3.2, -8.9}, 
    "CT/GT", {-6.1, -16.9},
    "GT/CT", {-7.4, -21.2},
    "AA/TC", {-1.6, -4.0}, 
    "AC/TA", {-1.8, -3.8}, 
    "CA/GC", {-2.6, -5.9},
    "CC/GA", {-2.7, -6.0}, 
    "GA/CC", {-5.0, -13.8},
    "GC/CA", {-3.2, -7.1},
    "TA/AC", {-2.3, -5.9},
    "TC/AA", {-2.7, -7.0},
    "AC/TT", {-0.9, -1.7},
    "AT/TC", {-2.3, -6.3}, 
    "CC/GT", {-3.2, -8.0},
    "CT/GC", {-3.9, -10.6}, 
    "GC/CT", {-4.9, -13.5}, 
    "GT/CC", {-3.0, -7.8},
    "TC/AT", {-2.5, -6.3}, 
    "TT/AC", {-0.7, -1.2},
    "AA/TG", {-1.9, -4.4}, 
    "AG/TA", {-2.5, -5.9}, 
    "CA/GG", {-3.9, -9.6},
    "CG/GA", {-6.0, -15.5}, 
    "GA/CG", {-4.3, -11.1}, 
    "GG/CA", {-4.6, -11.4},
    "TA/AG", {-2.0, -4.7}, 
    "TG/AA", {-2.4, -5.8},
    "AG/TT", {-3.2, -8.7}, 
    "AT/TG", {-3.5, -9.4}, 
    "CG/GT", {-3.8, -9.0},
    "CT/GG", {-6.6, -18.7}, 
    "GG/CT", {-5.7, -15.9}, 
    "GT/CG", {-5.9, -16.1},
    "TG/AT", {-3.9, -10.5}, 
    "TT/AG", {-3.6, -9.8}
};

// Dangling ends table {DNA}
// Bommarito et al. {2000}, Nucl Acids Res 28, 1929-1934
Therm_Param Dangling_End[] = {
    "AA/.T", {0.2, 2.3}, 
    "AC/.G", {-6.3, -17.1}, 
    "AG/.C", {-3.7, -10.0},
    "AT/.A", {-2.9, -7.6}, 
    "CA/.T", {0.6, 3.3},
    "CC/.G", {-4.4, -12.6},
    "CG/.C", {-4.0, -11.9},
    "CT/.A", {-4.1, -13.0},
    "GA/.T", {-1.1, -1.6},
    "GC/.G", {-5.1, -14.0},
    "GG/.C", {-3.9, -10.9},
    "GT/.A", {-4.2, -15.0},
    "TA/.T", {-6.9, -20.0},
    "TC/.G", {-4.0, -10.9},
    "TG/.C", {-4.9, -13.8},
    "TT/.A", {-0.2, -0.5},
    ".A/AT", {-0.7, -0.8},
    ".C/AG", {-2.1, -3.9},
    ".G/AC", {-5.9, -16.5},
    ".T/AA", {-0.5, -1.1},
    ".A/CT", {4.4, 14.9},
    ".C/CG", {-0.2, -0.1},
    ".G/CC", {-2.6, -7.4},
    ".T/CA", {4.7, 14.2},
    ".A/GT", {-1.6, -3.6},
    ".C/GG", {-3.9, -11.2},
    ".G/GC", {-3.2, -10.4},
    ".T/GA", {-4.1, -13.1},
    ".A/TT", {2.9, 10.4},
    ".C/TG", {-4.4, -13.1},
    ".G/TC", {-5.2, -15.0},
    ".T/TA", {-3.8, -12.6}
};

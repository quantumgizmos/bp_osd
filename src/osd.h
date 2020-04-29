/* OSD.H - Interface to decoding procedures. */


typedef struct osd_search_strategy
{
    char *encoding_op_input;
} osd_search_strategy;


typedef struct osd_struct
{
    int osd_methods_count;
    char *decoding;
    int order;
    char name[50];
    int index;
    unsigned long int strategy_size;
    osd_search_strategy *strategy;
    unsigned long int success;
    unsigned long int osd_logical_failure;
    double logical_failure_rate;
    double logical_failure_rate_error_bar;

} osd_struct;



void soft_decision_col_sort(double *soft_decisions,int *cols, int N);

void reorder_char(char *bit_char,
char *reordered_bit_char,
int *indices);

int cmp(const void *a, const void *b);

osd_struct *create_osd_e_struct(mod2sparse *H, int w, int k);
osd_struct *create_osd_cs_struct(mod2sparse *H, int w);



int mod2sparse_decomp_osd
( mod2sparse *A,	/* Input matrix, M by N */
  int R,		/* Size of sub-matrix to find LU decomposition of */
  mod2sparse *L,	/* Matrix in which L is stored, M by R */
  mod2sparse *U,	/* Matrix in which U is stored, R by N */
  int *rows,		/* Array where row indexes are stored, M long */
  int *cols,		/* Array where column indexes are stored, N long */
  mod2sparse_strategy strategy, /* Strategy to follow in picking rows/columns */
  int abandon_number,	/* Number of columns to abandon at some point */
  int abandon_when	/* When to abandon these columns */
);

// void osd_0
// ( mod2sparse *A,	/* Input matrix, M by N */
//   char *synd,
//   double *log_prob_ratios,
//   int R,		/* Size of sub-matrix to find LU decomposition of */
//   char *decoding,
//   int *rows,		/* Array where row indexes are stored, M long */
//   int *cols,		/* Array where column indexes are stored, N long */
//   mod2sparse_strategy strategy, /* Strategy to follow in picking rows/columns */
//   int abandon_number,	/* Number of columns to abandon at some point */
//   int abandon_when	/* When to abandon these columns */
// );


int osd_e
( mod2sparse *A,	/* Input matrix, M by N */
  char *synd,
  double *log_prob_ratios,
  int R,		/* Size of sub-matrix to find LU decomposition of */
  osd_struct *osd_data
);


 int osd_cs
( mod2sparse *A,	/* Input matrix, M by N */
  char *synd,
  double *log_prob_ratios,
  int R,		/* Size of sub-matrix to find LU decomposition of */
  osd_struct *osd_data
);








void LU_solve_osd
(mod2sparse *L,
mod2sparse *U,
int *rows,
int *cols,
char *z,
char *x);
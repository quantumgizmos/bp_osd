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






int osd_e
( mod2sparse *A,	/* Input matrix, M by N */
  char *synd,
  double *log_prob_ratios,
  int A_rank,		/* Size of sub-matrix to find LU decomposition of */
  osd_struct *osd_data
);


 int osd_cs
( mod2sparse *A,	/* Input matrix, M by N */
  char *synd,
  double *log_prob_ratios,
  int R,		/* Size of sub-matrix to find LU decomposition of */
  osd_struct *osd_data
);









void InitExprEval();				/* initialize the expression evaluator */
void DestroyExprEval();				/* clean up everything */
int ExprEval(char * expr, char *result);	/* evaluate an expression */
void DefineVar(char *name, double value);	/* define a variable */
extern char **var_names;			/* variables names. */
extern double *var_values;			/* variables values. */
extern int var_count;				/* Number of variables. */

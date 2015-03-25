void InitExprEval();				/* initialize the expression evaluator */
void DestroyExprEval();				/* clean up everything */
int ExprEval(char * expr, char *result);	/* evaluate an expression */
void DefineVar(char *name, double value);	/* define a variable */

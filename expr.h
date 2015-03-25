void InitMathEval();				/* initialize the expression evaluator */
void DestroyMathEval();				/* clean up everything */
void ExprEval(char * expr, char *result);	/* evaluate an expression */
void DefineVar(char *name, double value);	/* define a variable */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <matheval.h>
#define BUFFER_SIZE 256

char **var_names;			/* variables names. */
double *var_values;			/* variables values. */
int var_count;				/* Number of variables. */

void InitExprEval()
{
	var_count=0;
	var_names=malloc(sizeof(char *));
	var_values=malloc(sizeof(double));
}

void DestroyExprEval()
{
	int i;
	for (i=0;i<var_count;i++)
		free(var_names[i]);
	free(var_names);
	free(var_values);
}



int NameIndex(char *name)
{
	int nl, i;
	nl=strlen(name);
	for (i=0;i<var_count;i++)
	{
		if (nl==strlen(var_names[i]))
			if (strncmp(name, var_names[i], nl)==0)
				return i;
	}
	return i;
}


int IsDefined(char *name)
{
	return (var_count!=NameIndex(name));
}

void DefineVar(char *name, double value)
{
	int len;
	len=strlen(name)+1;
	if (IsDefined(name))
	{
		int i;
		fprintf(stderr, "Warning: redefining variable %s\n",name);
		i=NameIndex(name);
		var_names[i]=realloc(var_names[i],(len+1)*sizeof(char));
		strncpy(var_names[i], name, strlen(name)+1);	
		var_values[i]=value;		
	}
	else
	{
		var_count++;
		var_names=realloc(var_names, (var_count+1)*sizeof(char *));
		var_names[var_count-1]=malloc((len+1)*sizeof(char));
		strncpy(var_names[var_count-1], name, strlen(name)+1);
	
		var_values=realloc(var_values, (var_count+1)*sizeof(double));
		var_values[var_count-1]=value;
	}
}

int ExprEval(char * expr, char *result)
{
	int i;
#ifdef WITH_LIBMATHEVAL
	void *f;
	char **names;
	int count;
	f = evaluator_create (expr);
	if (!f)
	{
		fprintf(stderr, "ExprEval: invalid input\n");
		return 1;
	}
		
	
	evaluator_get_variables (f, &names, &count);
	for (i=0;i<count;i++)
	{
		if (!IsDefined(names[i]))
		{
			fprintf(stderr, "ExprEval: Error variable \"%s\" not defined\n", names[i]);
			evaluator_destroy (f);
			return 1;
		}
	}
	snprintf(result, BUFFER_SIZE, "%.14g", evaluator_evaluate(f, var_count, var_names, var_values));

	evaluator_destroy (f);
#else
	i=NameIndex(expr);
	if (i<var_count)
		snprintf(result, BUFFER_SIZE, "%g.14", var_values[i]);
	else
	{
		fprintf(stderr, "ExprEval: Error variable \"%s\" not defined\n", expr);
		return 1;
	}
#endif
	return 0;
}

/*
int main (int argc, char **argv)
{
	char *result;
	
	result=malloc(BUFFER_SIZE*sizeof(char));
	InitMathEval();
	DefineVar("a", 10);
	DefineVar("b", 2);
	DefineVar("c", 5);
	ExprEval("a", result);
	printf("a=%s\n",result);
	ExprEval("b", result);
	printf("b=%s\n",result);
	ExprEval("c", result);
	printf("c=%s\n",result);
	ExprEval("a*e^2+b*e+c+d", result);
	printf("a*e^2+b*e+c+d=%s\n",result);	
	DestroyMathEval();
	free(result);
}
 */

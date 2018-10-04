#include <stdio.h> 
#include <math.h>


// Function to create header file with beta values hard-coded
void calculateBetas(const int maxIt, const int use_double)
{

	int i;

	FILE *f;

	double alpha, alpha_new, beta;

	// create file
	f = fopen("include/betas.h", "w+");

	if (f == NULL)
	{
		printf("Error opening file!\n");
	}

	// declare vector of betas
	if (use_double == 1 )
	{
		fprintf(f, "real_t betas[%d] = {\n", maxIt);
	}
	else
	{
		fprintf(f, "real_t betas[%d] = {\n", maxIt);
	}

	alpha = 1.0;

	// calculate betas
	for (i=0;i<maxIt;i++)
	{

		alpha_new = (alpha/2.0)*(sqrt(alpha*alpha+4.0)-alpha);
        beta      = (alpha*(1.0-alpha))/(alpha*alpha + alpha_new);
        alpha     = alpha_new;

        if (i < maxIt-1 )
        {
			fprintf(f, "%f,\n", beta);	
		}
		else
		{
			fprintf(f, "%f };\n", beta);			
		}	
	}

	fclose(f);

}

// Main function for code-generation
int main(int argc, char const *argv[])
{

	int maxIt = 1000;
	int use_double = 1;

	calculateBetas(maxIt, use_double);

	return 0;
}


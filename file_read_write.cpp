# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>

FILE *io, *tamtesttime;
int NO_OF_CORES[5] = {32, 19, 28, 4, 10};
int NO_OF_TAMS[5];
int create_io_tam ( int TAM_WIDTH_MAX, int NO_OF_MODULE)
{
	
	//int module_no;
	int tam;
	long int testtime;
	for(int i = 1;i<= NO_OF_MODULE; i++)
	{
		FILE *read;
		char input_file [100]="sic_level_";
		char r[10];
		int core_count = 0;
		int total_tam = 0;
		sprintf (r, "%d",i);
		strcat (input_file, r);
		strcat (input_file, ".txt");
		read = fopen(input_file,"r");
		fprintf(io,"{");
		int count = 0;
		while(!feof(read))
		{
			
			fscanf(read,"%d\t%ld\n", &tam, &testtime);
			
			
			if(tam == 1)
			{
				if(count != 0)
				{
					fprintf(io,"%d",count);	
					fprintf(io,",");
					total_tam = total_tam + count;
				}

				count = 0;
				core_count = core_count + 1;
				
			}
			

			if(tam <= TAM_WIDTH_MAX)
			{
				fprintf(tamtesttime,"%d\t%ld\n", tam, testtime);
				count++;
			}

		}
		printf("%d\n", core_count);

		fprintf(io,"%d",count);	
		total_tam = total_tam + count;
		fprintf(io,"},");

		NO_OF_TAMS[i-1] = total_tam;

	}
	return 0;
}

void create_power ( int NO_OF_MODULE)
{
	
	//int module_no;
	int tam;
	long int testtime;
	for(int i = 1;i<= NO_OF_MODULE; i++)
	{
		FILE *read;
		char input_file [100]="sic_level_core_wise_power_";
		char r[10];
		double power;
		sprintf (r, "%d",i);
		strcat (input_file, r);
		strcat (input_file, ".txt");
		read = fopen(input_file,"r");
		// fprintf(io,"{");
	
		for (int j = 0; j< NO_OF_CORES[i-1]; j++)
		{	
			fscanf(read,"%lf\t", &power);
			fprintf(io,"%ld", (long int)power);
			if(j < NO_OF_CORES[i-1]-1)
			{
				fprintf(io, ",");
			}
		}
		fprintf(io, ",");
		

		
		// fprintf(io,"},");

	}
}



int main( int argc, char *argv [ ] )
{
		
	int TAM_WIDTH_MAX = atoi (argv[1]);
	int NO_OF_MODULE = atoi (argv [3]);
	int TSV_MAX = atoi (argv [2]);
	int HARD_DIE_TEST = atoi(argv[4]);
	int tamsum;
	


	io = fopen("iolist.h","w");
	tamtesttime = fopen("tam_testtime.txt", "w");
	fprintf(io,"#define TAM_WIDTH_MAX %d\n", TAM_WIDTH_MAX);
	if(TAM_WIDTH_MAX <= 16)
		fprintf(io,"#define BIT_LENGTH 4\n" );
	else if(TAM_WIDTH_MAX > 16 && TAM_WIDTH_MAX <= 32)
		fprintf(io,"#define BIT_LENGTH 5\n" );
	else if(TAM_WIDTH_MAX > 32 && TAM_WIDTH_MAX <= 64)
		fprintf(io,"#define BIT_LENGTH 6\n" );
	else if(TAM_WIDTH_MAX > 64)
		fprintf(io,"#define BIT_LENGTH 7\n" );
	fprintf(io,"#define TSV_MAX %d\n",TSV_MAX);
	fprintf(io,"#define NDies %d\n",NO_OF_MODULE);
	fprintf(io,"#define NCores {");
	int total_cores = 0;
	for(int i = 0; i < NO_OF_MODULE; i++)
	{
		fprintf(io,"%d", NO_OF_CORES[i]);
		total_cores = total_cores + NO_OF_CORES[i];
		if(i < NO_OF_MODULE-1)
			fprintf(io,",");	
	}
	fprintf(io,"}\n");
	
	fprintf(io,"#define SIZE %d\n",total_cores);

	fprintf(io,"#define HARD_DIE_TEST %d\n",HARD_DIE_TEST);

	fprintf(io,"#define TAM {");
	tamsum = create_io_tam ( TAM_WIDTH_MAX,NO_OF_MODULE);
	fprintf(io,"}\n");

	
	fprintf(io,"#define TAMSUM {");
	for(int i = 0; i < NO_OF_MODULE; i++)
	{
		fprintf(io,"%d", NO_OF_TAMS[i]);
		
		if(i < NO_OF_MODULE-1)
			fprintf(io,",");	
	}
	fprintf(io,"}\n");
	
	fprintf(io,"#define POWER {");
	create_power( NO_OF_MODULE);
	fseek(io,-1, SEEK_END);
	fprintf(io,"}\n");

	fclose(io);
	fclose(tamtesttime);

}

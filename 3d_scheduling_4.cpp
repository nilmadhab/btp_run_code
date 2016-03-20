# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "iolist.h"
#include <iostream>
#define LARGENUMBER 999999999
#define  NO_OF_PARTICLES 1000
#define  MAX_ITERATION 1000
#define ALPHA 0.1
#define BETA 0.10
#define  GAMMA 0.2
#define MAX_FREQS 5 //index 0 is normal frequency, 1 is 1/2 and so on
#define NEW_PARTICLE_SIZE BIT_LENGTH*SIZE+SIZE
#define FREQ_START_INDEX 2 //freq_factor[i] = pow(2,(FREQ_START_INDEX-freq_index[i]));
using namespace std;


typedef struct swapsequence{
int from;
int to;
} swapseq;

typedef struct particlestructure{
int info[NEW_PARTICLE_SIZE];
int lbest[NEW_PARTICLE_SIZE];
long int lbesttime;
long int time_fitness;
} particle;

typedef struct {
        int corenum;
        int tam_width ;
        long int starttime ;
        long int endtime ;
        double power;
        float freq;
        int tsv;
} schedulerInfo ;

typedef struct 
{
	int tam_list[TAM_WIDTH_MAX];
	long int testtime_list[TAM_WIDTH_MAX];
	int die_level;
	int no_of_tam;
	double peak_power;
} tam_testtime_information;

tam_testtime_information tam_testtime[SIZE];
schedulerInfo scheduler [ NO_OF_PARTICLES ] [ SIZE ];
particle partarray[NO_OF_PARTICLES];
particle globalbest;
schedulerInfo BestParticle [SIZE];

double power_budget;
double FIXED_POWER[] = POWER;

int TAM_INFO[][50] = TAM;
int Cores[] = NCores;
int globalbestIndx;
int particle_test = 1;
//long int temptime=LARGENUMBER;

//Hard die config for d695 as lowest die
int hard_tam_width[SIZE] = {7, 10, 12, 12, 15};
// long int hard_test_time[SIZE] = {94731, 605630, 556308, 1267313, 1835731};  //PSO 2D
long int hard_test_time[SIZE] = {96297, 669329, 651281, 1384949, 1947063};

void print_array(int length, int array[], char c)
{
    for ( int i = 0 ; i < length ; i ++ ) 
	{
		printf("%d ", array[i]);
	}
	printf("%c", c);
}
void init_array(int length, int array[])
{
    for ( int i = 0 ; i < length ; i ++ ) 
	{
		array[i] = 0;
	}
}
int bin_to_decimal(int binary)
{
	//int binary = 0;
	int decimal = 0;
	int remainder;
	int power = 1;
	
	while(binary!=0)
	{
        	remainder = binary % 10;
        	decimal = decimal + remainder * power;
        	power = power * 2;
        	binary=binary/10;
    	}	
	return decimal;
}

int tsv_required(int parallel_core_index[], int assigned_core_tam[], int layer_tsv[], int corelevel[])
{	

	int total_tsv = 0;
	int lowest_die_test = 0;
	int parallel_index[NDies];
	int assigned_tam[NDies];

	int offset = 0;
	//find parallel dies and total tam width assigned to each die
	
	for(int i = 0; i < NDies; i++)
	{
		assigned_tam[i] = 0;
		parallel_index[i] = 0;
	}
	/*
	for(int i = 0; i < SIZE; i++)
	{
		printf("%d ", assigned_core_tam[i]);
	}
	printf("\n");
	for(int i = 0; i < SIZE; i++)
	{
		printf("%d ", parallel_core_index[i]);
	}
	printf("\n");
*/
	for(int i = 0 ; i < SIZE; i++)
	{
		assigned_tam[corelevel[i]] = assigned_tam[corelevel[i]] + assigned_core_tam[i]*parallel_core_index[i];
		if(parallel_core_index[i])
			parallel_index[corelevel[i]] = parallel_core_index[i];
	}
	/*
	for(int i = 0; i < NDies; i++)
	{
		printf("%d ", assigned_tam[i]);
	}
	printf("\n");
	for(int i = 0; i < NDies; i++)
	{
		printf("%d ", parallel_index[i]);
	}
	printf("\n");
	*/
	//check if lowest die is being tested
	if(parallel_index[0])
	{
		lowest_die_test = 1;
		for(int i = 1; i < NDies; i++)
		{
			if(parallel_index[i])
			{
				lowest_die_test = 0;
			}
		}
	}

	//no TSV needed to test lowest die
	if(lowest_die_test)
	{
		return 0;
	}
	for(int i = 1; i < NDies; i++)
	{
		int temp[NDies]; 
		int count = 0;
		for(int jj = 0; jj < NDies; jj ++)
		{
			temp[jj] = 0;
		}
		for(int j = i; j < NDies; j++)
		{
			int parallel_tsv_count = 0;
			for(int k = j; k < NDies; k++)
			{
				parallel_tsv_count = parallel_tsv_count +  assigned_tam[k] * parallel_index[k];
			}
			
			if(parallel_tsv_count >= assigned_tam[j])
				temp[count] = parallel_tsv_count;
			else
				temp[count] = assigned_tam [j];
			count++;
		}

		
		int max = 0;
		for(int ii = 0; ii < NDies; ii++)
		{
			if(temp[ii] >= max)
			{
				max = temp[ii];
			}
		}
		layer_tsv[i-1] = max;
		//printf("Layer %d : %d ",i+1, max);

		total_tsv = total_tsv + max;
	}
	
	//printf("TSV used : %d\n", total_tsv);
	return total_tsv;
}

long int bin_packing(int *particle_info, int index, int genaration)
{
	// printf("test1\n");
	int binary;
	long int testtime[SIZE];
	int assigned_tam[SIZE];
	int assigned_tsv[SIZE];
	double area[SIZE];
	int schedule_index[SIZE];
	int no_of_core_scheduled = 0;
	int tam_index[SIZE];
	long int break_point[SIZE+2];
	int available_tam_width[SIZE+2];
	int insert_index;
	long int temp;
	double max_weightage;
	double min_weightage;
	int max_index;
	int voltage_index;
	int layer_tsv[NDies-1];
	int session_tsv[NDies-1];
	int session = 1;
	int parallel_session_index[NDies][NDies];
	int parallel_temp_index[NDies];
	int corelevel[SIZE];
	double assigned_peak_power[SIZE];
	double available_power_left[SIZE+2];
	int freq_index[SIZE];
	float freq_factor[SIZE];
	//printf("************Begin schedule************\n");

	for(int  i = 0; i< SIZE; i++)
		scheduler [index][i].starttime = -1;
	
	/*for(int i = 0; i< SIZE; i++)
	{
		tam_index[i] = (int)(particle_info[i] * tam_testtime[i].no_of_tam);
	}*/
	float voltage_levels[5] =  {1, 0.9, 0.8, 0.7, 0.6};

	for(int i = 0; i< SIZE; i++)
	{
		binary = 0;
		for(int k = 0; k < BIT_LENGTH ; k++)
			binary = binary + particle_info[(i* BIT_LENGTH) + k] * pow (10,BIT_LENGTH - 1 -k);
		tam_index[i] = bin_to_decimal(binary);
		freq_index[i] = particle_info[i+BIT_LENGTH*SIZE];
		freq_factor[i] = pow(2,(FREQ_START_INDEX-freq_index[i]));
	}		
	// printf("test2\n");

	/*for(int i = 0; i< SIZE; i++)
		printf("%d\t",tam_index[i]);
	printf("\n\n\n");	*/	
	if(HARD_DIE_TEST)
	{
		for(int i = 0; i< SIZE; i++)
		{
			assigned_tam[i] = hard_tam_width[i];
			testtime[i] = hard_test_time[i];
			area[i] = assigned_tam[i] * testtime[i];
		}
	}
	else
	{
		for(int i = 0; i< SIZE; i++)
		{

			assigned_tam[i] = tam_testtime[i].tam_list[tam_index[i]];
			testtime[i] = tam_testtime[i].testtime_list[tam_index[i]]/freq_factor[i];  //testtime twiced for 1/2 factor
			assigned_peak_power[i] = tam_testtime[i].peak_power * freq_factor[i]; //* pow(voltage_levels[freq_index[i]], 2);  //power halved for 1/2 factor
			corelevel[i] = tam_testtime[i].die_level;
			area[i] = assigned_tam[i] * testtime[i];

		}
	}

	//printf("\n i am here\n");
	for(int i = 0; i< SIZE; i++)
			schedule_index[i] = 0;
			
	for(int i = 0; i<= SIZE+1; i++)
		available_tam_width[i] = TAM_WIDTH_MAX;
	for(int i = 0; i<= SIZE+1; i++)
		available_power_left[i] = power_budget;		

	break_point[0] = 0;
	for(int i = 1; i<= SIZE+1; i++)
		break_point [i] = LARGENUMBER;
	int min_tsv[NDies-1];
	int max_tam = 0;
	init_array(NDies-1,  min_tsv);
	for(int i = 0; i < SIZE; i++)
	{
		if(assigned_tam[i] > max_tam)
			max_tam = assigned_tam[i];
		if(tam_testtime[i].die_level > 0 && min_tsv[tam_testtime[i].die_level-1] < assigned_tam[i])
		{
			min_tsv[tam_testtime[i].die_level-1] = assigned_tam[i];
		}
		

	}
	for(int i = 0; i< NDies-1; i++)
		session_tsv[i] = layer_tsv [i] = min_tsv[i];

	int t =0;
	
	long int schedule_start_time = 0;
	while (no_of_core_scheduled < SIZE)
	{
		int flag = 1;

		//If no breakpoints are left then the schedule is invalid due to TSV limit
		int number_of_brk = 0;

		while(break_point [number_of_brk] != LARGENUMBER)
			number_of_brk++;
		if(number_of_brk == 0 && schedule_start_time)
		{
			//printf("%d\n", -1);
			return -1;
		}
		while(flag ==1)
		{
			max_weightage = 0;
			voltage_index = 0;
			min_weightage = LARGENUMBER;
			flag = 0;
			int jj = 0;
			schedule_start_time = break_point[jj];

			
			int current_layer_tsv[NDies-1];
			if(session)
			{
				for(int i = 0; i< NDies-1; i++)
					current_layer_tsv [i] = 0;
			}
			for (int kk = 0; kk< SIZE; kk++)
			{

				
				if(schedule_index[kk] == 0 )
				{
					
					if (available_tam_width[jj] >= ((int)(assigned_tam[kk]*freq_factor[kk])) && available_power_left[jj] >= assigned_peak_power[kk] * pow(voltage_levels[voltage_index], 2))
					{
						int parallel_index[SIZE];
						int finish_index[SIZE];
						for(int aa = 0; aa < SIZE; aa++)
						{
							parallel_index[aa] = 0;
							finish_index[aa] = 0;
						}
						
						for(int zz = 0; zz < SIZE; zz++)
						{
							// Find parallel cores
							if(schedule_index[zz] == 1)
							{
								for(int yy = 0; yy < t; yy++)
								{
									if(scheduler[index][yy].corenum == zz)
									{
										if(schedule_start_time >= scheduler [index][yy].starttime && schedule_start_time < scheduler [index][yy].endtime)
											parallel_index[zz] = 1;	

										if(schedule_start_time >= scheduler [index][yy].endtime)
											finish_index[zz] = 1;
										
											
									}
								}
								
							}
							
							
							

						}

						parallel_index[kk] = 1; 

						
						int tsv_count = tsv_required(parallel_index, assigned_tam, current_layer_tsv, corelevel);
						
						/*for(int i = 0; i < SIZE; i++)
						{
							printf("%d ", finish_index[i]);
						}
						printf("\n");
						for(int i = 0; i < SIZE; i++)
						{
							printf("%d ", parallel_index[i]);
						}
						printf("\n");*/
						/*
						for(int i = 0; i < SIZE; i++)
						{
							//printf("%ld ", testtime[i]);
						}
						//printf("\n");
						for(int i = 0; i < SIZE; i++)
						{
							printf("%d ", assigned_tam[i]);
						}
						printf("\n");
						*/
						
						// print_array(SIZE, finish_index, '\n');
						// print_array(NDies-1, current_layer_tsv, '\n');
						// print_array(NDies-1, session_tsv, '\n');
						
						// printf("TSV used : %d\n", tsv_count);
						

						int cost = area[kk];// * (corelevel[kk]+1);// * (TSV_MAX - tsv_count);
						int temp_voltage_index = freq_index[kk];
						// printf("%d\n", cost);
						cost = break_point[jj]- testtime[kk];				
						int temp_layer_tsv[NDies-1];
						for(int pp = 0; pp < NDies-1; pp++)
						{
							temp_layer_tsv[pp] = 0;
						}
						
						int layer_check = 1;
						int sum = 0;
						for (int pp = 0; pp < NDies-1; pp++)
						{
							if(current_layer_tsv[pp] < session_tsv[pp])
							{
								temp_layer_tsv[pp] = session_tsv[pp];
							}
							else if(current_layer_tsv[pp] >= session_tsv[pp])
							{
								temp_layer_tsv[pp] = current_layer_tsv[pp];
							}
							sum = sum + temp_layer_tsv[pp];
						}
						// printf("TSV used : %d\n", sum);
						if(sum <= TSV_MAX || tsv_count == 0)
						{
							
							
							flag = 1;
							/*if(cost > max_weightage)
							{							
								
								max_weightage = cost;
								max_index = kk;
								for (int pp = 0; pp < NDies-1; pp++)
								{
									layer_tsv[pp] = temp_layer_tsv[pp]; 
								}
								if(temp_voltage_index < voltage_index)
								{
									//less index means higher voltage
									voltage_index = temp_voltage_index;
								}
								
							}*/

							if(cost < min_weightage)
							{							
								
								min_weightage = cost;
								max_index = kk;
								for (int pp = 0; pp < NDies-1; pp++)
								{
									layer_tsv[pp] = temp_layer_tsv[pp]; 
								}
								if(temp_voltage_index < voltage_index)
								{
									//less index means higher voltage
									voltage_index = temp_voltage_index;
								}
								
							}
										
					
						}

						
						
					}
				}
			
				
			}
			
			if (flag == 0)
			{
				
				break;
			}
			
			
			scheduler [index][t].corenum =  max_index;
			scheduler [index][t].tam_width = assigned_tam[max_index];
			scheduler [index][t].starttime = break_point[jj];
			scheduler [index][t].endtime = break_point[jj] + testtime [max_index];	
			scheduler [index][t].power = assigned_peak_power[max_index];
			scheduler [index][t].freq = freq_factor[max_index];
			
			int number_of_brk = 0;
			while(break_point [number_of_brk] != LARGENUMBER)
				number_of_brk++;
			if(scheduler [index][t].endtime > break_point[number_of_brk - 1])
			{
				//session based
				if(number_of_brk == 1)
					break_point[number_of_brk] = scheduler [index][t].endtime;
				else
					break_point[number_of_brk-1] = scheduler [index][t].endtime;
			}
			// break_point[number_of_brk] = scheduler [index][t].endtime;
			for (int i = 0; i<=(SIZE); i++)
			{
				for (int j = 0; j<=(SIZE-i); j++)
				{
					if(break_point[j]> break_point[j+1])
					{
						temp = break_point[j];
						break_point[j] = break_point[j+1];
						break_point[j+1] = temp;
					}
				}
			}
			for (int i = 0; i<= SIZE+1; i++)
			{
				if (break_point[i] >= scheduler[index][t].endtime  )
				// if (break_point[i] == scheduler[index][t].endtime  )
				{
					insert_index = i;
					//break;
				}	
			}
			//inserting breakpoint
			//insert_index = number_of_brk;
			for (int i = number_of_brk; i>= insert_index; i--)
			{
				available_tam_width[i] = available_tam_width[i-1];	
				available_power_left[i] = available_power_left[i-1];	
			}
			int a = jj;
			// if(!particle_test)
				// printf("Session voltage %d\n", voltage_index);
			voltage_index = 0;
			while(break_point[a]< scheduler[index][t].endtime )
			{		
				available_tam_width[a] = available_tam_width[a] - (int)(assigned_tam[max_index]*freq_factor[max_index]);
				available_power_left[a] = available_power_left[a] - assigned_peak_power[max_index]*pow(voltage_levels[voltage_index], 2);
				a++;
			}	
			schedule_index[scheduler [index][t].corenum] = 1;
			no_of_core_scheduled++;	
			t++;
		}
		
		
		int latest_schedule_cores[SIZE];
		long int latest_schedule_core_testtime[SIZE];
		int latest_schedule_core_tam[SIZE];	
		int k = 0;
		int no_of_latest_scheduled_core = 0;
		int temp1;
		long int temp2;
		for(int i = 0; i< SIZE; i++)
		{
			//if(schedule_index[scheduler [index][i].corenum] == 1)
			if (scheduler [index][i].starttime == break_point[0] )
			{
				latest_schedule_cores[k] = scheduler [index][i].corenum;
				latest_schedule_core_testtime[k] = scheduler [index][i].endtime;
				latest_schedule_core_tam[k] = scheduler [index][i].tam_width;
				k++;
				no_of_latest_scheduled_core++;	
			}		
		}
		for(int i = 0; i< no_of_latest_scheduled_core-1; i++)
		{
			for(int j = 0 ; j< no_of_latest_scheduled_core-i-1; j++)
			{
				if(latest_schedule_core_testtime[j] <latest_schedule_core_testtime[j+1])
				{
					temp1 = latest_schedule_cores[j];
					temp2 = latest_schedule_core_testtime[j];
					latest_schedule_cores[j] = latest_schedule_cores[j+1];
					latest_schedule_core_testtime[j] = latest_schedule_core_testtime[j+1];
					latest_schedule_cores[j+1] = temp1;
					latest_schedule_core_testtime[j+1] = temp2;
				}
			}
		}
	//if(genaration == MAX_ITERATION -1)
	//{
		int core_check = 0;
		if(available_tam_width[0] > 0)
		{	
			/*int kk = 0;
			while((no_of_latest_scheduled_core - core_check) > 0 && available_tam_width[0] > 0)
			{		
				int extra_tam = available_tam_width[0];
				int break_count = 0;
				int core_number;
				int core_index_in_scheduler;
				int position_in_break_point;
				int current_tam_index, next_tam_index;
				int current_tam, next_tam;
				long int current_testtime, next_testtime;
				while(break_point [break_count] != LARGENUMBER)
					break_count++;
				for (int i = 0; i< break_count; i++)
				{
					if(break_point[i] == latest_schedule_core_testtime[kk])
					{
						position_in_break_point = i;
						break;
					}	
				}
					
				for (int i = 0; i< SIZE; i++)
				{
					if(scheduler[index][i].corenum == latest_schedule_cores[kk])
					{
						core_number = latest_schedule_cores[kk];
						core_index_in_scheduler = i;
						break;
					}
				}
				for(int i = 0; i< SIZE; i++)
				{
					if(i == core_number)
					{
						current_tam_index = tam_index[i];
						break;
					}
				}
				current_tam = tam_testtime[core_number]. tam_list[current_tam_index];
				current_testtime = tam_testtime[core_number]. testtime_list[current_tam_index];
				next_tam_index = current_tam_index + 1;
				int check = 0;
				while(next_tam_index < tam_testtime[core_number].no_of_tam && tam_testtime[core_number].tam_list[next_tam_index]- current_tam <= extra_tam)
				{
					check = 1;
					next_tam_index ++;
				}
				if(check == 1)
				{
					tam_index[core_number] = next_tam_index-1;
					next_tam = tam_testtime[core_number].tam_list[next_tam_index-1];
					next_testtime = tam_testtime[core_number].testtime_list[next_tam_index-1];
					for(int i = position_in_break_point; i< SIZE+1; i++)
						break_point[i] = break_point[i+1];
					for(int i = 0; i< (position_in_break_point); i++)
					{
						available_tam_width[i] = available_tam_width[i]+ current_tam;
					}
					for(int i = position_in_break_point; i< SIZE+1; i++)
						available_tam_width[i] = available_tam_width[i+1];
					
					scheduler [index][core_index_in_scheduler].tam_width = next_tam;
					scheduler [index][core_index_in_scheduler].endtime = break_point[0] + next_testtime;	
					int number_of_brk1 = 0;
					while(break_point [number_of_brk1] != LARGENUMBER)
						number_of_brk1++;
					break_point[number_of_brk1] = scheduler [index][core_index_in_scheduler].endtime;
					long int temp3;
					for (int i = 0; i<=(SIZE); i++)
						for (int j = 0; j<=(SIZE-i); j++)
						{
							if(break_point[j]> break_point[j+1])
							{
								temp3 = break_point[j];
								break_point[j] = break_point[j+1];
								break_point[j+1] = temp3;
							}
						}	
					for (int i = 0; i<= SIZE+1; i++)
					{
						if (break_point[i] == scheduler [index][core_index_in_scheduler].endtime )
						{
							insert_index = i;
							//break;
						}	
					}	
					for (int i = number_of_brk1; i>=insert_index; i--)
						available_tam_width[i] = available_tam_width[i-1];
					int a = 0;
					while(break_point[a]< scheduler [index][core_index_in_scheduler].endtime )
					{				
						available_tam_width[a] = available_tam_width[a] - next_tam;
						a++;
					}	
				}
				
				core_check++;
				kk++;
			}*/
		}
		
	//}	

		for (int i =0; i< SIZE+1; i++)
		{
			break_point[i] = break_point[i+1];
			available_tam_width[i] = available_tam_width[i+1];
			available_power_left[i] = available_power_left[i+1];
		}
		session++;
		for(int i = 0; i < NDies-1; i++)
		{		
			session_tsv[i] = layer_tsv[i];
			for(int kk = 0; kk<SIZE; kk++)
			{
				if(corelevel[scheduler [index][kk].corenum] == NDies-1)
					scheduler [index][kk].tsv = 0;
				else if(corelevel[scheduler [index][kk].corenum] == i)
				{
					scheduler [index][kk].tsv = layer_tsv[i];
				}
			}
		}
		//printf("\n");
		break_point[SIZE+1] = LARGENUMBER;
		available_tam_width[SIZE+1] = TAM_WIDTH_MAX;	
		available_power_left[SIZE+1] = power_budget;	
	}
	
	schedulerInfo tmp;
        for ( int i = 0 ; i < SIZE ; i ++ )
        {
                for ( int j = 0 ; j < SIZE-1-i ; j ++ ) 
                {
                        if ( scheduler [index][j+1].starttime <  scheduler [index][j].starttime ) 
                        {
                                tmp.corenum = scheduler [index][j+1]. corenum ;
                                tmp.tam_width = scheduler [index][j+1]. tam_width ;
                                tmp.starttime = scheduler [index][j+1]. starttime;
                                tmp.endtime = scheduler [index][j+1]. endtime;
                                tmp.tsv = scheduler [index][j+1]. tsv;
                                tmp.power = scheduler [index][j+1]. power;
                                tmp.freq = scheduler [index][j+1]. freq;
                                
                                scheduler [index][j+1]. corenum = scheduler [index][j]. corenum ;
                                scheduler [index][j+1]. tam_width = scheduler [index][j]. tam_width ;
                                scheduler [index][j+1]. starttime = scheduler [index][j]. starttime ;
                                scheduler [index][j+1]. endtime = scheduler [index][j]. endtime ;
                                scheduler [index][j+1]. tsv = scheduler [index][j]. tsv ;
								scheduler [index][j+1]. freq = scheduler [index][j]. freq;	


                                scheduler [index][j]. corenum = tmp.corenum ;
                                scheduler [index][j]. tam_width = tmp.tam_width ;
                                scheduler [index][j]. starttime = tmp.starttime ;
                                scheduler [index][j]. endtime = tmp.endtime ;
                                scheduler [index][j]. tsv = tmp.tsv;
                                scheduler [index][j]. freq =  tmp.freq;
                                
                        }
                }
        }
	
	
	int flag = 0;
	int kk;
	for ( kk = 0; kk<= SIZE+1; kk++)
		if(break_point[kk] == LARGENUMBER)
		{	
			flag = 1;
			break;
		}
	if (flag ==1)
	{
	
	if (break_point[kk-1] == 64)
	{	
		printf("\n********%ld\n",break_point[kk-1]);
		for(int i = 0; i<= SIZE+1; i++)
					printf("%ld\t", break_point[i]);
				printf("\n\n");
				for(int i = 0; i<= SIZE+1; i++)
					printf("%d\t", available_tam_width[i]);
				printf("\n\n");
	exit(0);
	}
	return(break_point[kk-1] );
	
	}
	else
	{
	
	if (break_point[kk] == 64)
	{	
		printf("\n++++++%ld\n",break_point[kk]);
		for(int i = 0; i<= SIZE+1; i++)
					printf("%ld\t", break_point[i]);
				printf("\n\n");
				for(int i = 0; i<= SIZE+1; i++)
					printf("%d\t", available_tam_width[i]);
				printf("\n\n");
	exit(0);
	}
	
	return (break_point[kk]);
	
	}
	
				
}				




/*void swaptamindex(double *a, double *b, double *c, double prob1, double prob2) //here prob is ALPHA-probal with which lbest is updated
{
	double rand1, rand2, cost;
	for (int i = 0; i< SIZE;)
	{
		rand1 = drand48();
		rand2 = drand48();
		cost = (prob1 * rand1 * (b[i] - a [i])) + (prob2 * rand2 * (c[i] - a [i]));
		a[i] = a[i]+ cost;
		if ((a[i]) >= 0 && a[i] < 1)
			i++;	
	}
}
*/
void swapparticleinfo(int *a, int *b, double prob)
{
	int binary;
	int decimal;
	int flag = 1;
	while(flag)
	{
		int tsv_count = 0;
		int core_temp[SIZE]; 
		for(int jj = 0; jj < SIZE; jj ++)
		{
			core_temp[jj] = 0;
		}
		int temp[NDies];
		for(int i = 0; i < NDies; i++)
			temp[i] = 0;
		for(int i = 0;i < (SIZE);)
		{
			//if(drand48() <= prob)
				for(int j = (i* BIT_LENGTH); j < ((i* BIT_LENGTH) + BIT_LENGTH); j++)	
					if(drand48() <= prob)
						a[j] = b[j];
			binary = 0;	
			for(int k = 0; k < BIT_LENGTH ; k++)
				binary = binary + a[(i* BIT_LENGTH) + k] * pow (10,BIT_LENGTH - 1 -k);
			decimal	= bin_to_decimal(binary);

			if(drand48() <= prob)
			{
				a[i+BIT_LENGTH*SIZE] = b[i+BIT_LENGTH*SIZE];
			}

			int freq_index =  a[i+BIT_LENGTH*SIZE];
			float freq_factor = pow(2,(FREQ_START_INDEX-freq_index));
			//printf("%d %d %d %d %lf\n", i, freq_index, decimal, tam_testtime[i].tam_list[decimal] , tam_testtime[i].peak_power * freq_factor);

			if(decimal< tam_testtime[i]. no_of_tam && tam_testtime[i].tam_list[decimal] * freq_factor < TAM_WIDTH_MAX && tam_testtime[i].peak_power * freq_factor < power_budget)
			{
				core_temp[i] = tam_testtime[i].tam_list[decimal];
				i++;
			}	
		}

		

		//calculate max tam width in each die
		
		for(int i = 0; i < SIZE; i++)
		{
			
			if(temp[tam_testtime[i].die_level] < core_temp[i])
			{
				temp[tam_testtime[i].die_level] = core_temp[i];
			}
			

		}

		//calculating minimum number of total TSVs needed by calculating minimum tsv required at each layer. 
		//min tsv at each layer is the max tam width assigned to layers above
		for(int i = 1; i < NDies; i++)
		{			
			int max = 0;
			for(int ii = i; ii < NDies; ii++)
			{
				if(temp[ii] >= max)
				{
					max = temp[ii];
				}
			}
			
			tsv_count = tsv_count + max;
		}
		if(tsv_count <= TSV_MAX)
		{
			flag = 0;
			
		}
		
	}
	/*for(int i = 0;i < (SIZE);)
		{
			
			for(int k = 0; k < BIT_LENGTH ; k++)
				binary = binary + a[(i* BIT_LENGTH) + k] * pow (10,BIT_LENGTH - 1 -k);
			decimal	= bin_to_decimal(binary);
			int freq_index =  b[i+BIT_LENGTH*SIZE];
			float freq_factor = pow(2,(FREQ_START_INDEX-freq_index));
			if(tam_testtime[i].tam_list[decimal] * freq_factor < TAM_WIDTH_MAX && tam_testtime[i].peak_power * freq_factor < power_budget)
			{
				a[i+BIT_LENGTH*SIZE] = freq_index;
				i++;
			}
		}*/
}

void initialiseparticle(particle *object)
{
	int binary;
	int decimal;
	object->lbesttime = LARGENUMBER;
	
	int flag = 1;
	while(flag)
	{
		int tsv_count = 0;
		int core_temp[SIZE]; 
	

		for(int jj = 0; jj < SIZE; jj ++)
		{
			core_temp[jj] = 0;
		}
		int temp[NDies];
		for(int i = 0; i < NDies; i++)
			temp[i] = 0;
		for(int i = 0;i < (SIZE);)
		{
			
			for(int j = (i* BIT_LENGTH); j < ((i* BIT_LENGTH) + BIT_LENGTH); j++)	
				object->info[j]=rand()% 2;
			binary = 0;	
			for(int k = 0; k < BIT_LENGTH ; k++)
				binary = binary + object->info[(i* BIT_LENGTH) + k] * pow (10,BIT_LENGTH - 1 -k);
			decimal	= bin_to_decimal(binary);
			int freq_index = rand() % MAX_FREQS;
			float freq_factor = pow(2,(FREQ_START_INDEX-freq_index));

			if(decimal< tam_testtime[i]. no_of_tam && tam_testtime[i].tam_list[decimal] * freq_factor < TAM_WIDTH_MAX && tam_testtime[i].peak_power * freq_factor < power_budget)
			{
				core_temp[i] = tam_testtime[i].tam_list[decimal];
				object->info[i+BIT_LENGTH*SIZE] = freq_index;
				i++;
			}
		}
		

		

		//calculate max tam width in each die
		
		for(int i = 0; i < SIZE; i++)
		{
			
			if(temp[tam_testtime[i].die_level] < core_temp[i])
			{
				temp[tam_testtime[i].die_level] = core_temp[i];
			}
			

		}

		//calculating minimum number of total TSVs needed by calculating minimum tsv required at each layer. 
		//min tsv at each layer is the max tam width assigned to layers above
		for(int i = 1; i < NDies; i++)
		{			
			int max = 0;
			for(int ii = i; ii < NDies; ii++)
			{
				if(temp[ii] >= max)
				{
					max = temp[ii];
				}
			}
			
			tsv_count = tsv_count + max;
		}
		if(tsv_count <= TSV_MAX)
		{
			flag = 0;
			
			
		}
		
	}	

	/*for(int i = 0;i < (SIZE);)
		{
			
			for(int k = 0; k < BIT_LENGTH ; k++)
				binary = binary + object->info[(i* BIT_LENGTH) + k] * pow (10,BIT_LENGTH - 1 -k);
			decimal	= bin_to_decimal(binary);
			int freq_index = rand() % MAX_FREQS;
			float freq_factor = pow(2,(FREQ_START_INDEX-freq_index));
			printf("%d %d %d %d %lf\n", i, freq_index, decimal, tam_testtime[i].tam_list[decimal] , tam_testtime[i].peak_power * freq_factor);
			if(tam_testtime[i].tam_list[decimal] * freq_factor < TAM_WIDTH_MAX && tam_testtime[i].peak_power * freq_factor < power_budget)
			{
				object->info[i+BIT_LENGTH*SIZE] = freq_index;
				i++;
			}
		}*/


}

int particle_swarm_optimization()
{
	
	time_t seconds;
	time( & seconds ) ;
	long int gencount=0;
	int temp=0;
	long int temptime=LARGENUMBER;
	long int lasttesttime = LARGENUMBER;
	int samegencount=0;
	int generation = 0;
	cout << "initialiseparticle " << endl;
	for (int i=0;i<( NO_OF_PARTICLES  ); ){
			//printf("not ok\n");
    		initialiseparticle(&(partarray[i]));
    		if(i%50 == 0){
    			 	 		cout << "particle no " << i << endl;

            }
    		long int temp = bin_packing(partarray[i].info,i,generation);
    		if(i != 0 && temp != -1)
    		{
	    			partarray[i].time_fitness = temp;
	    			//printf("ok\n");
	    		
				for (int ii=0;ii< NEW_PARTICLE_SIZE;ii++)
		    	{
		    	        partarray[i].lbest[ii]=partarray[i].info[ii];
		   	 	}
				partarray[i].lbesttime=partarray[i].time_fitness;

			for ( int l = 0 ; l < i ; l ++ ) {
				int flag = 0 ;
				for ( int t = 0 ; t < NEW_PARTICLE_SIZE && flag == 0; t ++ ) {
					if ( partarray [ l ] . info [ t ] != partarray [ i ] . info [ t ] ) {
						flag = 1 ;
						
					}
				}
				if ( flag == 0 ) break ;
				if ( l == i - 1 && partarray[i].time_fitness != -1) {
					//printf("Generated %d \n", i);
					// print_array(NEW_PARTICLE_SIZE, partarray[i].info, '\n');
					i ++ ;
					break ;
				}
			}
		}
		if ( i == 0 && temp != -1) 
			{
				partarray[i].time_fitness = temp;
				//printf("Generated... %d \n", i);
				
				for (int ii=0;ii< NEW_PARTICLE_SIZE;ii++)
		    	{
		    	        partarray[i].lbest[ii]=partarray[i].info[ii];
		   	 	}
				partarray[i].lbesttime=partarray[i].time_fitness;
				i ++;
			}
		
	}
	// printf("test\n");
	// for (int i=0;i<NO_OF_PARTICLES;i++)
	//for(int kkk = 0; kkk< NO_OF_PARTICLES; kkk++)
        	//printf("%ld\t",partarray[kkk].time_fitness);

    particle_test = 0;
        // printf("\n\n");
	// int generation = 0;
	/*
	for (int i=0;i<NO_OF_PARTICLES;i++)
    	{
    		// printf("test\n");
    		partarray[i].time_fitness = bin_packing(partarray[i].info,i,generation);  //testtime returns maxtime for scheduling process
    		
    		for (int ii=0;ii< BIT_LENGTH * SIZE;ii++)
        	{
        	        partarray[i].lbest[ii]=partarray[i].info[ii];
       	 	}
    		partarray[i].lbesttime=partarray[i].time_fitness;
    	}
    */
    	int flag = 0;
    	for (int i=0;i<NO_OF_PARTICLES;i++)
    	{
      	 	if(partarray[i].time_fitness<temptime)
         	{
                	//Make it the globalbest
                	temp=i;
                	globalbestIndx = i;
                	temptime=partarray[i].time_fitness;
                	flag = 1;
            }
    	}
    	//printf("global best********%d\n",globalbestIndx);
   	 if(flag == 1)
   	 {
    		globalbest.time_fitness=temptime;
    		// Copy info of the best particle to the global best particle
    		for (int i=0;i< NEW_PARTICLE_SIZE;i++)
    			globalbest.info[i]=partarray[globalbestIndx].info[i];
    		for ( int pp = 0 ; pp < SIZE ; pp ++ )
        	{
       				BestParticle [ pp ] .corenum = scheduler [globalbestIndx][pp]. corenum;
          			BestParticle [ pp ] .tam_width = scheduler [globalbestIndx][pp]. tam_width;
               		BestParticle [ pp ] .starttime = scheduler [globalbestIndx][pp]. starttime;
                	BestParticle [ pp ] .endtime = scheduler [globalbestIndx][pp]. endtime;
                	BestParticle [ pp ] .tsv = scheduler [globalbestIndx][pp]. tsv;
                	BestParticle [ pp ] .freq = scheduler [globalbestIndx][pp]. freq;
                	//printf ("{ %d, %d, %d, %d }\n", scheduler [globalbestIndx][pp]. io, scheduler [globalbestIndx][pp]. corenum, scheduler [globalbestIndx][pp]. starttime, scheduler [globalbestIndx][pp]. endtime);
        	}
    	}
    	int count;
    	while(gencount<MAX_ITERATION)
    	{
        	cout << "current generation " << gencount << endl;
        	for(int i=0;i<NO_OF_PARTICLES;i++)
        	{
            	   	 /*for(int kkk = 0; kkk< SIZE; kkk++)
            	   	 	printf("%lf\t",partarray[i].info[kkk]);
            	   	 printf("\n\n");*/
            	   	 	if(i%100 == 0){
            	   	 		cout << "particle no " << i << endl;

            	   	 	}
            	   	 //swaptamindex(partarray[i].info,partarray[i].lbest,globalbest.info, ALPHA,BETA); //paticle's IOpair part is updated to pbest  with prob APLHA	
            		int temp_info[SIZE];
            		for(int k = 0; k < SIZE; k++)
            		{
            			temp_info[k] = partarray[i].info[k];
            		}
            		swapparticleinfo(partarray[i].info,partarray[i].lbest,ALPHA);
            		swapparticleinfo(partarray[i].info,globalbest.info,BETA);
            		/*for(int kkk = 0; kkk< SIZE; kkk++)
            	   	 	printf("%lf\t",partarray[i].info[kkk]);
            	   	 printf("\n\n");*/
            		//swapiopair(partarray[i].info,globalbest.info, BETA);     //paticle's IOpair part is updated to gbest  with prob BETA
    			
	    			// partarray[i].time_fitness = bin_packing(partarray[i].info,i,gencount);
	    			long int time_fitness = bin_packing(partarray[i].info,i,gencount);
	    			if(time_fitness == -1)
	    			{
	    				for(int k = 0; k < SIZE; k++)
	    				{
	    					partarray[i].info[k] = temp_info[k] ;
	    				}
	    			}
	    			else
	    			{
	    				partarray[i].time_fitness = time_fitness;
	    			
	    			if((partarray[i].time_fitness)<=(partarray[i].lbesttime))
	            		{    
	                		for (int ii=0;ii< NEW_PARTICLE_SIZE;ii++)
	                		{
	                			partarray[i].lbest[ii]=partarray[i].info[ii];
	                		}
	                		partarray[i].lbesttime=partarray[i].time_fitness;
	            		}
	            	}	
            }
               	int flag = 0;
        	for (int i=0;i<NO_OF_PARTICLES;i++)
        	{
           		 if(partarray[i].time_fitness<temptime) //temptime is holding globalbest time for a generation.
               		 {
                    		//Make it the globalbest
                   			temp=i;
                    		globalbestIndx = i;
                    		temptime=partarray[i].time_fitness;
                    		flag = 1;
                    
                	}
        	}
        	//printf("global best********%d\n",globalbestIndx);
        	if (flag == 1)
   			{
      			globalbest.time_fitness=temptime;
        		for (int i=0;i< NEW_PARTICLE_SIZE;i++)
        			globalbest.info[i]=partarray[globalbestIndx].info[i];
       			for ( int pp = 0 ; pp < SIZE ; pp ++ )
       			{
						BestParticle [ pp ] .corenum = scheduler [globalbestIndx][pp]. corenum;
               		 	BestParticle [ pp ] .tam_width = scheduler [globalbestIndx][pp]. tam_width;
                		BestParticle [ pp ] .starttime = scheduler [globalbestIndx][pp]. starttime;
                		BestParticle [ pp ] .endtime = scheduler [globalbestIndx][pp]. endtime;
                		BestParticle [ pp ] .tsv = scheduler [globalbestIndx][pp]. tsv;
                		BestParticle [ pp ] .freq = scheduler [globalbestIndx][pp]. freq;
                		//printf ("{ %d, %d, %d, %d }\n", scheduler [globalbestIndx][pp]. io, scheduler [globalbestIndx][pp]. corenum, scheduler [globalbestIndx][pp]. starttime, scheduler [globalbestIndx][pp]. endtime);
                }
   			}	
   	
   		gencount++;
   		//printf("\ngencount is %ld\n", gencount);
        if (temptime == lasttesttime)
   		{
   			samegencount++;
   			printf("\nlasttime is %ld for gen %ld",lasttesttime, gencount);
   			printf("\n\n the semgencount value is %d\n", samegencount);
   		}
   		else
   		{	
   			samegencount = 0;
   			lasttesttime= temptime;
   		}	
   		if (samegencount == 150)
   			break;	
   	}
   	printf("\ngencount is %ld\n", gencount);

    	return temptime;		
}


void print_parallel(long int mint, schedulerInfo BestParticle [])
{
	int breakpoint[SIZE*2];
	init_array(SIZE*2, breakpoint);
	int count = 1;
	int tam[SIZE];
	int level[SIZE];
    int finallayer[NDies];

    for ( int i = 1 ; i < SIZE + 1; i ++ ) 
    {
    	tam[BestParticle[i-1].corenum] = BestParticle[i-1].tam_width;
    	level[BestParticle[i-1].corenum] = tam_testtime[BestParticle[i-1].corenum].die_level;
    	finallayer[tam_testtime[BestParticle[i-1].corenum].die_level] = BestParticle[i-1].tsv;

    	int bp = BestParticle[i-1].endtime; 
    	int j = count;
    	while(breakpoint[j-1] > bp)
    	{
    		breakpoint[j] = breakpoint[j-1];
    		j--;
    	} 
		breakpoint[j] = bp;
		count++;
    	


    	
    }
    
    // print_array(count, breakpoint, '\n');
    for ( int i = 0 ; i < count ; i ++ ) 
    {
    	int parallel[NDies];
    	int dietam[NDies];
    	int layer[NDies-1];
    	int parallel_cores[SIZE];
    	init_array(NDies, parallel);
    	init_array(SIZE, parallel_cores);
    	init_array(NDies, dietam);
    	for(int j = 0; j < SIZE; j++)
    	{
	    	if(BestParticle[j].starttime <= breakpoint[i] && breakpoint[i] < BestParticle[j].endtime)
	    	{
	    		parallel[tam_testtime[BestParticle[j].corenum].die_level] = 1;
	    		dietam[tam_testtime[BestParticle[j].corenum].die_level] += tam[BestParticle[j].corenum];
	    		
	    		parallel_cores[BestParticle[j].corenum] = 1;

	    	}
	    }
	    printf ("%d ;",breakpoint[i]);
    	for(int j = 0; j < NDies; j++)
    	{
    		if(parallel[j])
				printf ("%d ", j + 1);
    	}

	    printf(";");
	    for(int j = 0; j < SIZE; j++)
    	{
    		if(parallel_cores[j])
				printf ("%d ", j + 1);
    	}
	    printf(";");

       
        print_array(NDies, dietam, ';');
    	int tsv = tsv_required(parallel_cores, tam, layer, level);
    	print_array(NDies-1, layer, ';');
    	print_array(NDies, finallayer, ';');

	    printf("\n");


    }
    printf("\n");
}

void print_csv(long int mint, schedulerInfo BestParticle [])
{
	/*Print Parallelism, TAM width and TSV layer in CSV*/
	printf("%d;%ld;",TAM_WIDTH_MAX,mint);
    int prev = 0;
    for ( int i = 0 ; i < SIZE ; i ++ ) {
    	if(i == 0)
    		printf ("%d", BestParticle[i].corenum+1);
    	else if(BestParticle[i]. starttime == prev)
            printf ("|%d", BestParticle[i].corenum+1);
        else
            printf (",%d", BestParticle[i].corenum+1);
        prev = BestParticle[i].starttime;
    }
	printf(";");
    for ( int i = 0 ; i < SIZE ; i ++ ) {
    	for ( int j = 0 ; j < SIZE ; j ++ ) {
    		if(BestParticle[j].corenum == i)
    		{
    			printf("%d",BestParticle[j].tam_width);
    			if(i != SIZE-1)
					printf(",");
			}        			
    	}
    }
	printf(";");
	for ( int i = 0 ; i < SIZE ; i ++ ) {
    	for ( int j = 0 ; j < SIZE ; j ++ ) {
    		if(BestParticle[j].corenum == i)
    		{
    			printf("%d",BestParticle[j].tsv);
        		if(i != SIZE-1)
					printf(",");
    		}
    	}
    }
	printf(";\n");
}




int main(int argc, char *argv [ ])
{
	long int ppp=1111990911;
	
	
	time_t seconds;
	time( & seconds ) ;
	srand ( seconds ) ;
	srand48 ( seconds ) ; //srand generates a sequence of random numbers depending upon the seed value
	FILE *p1;

	//p1 = fopen ( argv [ 1 ],"w" ) ;
	int NoIteration = atoi ( argv [ 1 ] ) ;
	power_budget = atof (argv[2]);
        //TAM_WIDTH_MAX = atoi ( argv [ 3 ] ) ;
        int iteration;

	long int t,mint;
	
	FILE *tamread;
	tamread = fopen("tam_testtime.txt","r");
	int kk = 0;
	 cout << power_budget << "is the power_budget" << endl;

    cout << TSV_MAX << " is TSV_MAX " << endl;

    cout << TAM_WIDTH_MAX << "is TAM_WIDTH_MAX " << endl;
	for(int ii = 0; ii< NDies; ii++)
	{
		for(int jj = 0; jj< Cores[ii]; jj++)
		{
			tam_testtime[kk].no_of_tam =  TAM_INFO[ii][jj];
			tam_testtime[kk].die_level =  ii;
			kk = kk + 1;
		}
	}
	for(int ii = 0; ii< SIZE; ii++)
		tam_testtime[ii].peak_power = FIXED_POWER[ii];	

	for(int ii = 0; ii< SIZE; ii++)
	{
		for(int jj = 0; jj< tam_testtime[ii].no_of_tam; jj++)
		{
			fscanf(tamread,"%d\t%ld", &tam_testtime[ii].tam_list[jj], &tam_testtime[ii].testtime_list[jj]); 
		}
	}
	fclose (tamread);
	// printf("test\n");
	/*temp_tam_index[0] = TAM_INFO[0];
	for(int i = 1; i< SIZE; i++)
	{
		temp_tam_index[i] = temp_tam_index[i-1]+TAM_INFO[i];
	}
	*/

	if(HARD_DIE_TEST)
	{
		initialiseparticle(&(partarray[0]));
		partarray[0].time_fitness = bin_packing(partarray[0].info,0,0);
		// printf("%ld \n", partarray[0].time_fitness);
		for ( int i = 0 ; i < SIZE ; i ++ ) {
		        //printf ("{ %d, %d, %ld, %ld, %d }\n", scheduler[0][i].corenum, scheduler[0][i].tam_width, scheduler[0][i]. starttime, scheduler[0][i]. endtime,  scheduler[0][i]. tsv);
		}
		print_csv(partarray[0].time_fitness, scheduler[0]);
		return 0;
	}
	
	////////////////////////////////////////////////////

	// modified by nil to check
		
	cout << NDies << " " <<"nil is here" << endl;
	kk = 0;
	cout << "no of tam " << " " << " die level" <<" " << "peak_power" << endl;
 	 for(int ii = 0; ii< NDies; ii++)
	{
		for(int jj = 0; jj< Cores[ii]; jj++)
		{
			cout << " ************* " << tam_testtime[kk].no_of_tam << " " ;
			cout << tam_testtime[kk].die_level << " ";
			cout << tam_testtime[kk].peak_power << endl;
			for(int jj = 0; jj< tam_testtime[kk].no_of_tam; jj++)
			{
				cout << tam_testtime[kk].tam_list[jj] << " " << tam_testtime[kk].testtime_list[jj] << endl; 
			}
			kk = kk + 1;
		}
	}

	//end modified by nil

	/////////////////////////////////////////////

	for(iteration=0;iteration < NoIteration;iteration++)
	{
		globalbestIndx  = 0;
		// printf("test\n");
		mint=particle_swarm_optimization();
		//printf ("Globalbestindx: %d\n",globalbestIndx);
		/*fo/r ( int pp = 0 ; pp < SIZE ; pp ++ ){
                 	       BestParticle [ pp ] .tam_width = scheduler [globalbestIndx][pp]. tam_width;
                        BestParticle [ pp ] .corenum = scheduler [globalbestIndx][pp]. corenum;
                        BestParticle [ pp ] .starttime = scheduler [globalbestIndx][pp]. starttime;
                        BestParticle [ pp ] .endtime = scheduler [globalbestIndx][pp]. endtime;
                }*/
		for(int i=1;i<1;i++)
		{
    			   printf("RunForrestRun: %d\n",i);
    				globalbestIndx = 0;
    			
    			t=particle_swarm_optimization();
    			/*for ( int ii = 0 ; ii < SIZE ; ii ++ ) {
                printf ("{ %d, %d, %ld, %ld }\n", scheduler [globalbestIndx][ii].corenum, scheduler [globalbestIndx][ii].tam_width, scheduler [globalbestIndx][ii]. starttime, scheduler [globalbestIndx][ii]. endtime);
        }*/

			if ( t < mint ) {
				mint = t ;
				/*for ( int pp = 0 ; pp < SIZE ; pp ++ ){
                                        BestParticle [ pp ] .corenum = scheduler [globalbestIndx][pp]. corenum;
                                        BestParticle [ pp ] .tam_width = scheduler [globalbestIndx][pp]. tam_width;
                                        BestParticle [ pp ] .starttime = scheduler [globalbestIndx][pp]. starttime;
                                        BestParticle [ pp ] .endtime = scheduler [globalbestIndx][pp]. endtime;
                                }*/
			}
			
			printf("%ld \n",t);
		}
	
		// fprintf(p1,"%ld\n",mint);
		 printf(" best time is %ld\n",mint);
		
	}
        /*for ( int i = 0 ; i < SIZE ; i ++ ) {
                printf ("{ %d, %d, %ld, %ld, %d , %lf}\n", BestParticle[i].corenum, BestParticle[i].tam_width, BestParticle[i]. starttime, BestParticle[i]. endtime,  BestParticle[i]. tsv, BestParticle[i]. freq);
        }*/

        cout << "BestParticle[i].corenum, BestParticle[i].tam_width, BestParticle[i]. starttime, BestParticle[i]. endtime,  BestParticle[i]. tsv, BestParticle[i]. freq" << endl;
        for ( int i = 0 ; i < SIZE ; i ++ ) {
                printf ("{ %d, %d, %ld, %ld, %d , %lf}\n", BestParticle[i].corenum, BestParticle[i].tam_width, BestParticle[i]. starttime, BestParticle[i]. endtime,  BestParticle[i]. tsv, BestParticle[i]. freq);
        }

    // print_csv(mint, BestParticle);
    //print_parallel(mint, BestParticle);
	

    return 0;
}


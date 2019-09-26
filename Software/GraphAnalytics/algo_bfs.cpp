#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <time.h>
#include <omp.h>
#include <stdint.h>
#include <iostream>

#define  final_run 							1024
#define  control_1							256
#define  control_2							512
#define  call_counter_mask			1023
#define  FPGA_no_of_writes_mask		524287
#define 	FILE_NAME  					"soc.txt"
#define 	V  		 							92436*2
#define 	E			  							4048860
#define 	I			  							(256*1024)
#define 	P			  							10
#define 	CL_PER_INTERBAL			4096
#define	CL_ID_MASK					15
#define	MaxUpdateBufferSize	2520700
#define 		ThreadNum						16
#define		LocalUpdateBufferSize		1024
#define		workload_of_FPGA		3  // fpga start work from
#define		r_threshold					0.02

#ifndef CL
# define CL(x)                     ((x) * 64)
#endif // CL
#ifndef LOG2_CL
# define LOG2_CL                   6
#endif // LOG2_CL
#ifndef MB
# define MB(x)                     ((x) * 1024 * 1024)
#endif // MB

#define LPBK1_DSM_SIZE           MB(4)


///////////////////////////////////////////

typedef int32_t bt32bitInt;
typedef uint32_t btUnsigned32bitInt;
typedef uint64_t btUnsigned64bitInt;
typedef uint64_t btWSSize;     ///< Workspace size type.

typedef unsigned char UCHAR;
typedef UCHAR *PUCHAR;
typedef PUCHAR btVirtAddr;

struct VAFU2_CNTXT {
   union {
      btUnsigned64bitInt          qword0[8];       // make it a whole cacheline
      struct {
         union {                                   // first qword
            btUnsigned64bitInt    dword0;
            struct {
               btUnsigned64bitInt rsvd0:   32;
               btUnsigned64bitInt delay:   16;     // undefined, but in the structure definition, set to 0
            };
         };
         void          *pSource;                   ///< Source user-mode virtual address, cache-line aligned.
         void          *pDest;                     ///< Destination user-mode virtual address, cache-line aligned.
         btUnsigned32bitInt       num_cl;          ///< Number of cache lines to copy.
      };
   };
   union {
      btUnsigned64bitInt          qword1[8];       // make it a whole cacheline
      struct {
         btUnsigned32bitInt       Status;          ///< Bit0 true if done
#define  VAFU2_CNTXT_STATUS_DONE   0x00000001      ///< Bit0 selector
      };
   };
}; // struct VAFU2_CNTXT

////////////////////////////////////////////////////////////

typedef struct {
	  unsigned char interval_id;
	  unsigned int 	edge_start_cl;
	  unsigned int 	edge_end_cl;
	  unsigned int  edge_start_offset;
	  unsigned int	edge_end_offset;
	  unsigned int 	no_of_active_vertex;
	  unsigned int * update_buffer;
	  unsigned int  update_buffer_counter;
} interval;

typedef struct {  
	  unsigned int 	offset;
	  unsigned int 	count;
} vertex;


void bfs(){
	//=============================
	// Now we have the NLB Service
	//   now we can use it
	//=============================

	// btVirtAddr         pWSUsrVirt = m_pWkspcVirt; // Address of Workspace
	// const btWSSize     WSLen      = m_WkspcSize; // Length of workspace in bytes

	int i, j, k, o, p, q;

	// Number of bytes in each of the source and destination buffers (4 MiB in this case)
	btUnsigned32bitInt a_num_bytes= 1024*1024*8; //8MiB guessed
	btUnsigned32bitInt a_num_cl   = a_num_bytes / CL(1);  // number of cache lines in buffer

	// VAFU Context is at the beginning of the buffer
	btVirtAddr 			pWSUsrVirt = (btVirtAddr) malloc((2*627700 * CL(1) + 40960 * CL(1))*sizeof(UCHAR)); // Char Memory : Large Buffer
	VAFU2_CNTXT       	*pVAFU2_cntxt = reinterpret_cast<VAFU2_CNTXT *>(pWSUsrVirt);
	btVirtAddr         	pSource_V = pWSUsrVirt;

	btVirtAddr         pSource_E = pSource_V + 40960 * CL(1);

	// The destination buffer is right after the source buffer
	btVirtAddr         pDest   = pSource_E + 627700* CL(1);

	omp_set_nested(1);

	int						*LocalUpdateBufferCounter[ThreadNum]; 
	int 						**UpdateBufferCache[ThreadNum];
	
	for(i=0; i<ThreadNum; i++){
		LocalUpdateBufferCounter[i] = (int*) malloc(sizeof(int)*P);	
		UpdateBufferCache[i] = (int**) malloc(sizeof(int*)*P);
		for(j=0;	j<P; j++){
			LocalUpdateBufferCounter[i][j]=0;
			UpdateBufferCache[i][j] = (int*) malloc(sizeof(int)*LocalUpdateBufferSize);	
		}
	}

      bt32bitInt delay(0.001);   
      
	  int test_source = 2;	
	  int root = 2;
	  int current_level =0;	
	  int call_counter = 15;
	  int interval_id;
	  int FrontierSize = 0;
	  int FPGA_no_of_writes;
	  int counter=0;
	  int have_update=1;
	  volatile bt32bitInt done = 0;

  	// read files and initialization
	interval* 		Intervals 	= (interval*) malloc(sizeof(interval)*10);		
	vertex* 			vertex_set 	= (vertex*) malloc(sizeof(vertex)*P*I);   // this set is only used by cpu

	for(i=0; i<P*I; i++){
		vertex_set[i].offset = 0;
		vertex_set[i].count = 0;
	}

	   ::memset( pSource_V, 	0xff,  	40960*CL(1) );
      ::memset( pDest,   			0x00, 	320000*CL(1) );      	         
      *((unsigned char*)(pSource_V)+root) = current_level;
	  	for(i=0; i<40960; i++){	  		 		 
		 	*((unsigned short*)(pSource_V)+31+32*i) = (*((unsigned short*)(pSource_V)+31+32*i) & CL_ID_MASK) + ((i&4095)<<4);
		} 
	  
	 	 FILE *fp;  unsigned int u1, u2;

		if ((fp=fopen(FILE_NAME,"r"))==NULL) printf("Cannot open file. Check the name.\n"); 
		else {
			if(fscanf(fp,"%d %d\n",&u1,&u2)!=EOF){				
					*((unsigned int*)(pSource_E)+1) = u1;
					*((unsigned int*)(pSource_E)+0) = u2;					
			}
			for(i=1; i<E; i++){
				if(fscanf(fp,"%d %d\n",&u1,&u2)!=EOF){
					*((unsigned int*)(pSource_E)+i*2+1)=u1;			// source 
					*((unsigned int*)(pSource_E)+i*2+0)=u2;         // destination 																																				
				}
				if(*((unsigned int*)(pSource_E)+i*2+1) != *((unsigned int*)(pSource_E)+(i-1)*2+1)){
					vertex_set[*((unsigned int*)(pSource_E)+i*2+1)].offset = i;	
					vertex_set[*((unsigned int*)(pSource_E)+(i-1)*2+1)].count=counter+1;
					counter=0;
				} else{
					counter++;
				}
				// printf("Edges: %d %d\n",u1,u2);
			}
			fclose(fp);
		}
		
		for(i=0; i<P; i++){
			for(j=I*i; j<(i+1)*I-1; j++){
				if(vertex_set[j].count !=0) {
						Intervals[i].edge_end_offset = vertex_set[j].offset+vertex_set[j].count;
						Intervals[i].edge_end_cl = Intervals[i].edge_end_offset/8;		
				}		
			}	
			for(j=(i+1)*I-1; j>i*I; j--){
				if(vertex_set[j].count !=0) {
						Intervals[i].edge_start_offset = vertex_set[j].offset;
						Intervals[i].edge_start_cl = Intervals[i].edge_start_offset/8;	
				}		
			}	
			Intervals[i].interval_id = i;
			Intervals[i].no_of_active_vertex = 0;
			Intervals[i].update_buffer_counter = 0;
			Intervals[i].update_buffer = (unsigned int *) malloc(sizeof(unsigned int)*MaxUpdateBufferSize);
			for(j=0; j<MaxUpdateBufferSize; j++)
				Intervals[i].update_buffer[j] = 0;
		}		
		Intervals[0].no_of_active_vertex=1;
		
		//for(i=0; i<P; i++){cout<<Intervals[i].edge_start_cl << " "<<Intervals[i].edge_end_cl<<" "<<Intervals[i].edge_start_offset<<" "<<Intervals[i].edge_end_offset<<endl;}
	
		//------------------------------------------------------------------
		
	  ::memset(pVAFU2_cntxt, 0, sizeof(VAFU2_CNTXT));
      //pVAFU2_cntxt->num_cl  = 4096;
      //pVAFU2_cntxt->pSource = pSource_V;
      //pVAFU2_cntxt->pDest   = pDest;
	  pVAFU2_cntxt->dword0  = ((call_counter-1)<<11);	  		  	 
	  struct timespec start, stop; 
	  double exe_time;	  
	  double total_time = 0;
	  int thread_id;
	  int edge_centric_no = 0;
	  while(have_update){

	 		have_update = 0;
	 		// if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
	 		edge_centric_no = 0;
	 		for(p=0; p<P; p++){
	 			if(Intervals[p].no_of_active_vertex>r_threshold*I)
	 					edge_centric_no++;
	 		}		 			
		 //----------------------------scatter------------------------------
		if(edge_centric_no<P){				
			for(p=0; p<P; p++){
				if(Intervals[p].no_of_active_vertex	!= 0){
					if(Intervals[p].no_of_active_vertex>r_threshold*I){  // edge centric goes here
						std::cout<<"edge centric"<<std::endl;
						#pragma omp parallel num_threads(ThreadNum) shared(vertex_set, pSource_V, pSource_E, Intervals, p) private(i, j, k, interval_id, thread_id)
						{
							thread_id = omp_get_thread_num();					
							#pragma omp for schedule(static) 					
							for(k=Intervals[p].edge_start_offset; k<Intervals[p].edge_end_offset; k++){
								if((*((unsigned char*)(pSource_V)+*((unsigned int*)(pSource_E)+k*2+1))==current_level) & (*((unsigned int*)(pSource_E)+k*2+1)%64<61)){											
									interval_id=*((unsigned int*)(pSource_E)+k*2+0)/I;	
									UpdateBufferCache[thread_id][interval_id][LocalUpdateBufferCounter[thread_id][interval_id]] = *((unsigned int*)(pSource_E)+k*2+0);
									LocalUpdateBufferCounter[thread_id][interval_id]++;								
									if(LocalUpdateBufferCounter[thread_id][interval_id] == LocalUpdateBufferSize){															
										#pragma omp critical 
										{
											for(j=Intervals[interval_id].update_buffer_counter; j<Intervals[interval_id].update_buffer_counter+LocalUpdateBufferSize; j++){											
												 Intervals[interval_id]. update_buffer[j]	 = UpdateBufferCache[thread_id][interval_id][j-Intervals[interval_id].update_buffer_counter];
											}																					
											Intervals[interval_id].update_buffer_counter+=LocalUpdateBufferSize;
										}
										LocalUpdateBufferCounter[thread_id][interval_id]=0;														
									}	
								}	
							}												
						}					
					} else {   // vertex centric goes here
						std::cout<<"vertex centric"<<std::endl;
						#pragma omp parallel num_threads(ThreadNum) shared(vertex_set, pSource_V, pSource_E, Intervals, p) private(i, j, k, interval_id, thread_id)
						{
							thread_id = omp_get_thread_num();					
							#pragma omp for schedule(static, I/ThreadNum) 
							for(i=p*I; i<(p+1)*I-1; i++){
								if((*((unsigned char*)(pSource_V)+i)==current_level) & (i%64<61)){								
									for(k=vertex_set[i].offset; k<vertex_set[i].offset+vertex_set[i].count; k++){
										interval_id=*((unsigned int*)(pSource_E)+k*2+0)/I;	
										UpdateBufferCache[thread_id][interval_id][LocalUpdateBufferCounter[thread_id][interval_id]] = *((unsigned int*)(pSource_E)+k*2+0);
										LocalUpdateBufferCounter[thread_id][interval_id]++;								
										if(LocalUpdateBufferCounter[thread_id][interval_id] == LocalUpdateBufferSize){															
											#pragma omp critical 
											{
												for(j=Intervals[interval_id].update_buffer_counter; j<Intervals[interval_id].update_buffer_counter+LocalUpdateBufferSize; j++){											
												 	 Intervals[interval_id]. update_buffer[j]	 = UpdateBufferCache[thread_id][interval_id][j-Intervals[interval_id].update_buffer_counter];
												}																					
												Intervals[interval_id].update_buffer_counter+=LocalUpdateBufferSize;
											}
											LocalUpdateBufferCounter[thread_id][interval_id]=0;														
										}
									}
								}	
							}
						}										
					}
				}
			}		
		} else {
					#pragma omp parallel num_threads(2) shared(Intervals, pSource_V, pSource_E, vertex_set, current_level) private(i, j, k, p,q, interval_id) 
					{	
						#pragma omp sections private(i, p, q, interval_id)
						{  
							#pragma omp section
							{		
								for(p=workload_of_FPGA; p<10; p++){ 
									  pVAFU2_cntxt->num_cl  = 4096;
							      	  pVAFU2_cntxt->pSource = pSource_V+p*4096*64;
									  pVAFU2_cntxt->dword0  = (call_counter<<11)+control_1+current_level;
								       done = ((pVAFU2_cntxt->Status)>>1)  & call_counter_mask;
								      while ((done !=call_counter)) {
								         done = ((pVAFU2_cntxt->Status)>>1)  & call_counter_mask;
								      }    
								      call_counter++;           
								     
									  // read edges	  				  			      
								      pVAFU2_cntxt->pDest     = pDest;
								      pVAFU2_cntxt->pSource = pSource_E+Intervals[p].edge_start_cl*64;				   
									  pVAFU2_cntxt->num_cl  = Intervals[p].edge_end_cl-Intervals[p].edge_start_cl+20;
								      pVAFU2_cntxt->dword0  = (call_counter<<11)+control_2+current_level;
								      
									  //cout <<  pVAFU2_cntxt->num_cl <<endl;     
								       done = ((pVAFU2_cntxt->Status)>>1)  & call_counter_mask;
								      while ((done!=call_counter)) {
								         done = ((pVAFU2_cntxt->Status)>>1)  & call_counter_mask;
								      }            
								      call_counter++;
								     
								      FPGA_no_of_writes = ((pVAFU2_cntxt->Status)>>11)  & FPGA_no_of_writes_mask;  		
										
										#pragma omp critical 
										{
									      if(FPGA_no_of_writes !=0){
												for(q=0; q<16*FPGA_no_of_writes; q++){
													interval_id =  *((unsigned int*)(pDest)+q)/I ;
													Intervals[interval_id].update_buffer[Intervals[interval_id].update_buffer_counter] = *((unsigned int*)(pDest)+q);
													Intervals[interval_id].update_buffer_counter++;	
												}	
											}
										}	
							   }
							}
							#pragma omp section
							{								
								for(p=0; p<workload_of_FPGA; p++){
									#pragma omp parallel num_threads(ThreadNum) shared(vertex_set, pSource_V, pSource_E, Intervals, p) private(i, j, k, interval_id, thread_id)
									{
										thread_id = omp_get_thread_num();					
										#pragma omp for schedule(static) 					
										for(k=Intervals[p].edge_start_offset; k<Intervals[p].edge_end_offset; k++){
											if((*((unsigned char*)(pSource_V)+*((unsigned int*)(pSource_E)+k*2+1))==current_level) & (*((unsigned int*)(pSource_E)+k*2+1)%64<61)){												
												interval_id=*((unsigned int*)(pSource_E)+k*2+0)/I;	
												UpdateBufferCache[thread_id][interval_id][LocalUpdateBufferCounter[thread_id][interval_id]] = *((unsigned int*)(pSource_E)+k*2+0);
												LocalUpdateBufferCounter[thread_id][interval_id]++;								
												if(LocalUpdateBufferCounter[thread_id][interval_id] == LocalUpdateBufferSize){															
													#pragma omp critical 
													{
														for(j=Intervals[interval_id].update_buffer_counter; j<Intervals[interval_id].update_buffer_counter+LocalUpdateBufferSize; j++){											
															 Intervals[interval_id]. update_buffer[j]	 = UpdateBufferCache[thread_id][interval_id][j-Intervals[interval_id].update_buffer_counter];
														}																					
														Intervals[interval_id].update_buffer_counter+=LocalUpdateBufferSize;
													}
													LocalUpdateBufferCounter[thread_id][interval_id]=0;														
												}	
											}	
										}												
									}
								}												
							}
					}
			}				
		}
		edge_centric_no = 0;
		
		// Serial flush	buffer cache					
			
			for(thread_id=0; thread_id<ThreadNum; thread_id++){
				for(interval_id=0;	interval_id<P; interval_id++){			
					if(LocalUpdateBufferCounter[thread_id][interval_id]!=0 && LocalUpdateBufferCounter[thread_id][interval_id]<LocalUpdateBufferSize){											
						for(i=Intervals[interval_id].update_buffer_counter; i<Intervals[interval_id].update_buffer_counter+LocalUpdateBufferCounter[thread_id][interval_id]; i++){
						 	Intervals[interval_id].update_buffer[i] = 	UpdateBufferCache[thread_id][interval_id][i-Intervals[interval_id].update_buffer_counter];
						}
						Intervals[interval_id].update_buffer_counter+= LocalUpdateBufferCounter[thread_id][interval_id];
					}											
					LocalUpdateBufferCounter[thread_id][interval_id]=0;
				}
			}
			
			
		for(p=0; p<P; p++){Intervals[p].no_of_active_vertex=0;}			
		 //----------------------------gather------------------------------		
			#pragma omp parallel for num_threads(ThreadNum)  shared(Intervals, pSource_V) private(p,j) schedule(dynamic) reduction(+:FrontierSize)	   		
			for(p=0; p<P; p++){					
				if(Intervals[p].update_buffer_counter!=0){					
					for(j=0; j<Intervals[p].update_buffer_counter; j++){
						if(Intervals[p].update_buffer[j]%64 <61){
							if(*((unsigned char*)(pSource_V)+Intervals[p].update_buffer[j])  > current_level+1){
								*((unsigned char*)(pSource_V)+Intervals[p].update_buffer[j]) = current_level+1;							
								Intervals[p].no_of_active_vertex ++;
								FrontierSize++;								
							}
						}
					}	
				}
				Intervals[p].update_buffer_counter=0;									
			}			
			
			current_level++;												
			
			have_update = FrontierSize; 
					 
			
			//if( clock_gettime(CLOCK_REALTIME, &stop) == -1) { perror("clock gettime");}
			//exe_time = (stop.tv_sec - start.tv_sec)+ (double)(stop.tv_nsec - start.tv_nsec)/1e9;				
			//total_time = total_time + exe_time;
			FrontierSize=0;	
	} 					
	//printf("total time is  %f sec\n", total_time);



}

int main() {
	bfs();
	return 0;
}
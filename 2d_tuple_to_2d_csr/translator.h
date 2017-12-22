#include <mpi.h>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <cstring>
#include <unistd.h>

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

typedef struct packed_edge {
  int64_t v0;
  int64_t v1;
} packed_edge;


#ifndef TIME_H
#define TIME_H

inline double wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1.e-6;

}
#endif


typedef uint64_t vertex_t;
typedef uint64_t index_t;
typedef uint64_t data_t;
off_t fsize(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1;
}

int compare(const void *a, const void *b)
{
	if(*(vertex_t*)a>*(vertex_t*)b) return 1;
	else return -1;
}

template<typename index_t>
void progress(
		index_t &report,
		index_t beg_time,
		int			tid,
		index_t prc_line
){
	if(prc_line>report && tid==0){
			std::cout<<prc_line<<" lines\t"
					<<"time: "<<wtime()-beg_time<<" seconds\n";
		report<<=1;
	}
}

std::string numtostr(index_t num)
{
	std::stringstream ss;
	ss.clear();
	ss.str("");
	ss<<num;
	return ss.str();
}

//using 3 64bits to store 4 48bits verts.
//the first 3 is stored lower bits
//the 4th is stored in the higher bits across 3 64bits.
//
#ifdef COMPRESS
void adj_store(index_t off, vertex_t vert, vertex_t *adj_list){
  index_t base_off=(off>>2)*3;
  index_t group_off=off&0x03;
  if(group_off<3)
	{
		/*clear your bits*/
		adj_list[base_off+group_off]&=0xffff000000000000;
		
		/*set your bits*/
		adj_list[base_off+group_off]|=vert;
	}
  else
	{//group_off==3
		/*clear your bits*/
		adj_list[base_off]&=0x0000ffffffffffff;
		adj_list[base_off+1]&=0x0000ffffffffffff;
		adj_list[base_off+2]&=0x0000ffffffffffff;
    
		/*set your bits*/
		adj_list[base_off]|=((vert&0xffff)<<48);
    adj_list[base_off+1]|=((vert&0xffff0000)<<32);
    adj_list[base_off+2]|=((vert&0xffff00000000)<<16);
  }
}

vertex_t adj_load(index_t off, const vertex_t *adj_list){
  vertex_t res;
  index_t base_off=(off>>2)*3;
  index_t group_off=off&0x03;
  if(group_off<3) res=adj_list[base_off+group_off]&0xffffffffffff;
  else{
    res=(adj_list[base_off]>>48)+
        ((adj_list[base_off+1]&0xffff000000000000)>>32)+
        ((adj_list[base_off+2]&0xffff000000000000)>>16);
  }
  return res; 
}
#else
void adj_store(index_t off, vertex_t vert, vertex_t *adj_list){
	adj_list[off]=vert;
}

vertex_t adj_load(index_t off, const vertex_t *adj_list){
  return adj_list[off]; 
}
#endif


template<typename data_t, typename index_t>
bool translator( 
				int 			argc,
				char**		argv,
				data_t*		&ranger_adj_list,
				index_t*	&ranger_beg_pos,
				index_t*	&ranger_adj_card)
{
  int size, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(rank==0) std::cout<<"Input: /path/to/exe scale degree row-partition column-partition generation-process-count\n";
	if(argc!=6){std::cout<<"Wrong input\n";exit(-1);}
	
	index_t scale=atoi(argv[1]);
	index_t degree=atoi(argv[2]);
	index_t row_par=atoi(argv[3]);
	index_t col_par=atoi(argv[4]);
	index_t gen_sz=atoi(argv[5]);
	
	index_t my_row=rank/col_par;
	index_t my_col=rank%col_par;

	if(size!=row_par*col_par)
	{
		std::cout<<"Wrong input. Each process responsible for one partition\n";
		exit(-1);
	}
	
	index_t vert_count=((int64_t)1)<<scale;
	index_t ranger_vert_count=vert_count/row_par;
	if(ranger_vert_count*row_par!=vert_count)
	{
		std::cout<<"Wrong partition, should be even!\n";
		exit(-1);
	}
	vertex_t my_vert_off=my_row*ranger_vert_count;
	ranger_beg_pos	= new index_t[ranger_vert_count+1];	
	ranger_adj_card	= new index_t[ranger_vert_count];	
	index_t *ranger_offset = new index_t[ranger_vert_count];
	memset(ranger_offset,0,sizeof(index_t)*ranger_vert_count);


	//get the total ranger_edge_count
	index_t ranger_edge_count=0;
	char filename[1024];
	for(index_t i=0;i<gen_sz;i++)
	{
			sprintf(filename,"scale-%lu-degree-%lu-rank-%lu-of-%lu-par-rowcol-%lu-%lu.bin",
						scale,degree,i,gen_sz,my_row,my_col);
			FILE *fid=fopen(filename,"rb");
			if(fid==NULL){fprintf(stderr, "Wrong fopen %s\n",filename);}
			off_t size_offset = fsize((const char *)filename);
			ranger_edge_count+=size_offset/sizeof(packed_edge);
			packed_edge *check_res = (packed_edge*) malloc(size_offset);
			fread(check_res,sizeof(packed_edge),size_offset/sizeof(packed_edge),fid);
			
			for(index_t k=0;k<size_offset/sizeof(packed_edge);k++)
				ranger_adj_card[check_res[k].v0-my_vert_off]++;
			
			fclose(fid);
			free(check_res);
	}

	/*alloc space*/
	#ifdef COMPRESS
	index_t ranger_edge_comp_count=((ranger_edge_count)>>2)*3+((ranger_edge_count)&0x03);
	#else
	index_t ranger_edge_comp_count=ranger_edge_count;
	#endif

	ranger_adj_list	= new data_t[ranger_edge_comp_count];

	/*compute beg_pos*/
	double tm=wtime();
	ranger_beg_pos[0]=0;
	for(index_t i=0;i<ranger_vert_count;i++)
		ranger_beg_pos[i+1]=ranger_beg_pos[i]+ranger_adj_card[i];
	
	if(ranger_beg_pos[ranger_vert_count]!=ranger_edge_count)
	{
		std::cout<<"Wrong edge count\n";
		std::cout<<ranger_vert_count<<"\t"<<ranger_beg_pos[ranger_vert_count]
			<<" vs "<<ranger_edge_count<<"\n\n\n\n";
		exit(-1);
	}

	if(rank==0) std::cout<<"beg_pos: "<<wtime()-tm<<" second(s)\n";
	
	/*read and process adj_list*/
	for(index_t i=0;i<gen_sz;i++)
	{
			sprintf(filename,"scale-%lu-degree-%lu-rank-%lu-of-%lu-par-rowcol-%lu-%lu.bin",
						scale,degree,i,gen_sz,my_row,my_col);
			FILE *fid=fopen(filename,"rb");
			if(fid==NULL){fprintf(stderr, "Wrong fopen\n");}
			off_t size_offset = fsize((const char *)filename);
			packed_edge *check_res = (packed_edge*) malloc(size_offset);
			fread(check_res,sizeof(packed_edge),size_offset/sizeof(packed_edge),fid);
			
			for(index_t k=0;k<size_offset/sizeof(packed_edge);k++)
			{
				adj_store(ranger_beg_pos[check_res[k].v0-my_vert_off]+
									ranger_offset[check_res[k].v0-my_vert_off],
									check_res[k].v1,ranger_adj_list);
				ranger_offset[check_res[k].v0-my_vert_off]++;
			}
			fclose(fid);
			free(check_res);

			//delete the file
			unlink(filename);
	}

	tm=wtime();
	sprintf(filename,"csr_%lu_%lu.%lu_%lu_of_%lu_%lu.bin",
					scale,degree,my_row,my_col,row_par, col_par);
	FILE *fid=fopen(filename, "wb");
	index_t ret=fwrite(ranger_adj_list,sizeof(data_t),ranger_edge_comp_count,fid);
	if(ret !=ranger_edge_comp_count)
	{
		std::cout<<"Dump adj-list wrong"<<rank<<"\n";
		exit(-1);
	}
	fclose(fid);
	if(rank==0) std::cout<<"Dump adj-list time: "<<wtime()-tm<<" second(s)\n";

	tm=wtime();
	sprintf(filename,"beg_%lu_%lu.%lu_%lu_of_%lu_%lu.bin",
					scale,degree,my_row,my_col,row_par, col_par);
	fid=fopen(filename, "wb");
	ret=fwrite(ranger_beg_pos,sizeof(index_t),ranger_vert_count+1,fid);
	if(ret !=ranger_vert_count+1)
	{
		std::cout<<"Dump beg-pos wrong"<<rank<<"\n";
		exit(-1);
	}
	fclose(fid);
	if(rank==0) std::cout<<"Dump beg-pos time: "<<wtime()-tm<<" second(s)\n";

	/*dump adj_list: text format*/
//  std::ofstream outfile("kron.sorted");
//  for(index_t i=0; i<vert_count; i++)
//	{
//    outfile<<i<<" "<<cbeg_pos[i+1]-cbeg_pos[i]<<" ";
//    for(index_t j=cbeg_pos[i]; j<cbeg_pos[i+1]; j++)
//      outfile<<adj_load(j,adj_list)<<" ";
//    outfile<<"\n";
//  }
//  outfile.close();
	 MPI_Finalize();	
	return true;
}

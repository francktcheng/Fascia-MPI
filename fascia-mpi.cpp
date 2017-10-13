// Copyright (c) 2013, The Pennsylvania State University.
// All rights reserved.
// 
// See COPYING for license.

using namespace std;

#include <mpi.h>
#include <omp.h>

int rank, nprocs;
bool verbose, debug;

#include <stdio.h>
#include <cstdlib>
#include <assert.h>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <unistd.h>
#include <climits>

#include "fascia.h"
#include "mpi_wrappers.hpp"
#include "xs1024star.hpp"
#include "graph.hpp"
#include "util.hpp"
#include "output.hpp"
#include "dynamic_table.hpp"
#include "dynamic_table_array.hpp"
#include "dynamic_table_part.hpp"
#include "partitioner.hpp"
#include "colorcount.hpp"
#include "colorcount_part.hpp"

bool timing = false;

void print_info_short(char* name)
{
  printf("\nTo run: %s [-g graphfile] [-t template || -b batchfile] [options]\n", name);
  printf("Help: %s -h\n\n", name);
}

void print_info(char* name)
{
  printf("\nTo run: %s [-g graphfile] [-t template || -b batchfile] [options]\n\n", name);

  printf("\tgraphfile = \n");
  printf("\t\tn\n");
  printf("\t\tm\n");
  printf("\t\tv0 v1\n");
  printf("\t\tv0 v2\n");
  printf("\t\t...\n");
  printf("\t\t(zero indexed)\n\n");

  printf("\tgraphfile (if labeled) = \n");
  printf("\t\tn\n");
  printf("\t\tm\n");
  printf("\t\tlabel_v0\n");
  printf("\t\tlabel_v1\n");
  printf("\t\t...\n");
  printf("\t\tv0 v1\n");
  printf("\t\tv0 v2\n");
  printf("\t\t...\n");
  printf("\t\t(zero indexed)\n\n"); 

  printf("\ttemplate =\n");
  printf("\t\tsame format as graphfile\n\n");

  printf("\tbatchfile =\n");
  printf("\t\ttemplateFile1\n");
  printf("\t\ttemplateFile2\n");
  printf("\t\t...\n");
  printf("\t\t(must supply only one of template file or batchfile)\n\n");

  printf("\toptions = \n");
  printf("\t\t-p  Perform partitioned count, default is distributed\n");
  printf("\t\t-o  Use outerloop parallelization\n");
  printf("\t\t-m  [#], compute counts for motifs of size #\n");
  printf("\t\t-l  Graph and template are labeled\n");
  printf("\t\t-i  [# iterations], default: 1\n");
  printf("\t\t-c  Output per-vertex counts to [template].vert\n");
  printf("\t\t-d  Output graphlet degree distribution to [template].gdd\n");
  printf("\t\t-a  Do not calculate automorphism of template\n");
  printf("\t\t\t(recommended when template size > 10)\n");
  printf("\t\t-r  Report runtime\n");
  printf("\t\t-v  Verbose output\n");
  printf("\t\t-h  Print this\n\n");
}


// void read_in_graph(Graph& g, char* graph_file, bool labeled, int*& labels_g)
// {  
//   ifstream file_g;
//   string line;
//   int* srcs_g;
//   int* dsts_g;
//   int n_g;
//   unsigned m_g;
//     
//   if (rank == 0)
//   {
//     file_g.open(graph_file);
//
//     // Read in labels and vertices for graph
//     getline(file_g, line);
//     n_g = atoi(line.c_str());
//     getline(file_g, line);
//     m_g = strtoul(line.c_str(), NULL, 10);
//
//     if (verbose) 
//       if (rank == 0) printf("n %d, m %u, p %d\n", n_g, m_g, nprocs);
//
//     srcs_g = new int[m_g];
//     dsts_g = new int[m_g];
//
//     if (labeled) {
//       labels_g = new int[n_g];
//       for (int i = 0; i < n_g; ++i)
//       {
//         getline(file_g, line);  
//         labels_g[i] = atoi(line.c_str());
//       }
//     }
//
//     for (unsigned i = 0; i < m_g; ++i)
//     {
//       getline(file_g, line, ' ');    
//       srcs_g[i] = atoi(line.c_str());
//       getline(file_g, line);  
//       dsts_g[i] = atoi(line.c_str());
//     }
//
//     file_g.close();
//
//     if (nprocs > 1)
//     {
//       MPI_Bcast(&n_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
//       MPI_Bcast(&m_g, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//       MPI_Bcast_chunk(srcs_g, m_g, 0, rank);
//       MPI_Bcast_chunk(dsts_g, m_g, 0, rank);
//
//       if (labeled)
//         MPI_Bcast_chunk(labels_g, n_g, 0, rank);
//     }
//   }
//   else
//   {
//     MPI_Bcast(&n_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     MPI_Bcast(&m_g, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//     srcs_g = new int[m_g];
//     dsts_g = new int[m_g];
//     MPI_Bcast_chunk(srcs_g, m_g, 0, rank);
//     MPI_Bcast_chunk(dsts_g, m_g, 0, rank);
//
//     if (labeled) {
//       labels_g = new int[n_g];
//       MPI_Bcast_chunk(labels_g, n_g, 0, rank);
//     }
//   }
//
//   g.init(n_g, m_g, srcs_g, dsts_g);
//   delete [] srcs_g;
//   delete [] dsts_g;
// }

/**
 * @brief Another read in graph function to deal with non-consecutive 
 * file format 
 *
 * @param g
 * @param graph_file
 * @param labeled
 * @param labels_g
 */
void read_in_graph(Graph& g, char* graph_file, bool labeled, int*& labels_g)
{  
  ifstream file_g;
  string line;
  int* srcs_g;
  int* dsts_g;
  int n_g;
  unsigned m_g;
    
  if (rank == 0)
  {
    file_g.open(graph_file);

    // Read in labels and vertices for graph
    getline(file_g, line);
    n_g = atoi(line.c_str());
    getline(file_g, line);
    m_g = strtoul(line.c_str(), NULL, 10);

    if (verbose) 
      if (rank == 0) printf("n %d, m %u, p %d\n", n_g, m_g, nprocs);

    srcs_g = new int[m_g];
    dsts_g = new int[m_g];

    if (labeled) {
      labels_g = new int[n_g];
      for (int i = 0; i < n_g; ++i)
      {
        getline(file_g, line);  
        labels_g[i] = atoi(line.c_str());
      }
    }

    int max_id = 0;

    for (unsigned i = 0; i < m_g; ++i)
    {
      getline(file_g, line, ' ');    
      srcs_g[i] = atoi(line.c_str());

      max_id = (srcs_g[i] > max_id ) ? srcs_g[i] : max_id;

      getline(file_g, line);  
      dsts_g[i] = atoi(line.c_str());

      max_id = (dsts_g[i] > max_id ) ? dsts_g[i] : max_id;

    }

    file_g.close();

    if (max_id != n_g - 1)
    {
        //remove the "holes" in vertices ids
        if (verbose) 
           printf("Start remove holes; max_id: %d, n_g: %d\n", max_id, n_g);

        int* v_id = new int[max_id+1];
        std::memset(v_id, 0, (max_id+1)*sizeof(int));

        for(int i=0;i<m_g;i++)
        {
            v_id[srcs_g[i]] = 1;
            v_id[dsts_g[i]] = 1;
        }

        int itr = 0;
        for(int i=0;i<max_id+1;i++)
        {
            if (v_id[i] == 1)
                v_id[i] = (itr++);
        }

        for(int i=0;i<m_g;i++)
        {
            srcs_g[i] = v_id[srcs_g[i]];
            dsts_g[i] = v_id[dsts_g[i]];
        }

        if (verbose) 
          printf("Finish remove holes\n");

        delete[] v_id;
    }

    if (nprocs > 1)
    {
      MPI_Bcast(&n_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&m_g, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      MPI_Bcast_chunk(srcs_g, m_g, 0, rank);
      MPI_Bcast_chunk(dsts_g, m_g, 0, rank);

      if (labeled)
        MPI_Bcast_chunk(labels_g, n_g, 0, rank);
    }
  }
  else
  {
    MPI_Bcast(&n_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m_g, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    srcs_g = new int[m_g];
    dsts_g = new int[m_g];
    MPI_Bcast_chunk(srcs_g, m_g, 0, rank);
    MPI_Bcast_chunk(dsts_g, m_g, 0, rank);

    if (labeled) {
      labels_g = new int[n_g];
      MPI_Bcast_chunk(labels_g, n_g, 0, rank);
    }
  }

  g.init(n_g, m_g, srcs_g, dsts_g);
  delete [] srcs_g;
  delete [] dsts_g;
}

void read_in_graph(Graph& g, char* graph_file, bool labeled, int*& labels_g,
  int*& part_offsets)
{  

   read_in_graph(g, graph_file, labeled, labels_g);

  if (rank == 0)
  {
    part_offsets = new int[nprocs+1];    
    for (int i = 0; i < nprocs; ++i) {
      part_offsets[i] = i * (g.num_vertices() / nprocs);
    }
    part_offsets[nprocs] = g.num_vertices();

    MPI_Bcast(&nprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(part_offsets, nprocs+1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else
  {
    MPI_Bcast(&nprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    part_offsets = new int[nprocs+1];
    MPI_Bcast(part_offsets, nprocs+1, MPI_INT, 0, MPI_COMM_WORLD);
  }
}

void run_dist(char* graph_file, char* template_file, bool labeled,
  bool do_vert, bool do_gdd, int iterations, bool calc_auto, int omp_thds)
{
  Graph g;
  Graph t;
  int* labels_g = NULL;
  int* labels_t = NULL;
  int* part_offsets = NULL;
  char* vert_file = new char[1024];
  char* gdd_file = new char[1024];

  if (do_vert) {
    strcat(vert_file, template_file);
    strcat(vert_file, ".vert");
  }
  if (do_gdd) {
    strcat(gdd_file, template_file);
    strcat(gdd_file, ".gdd");
  }

  double elt = 0.0;
  if (timing || verbose) {
    elt = timer();
  }

  read_in_graph(g, graph_file, labeled, labels_g, part_offsets);
  read_in_graph(t, template_file, labeled, labels_t);

  colorcount_part graph_count;
  graph_count.init(g, part_offsets, nprocs, labels_g, labeled,
                    calc_auto, do_gdd, do_vert, omp_thds);
  double full_count = graph_count.do_full_count(&t, labels_t, iterations);

  if (do_gdd || do_vert)
  {
    double* vert_counts = new double[g.num_vertices()];
    double* local_vert_counts = graph_count.get_vert_counts();
    for (int i = 0; i < g.num_vertices(); ++i)
      vert_counts[i] = 0.0;
    for (int i = part_offsets[rank]; i < part_offsets[rank+1]; ++i)
      vert_counts[i] = local_vert_counts[i - part_offsets[rank]];

    output out(vert_counts, g.num_vertices());

    double* full_vert_counts = new double[g.num_vertices()];
    for (int i = 0; i < g.num_vertices(); ++i)
      full_vert_counts[i] = 0.0;

    MPI_Reduce(out.vert_counts, full_vert_counts, g.num_vertices(),
      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    out.vert_counts = full_vert_counts;

    if (rank == 0) {
      if (do_gdd) {
        out.output_gdd(gdd_file);
        free(gdd_file);
      } 
      if (do_vert) {        
        out.output_verts(vert_file);
        free(vert_file);
      }
    }

    delete [] vert_counts;
    delete [] full_vert_counts;
  }

  if (rank == 0) printf("Count:\n\t%e\n", full_count);

if (timing || verbose) {
  elt = timer() - elt;
  if (rank == 0) printf("Total time:\n\t%9.6lf seconds\n", elt);
}
  
  if (labels_g != NULL) delete [] labels_g;
  if (labels_g != NULL) delete [] labels_t;
}

void run_single(char* graph_file, char* template_file, bool labeled,
                bool do_vert, bool do_gdd,
                int iterations, 
                bool do_outerloop, bool calc_auto, int omp_thds)
{
  Graph g;
  Graph t;
  int* labels_g = NULL;
  int* labels_t = NULL;
  char* vert_file = new char[1024];
  char* gdd_file = new char[1024];

  if (do_vert) {
    strcat(vert_file, template_file);
    strcat(vert_file, ".vert");
  }
  if (do_gdd) {
    strcat(gdd_file, template_file);
    strcat(gdd_file, ".gdd");
  }

  double elt = 0.0;
  if (timing || verbose) {
    elt = timer();
  }

  read_in_graph(g, graph_file, labeled, labels_g);
  read_in_graph(t, template_file, labeled, labels_t);

  int task_iterations = int( (double)iterations / (double)nprocs + 0.5);
  if (task_iterations < 1 || iterations == 1) task_iterations = 1;
  if (verbose && rank == 1)
    printf("Running %d iteration(s), %d iteration(s) per task\n", iterations, task_iterations);

  double task_count = 0.0;
  double full_count = 0.0;

  if (do_outerloop)
  {
    // int num_threads = omp_get_max_threads();
    int num_threads = omp_thds;
    int iter = int( (double)task_iterations / (double)num_threads + 0.5);
    if (iter < 1)
      iter = 1;
    
    colorcount* graph_count = new colorcount[num_threads];
    for (int tid = 0; tid < num_threads; ++tid) {
      graph_count[tid].init(g, labels_g, labeled, 
                            calc_auto, do_gdd, do_vert, num_threads);
    }

    double** vert_counts;
    if (do_gdd || do_vert)
      vert_counts = new double*[num_threads];

#pragma omp parallel num_threads(num_threads) reduction(+:task_count)
{
    int tid = omp_get_thread_num();
    task_count += graph_count[tid].do_full_count(&t, labels_t, iter);
    if (do_gdd || do_vert)
      vert_counts[tid] = graph_count[tid].get_vert_counts();
}   

    task_count /= (double)num_threads;
    MPI_Reduce(&task_count, &full_count, 1, 
      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    full_count /= (double)nprocs;

    if (do_gdd || do_vert)
    {
      output out(vert_counts, num_threads, g.num_vertices());

      double* full_vert_counts = new double[g.num_vertices()];
      for (int i = 0; i < g.num_vertices(); ++i)
        full_vert_counts[i] = 0.0;

      MPI_Reduce(out.vert_counts, full_vert_counts, g.num_vertices(),
        MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      for (int i = 0; i < g.num_vertices(); ++i)
        full_vert_counts[i] /= (double)nprocs;
      out.vert_counts = full_vert_counts;

      if (rank == 0) {
        if (do_gdd) {
          out.output_gdd(gdd_file);
          free(gdd_file);
        } 
        if (do_vert) {        
          out.output_verts(vert_file);
          free(vert_file);
        }
      }

      delete [] vert_counts;
      delete [] full_vert_counts;
    }
  }
  else
  {
    colorcount graph_count;
    graph_count.init(g, labels_g, labeled, 
                      calc_auto, do_gdd, do_vert, omp_thds);

    task_count += graph_count.do_full_count(&t, labels_t, task_iterations);

    MPI_Reduce(&task_count, &full_count, 1, 
      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    full_count /= (double)nprocs;

    if (do_gdd || do_vert)
    {
      double* vert_counts = graph_count.get_vert_counts();
      output out(vert_counts, g.num_vertices());

      double* full_vert_counts = new double[g.num_vertices()];
      for (int i = 0; i < g.num_vertices(); ++i)
        full_vert_counts[i] = 0.0;

      MPI_Reduce(out.vert_counts, full_vert_counts, g.num_vertices(),
        MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      for (int i = 0; i < g.num_vertices(); ++i)
        full_vert_counts[i] /= (double)nprocs;
      out.vert_counts = full_vert_counts;

      if (rank == 0) {
        if (do_gdd) {
          out.output_gdd(gdd_file);
          free(gdd_file);
        } 
        if (do_vert) {        
          out.output_verts(vert_file);
          free(vert_file);
        }
      }

      delete [] full_vert_counts;
    }
  }

  if (rank == 0) printf("Count:\n\t%e\n", full_count);

if (timing || verbose) {
  elt = timer() - elt;
  if (rank == 0) printf("Total time:\n\t%9.6lf seconds\n", elt);
}

  if (labels_g != NULL) delete [] labels_g;
  if (labels_g != NULL) delete [] labels_t;
}


void run_batch(char* graph_file, char* batch_file, bool labeled,
                    bool do_vert, bool do_gdd,
                    int iterations, 
                    bool do_outerloop, bool calc_auto, int omp_thds)
{
  Graph g;
  Graph t;
  int* labels_g = NULL;
  char* vert_file = NULL;
  char* gdd_file = NULL;

  read_in_graph(g, graph_file, labeled, labels_g);

  double elt = 0.0;
  if (timing || verbose) {
    elt = timer();
  }

  ifstream if_batch;
  string line;
  if_batch.open(batch_file);
  while (getline(if_batch, line))
  {   
    int* labels_t = NULL;
    char* template_file = strdup(line.c_str());
    read_in_graph(t, template_file, labeled, labels_t);

    int task_iterations = int( (double)iterations / (double)nprocs + 0.5);
    if (task_iterations < 1 || iterations == 1) task_iterations = 1;
    if (verbose && rank == 0)
      printf("Running %d iteration(s), %d iteration(s) per task\n", iterations, task_iterations);

    double task_count = 0.0;
    double full_count = 0.0;

    if (do_outerloop)
    {
      int num_threads = omp_get_max_threads();
      int iter = int( (double)task_iterations / (double)num_threads + 0.5);
      if (iter < 1) iter = 1;
      
      colorcount* graph_count = new colorcount[num_threads];
      for (int tid = 0; tid < num_threads; ++tid) {
        graph_count[tid].init(g, labels_g, labeled, 
                              calc_auto, do_gdd, do_vert, num_threads);
      }

      double** vert_counts;
      if (do_gdd || do_vert)
        vert_counts = new double*[num_threads];

#pragma omp parallel num_threads(num_threads) reduction(+:task_count)
{
      int tid = omp_get_thread_num();
      task_count += graph_count[tid].do_full_count(&t, labels_t, iter);
      if (do_gdd || do_vert)
        vert_counts[tid] = graph_count[tid].get_vert_counts();
}      

      task_count /= (double)num_threads;
      MPI_Reduce(&task_count, &full_count, 1, 
        MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      full_count /= (double)nprocs;

      if (do_gdd || do_vert)
      {
        output out(vert_counts, num_threads, g.num_vertices());

        double* full_vert_counts = new double[g.num_vertices()];
        for (int i = 0; i < g.num_vertices(); ++i)
          full_vert_counts[i] = 0.0;

        MPI_Reduce(out.vert_counts, full_vert_counts, g.num_vertices(),
          MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        for (int i = 0; i < g.num_vertices(); ++i)
          full_vert_counts[i] /= (double)nprocs;
        out.vert_counts = full_vert_counts;

        if (rank == 0) {
          if (do_gdd) {
            gdd_file = strdup(template_file);
            strcat(gdd_file, ".gdd");
            out.output_gdd(gdd_file);
            free(gdd_file);
          } 
          if (do_vert) {        
            vert_file = strdup(template_file);
            strcat(vert_file, ".vert");
            out.output_verts(vert_file);
            free(vert_file);
          }
        }

        delete [] vert_counts;
        delete [] full_vert_counts;
      }
    }
    else
    {
      colorcount graph_count;
      graph_count.init(g, labels_g, labeled, 
                        calc_auto, do_gdd, do_vert, omp_thds);
      task_count += graph_count.do_full_count(&t, labels_t, task_iterations);

      MPI_Reduce(&task_count, &full_count, 1, 
        MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      full_count /= (double)nprocs;

      if (do_gdd || do_vert)
      {
        double* vert_counts = graph_count.get_vert_counts();
        output out(vert_counts, g.num_vertices());

        double* full_vert_counts = new double[g.num_vertices()];
        for (int i = 0; i < g.num_vertices(); ++i)
          full_vert_counts[i] = 0.0;

        MPI_Reduce(out.vert_counts, full_vert_counts, g.num_vertices(),
          MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        for (int i = 0; i < g.num_vertices(); ++i)
          full_vert_counts[i] /= (double)nprocs;
        out.vert_counts = full_vert_counts;

        if (rank == 0) {
          if (do_gdd) {
            gdd_file = strdup(template_file);
            strcat(gdd_file, ".gdd");
            out.output_gdd(gdd_file);
            free(gdd_file);
          } 
          if (do_vert) {        
            vert_file = strdup(template_file);
            strcat(vert_file, ".vert");
            out.output_verts(vert_file);
            free(vert_file);
          }
        }

        delete [] full_vert_counts;
      }
    }

    if (rank == 0) printf("%e\n", full_count);  

    if (labels_t != NULL) delete [] labels_t;
    delete [] template_file;
  }

if (timing || verbose) {
  elt = timer() - elt;
  printf("Total time:\n\t%9.6lf seconds\n", elt);
}

  delete [] labels_g;
}


void run_motif(char* graph_file, int motif, 
                bool do_vert, bool do_gdd, 
                int iterations, 
                bool do_outerloop, bool calc_auto, int omp_thds)
{
  char* motif_batchfile = NULL;

  switch(motif)
  {
    case(3): motif_batchfile = strdup("motif/graphs_n3_1/batchfile"); break;
    case(4): motif_batchfile = strdup("motif/graphs_n4_2/batchfile"); break;
    case(5): motif_batchfile = strdup("motif/graphs_n5_3/batchfile"); break;
    case(6): motif_batchfile = strdup("motif/graphs_n6_6/batchfile"); break;
    case(7): motif_batchfile = strdup("motif/graphs_n7_11/batchfile"); break;
    case(8): motif_batchfile = strdup("motif/graphs_n8_23/batchfile"); break;
    case(9): motif_batchfile = strdup("motif/graphs_n9_47/batchfile"); break;
    case(10): motif_batchfile = strdup("motif/graphs_n10_106/batchfile"); break;
    default: break;
  }

  run_batch(graph_file, motif_batchfile, false,
            do_vert, do_gdd,
            iterations, 
            do_outerloop, calc_auto, omp_thds);
}


int main(int argc, char** argv)
{ 
  srand(time(0));
  setbuf(stdout, NULL);

  MPI_Init(&argc, &argv);
  // int provided;
  // MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  // if (provided < MPI_THREAD_MULTIPLE) {
  //     // Error - MPI does not provide needed threading level
  // }

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char* graph_file = NULL;
  char* template_file = NULL;
  char* batch_file = NULL;
  int iterations = 1;
  bool do_outerloop = false;
  bool calculate_automorphism = true;
  bool labeled = false;
  bool do_gdd = false;
  bool do_vert = false;
  verbose = false;
  bool distributed_count = true;
  bool partitioned_count = false;
  int omp_thds = omp_get_max_threads();
  int motif = 0;

  char c;
  while ((c = getopt (argc, argv, "g:t:b:i:s:m:acdvrohlp")) != -1)
  {
    switch (c)
    {
      case 'h': print_info(argv[0]); break;
      case 'l': labeled = true; break;
      case 'g': graph_file = strdup(optarg); break;
      case 't': template_file = strdup(optarg); break;
      case 'b': batch_file = strdup(optarg); break;
      case 'i': iterations = atoi(optarg); break;
      case 'm': motif = atoi(optarg); break;
      case 's': omp_thds = atoi(optarg); break;
      case 'a': calculate_automorphism = false; break;
      case 'c': do_vert = true; break;
      case 'd': do_gdd = true; break;
      case 'o': do_outerloop = true; break;
      case 'v': verbose = true; break;
      case 'r': timing = true; break;
      case 'p': partitioned_count = true; distributed_count = false; break;
      case '?':
        if (optopt == 'g' || optopt == 't' || optopt == 'b' || optopt == 'i' || optopt == 'm')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n",
      optopt);
        print_info(argv[0]);
      default:
        abort();
    }
  } 

  if(argc < 3) {
    if (rank == 0) print_info_short(argv[0]);
    MPI_Finalize(); exit(0);
  }
  if (motif && (motif < 3 || motif > 10)) {
    if (rank == 0) {
      printf("\nMotif option must be between [3,10]\n");
      print_info(argv[0]);
    }
    MPI_Finalize(); exit(0);
  }
  if (graph_file == NULL) { 
    if (rank == 0) {
      printf("\nMust supply graph file\n");    
      print_info(argv[0]);
    }
    MPI_Finalize(); exit(0);
  }
  if (template_file == NULL && batch_file == NULL && !motif) {
    if (rank == 0) {
      printf("\nMust supply template XOR batchfile or -m option\n");
      print_info(argv[0]);
    }
    MPI_Finalize(); exit(0);
  }
  if (template_file != NULL && batch_file != NULL) {
    if (rank == 0) {
      printf("\nMust only supply template file XOR batch file\n");
      print_info(argv[0]);
    }
    MPI_Finalize(); exit(0);
  }
  if (iterations < 1) {
    if (rank == 0) {
      printf("\nNumber of iterations must be positive\n");    
      print_info(argv[0]);
    }
    MPI_Finalize(); exit(0);
  }

  if (partitioned_count && nprocs == 1) {
    if (rank == 0) {
      printf("\nWarning: requesting partitioned counting with only 1 proc\n");
      printf("==> Defaulting to thread parallel method\n\n");
    }
    partitioned_count = false;
  }
  if (partitioned_count && do_outerloop) {
    if (rank == 0) {
      printf("\nWarning: requesting outer loop parallelism for partitioned counting\n");
      printf("==> Defaulting to inner loop parallelism\n\n");
    }
    do_outerloop = false;
  } 
  if ((template_file != NULL || batch_file != NULL) && motif) {
    if (rank == 0) {
      printf("\nWarning: requesting motif count with batch/template file count\n");
      printf("==> Defaulting to use batch/template file\n\n");
    }
    motif = 0;
  }

  if (partitioned_count) {
    run_dist(graph_file, template_file, labeled, 
      do_vert, do_gdd, 
      iterations, calculate_automorphism, omp_thds);
  }
  else {
    if (motif) {
      run_motif(graph_file, motif, 
                do_vert, do_gdd, 
                iterations, do_outerloop, calculate_automorphism, omp_thds);
    }
    else if (template_file != NULL) {
      run_single(graph_file, template_file, labeled,                
                  do_vert, do_gdd,
                  iterations, do_outerloop, calculate_automorphism, omp_thds);
    }
    else if (batch_file != NULL) {
      run_batch(graph_file, batch_file, labeled,
                  do_vert, do_gdd,
                  iterations, do_outerloop, calculate_automorphism, omp_thds);
    }
  }

  MPI_Finalize();

  return 0;
}

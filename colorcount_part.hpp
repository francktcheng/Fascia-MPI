
using namespace std;

class colorcount_part{
public:
    
  colorcount_part() {};  
  ~colorcount_part()
  {  }
  
/*
'####:'##::: ##:'####:'########:
. ##:: ###:: ##:. ##::... ##..::
: ##:: ####: ##:: ##::::: ##::::
: ##:: ## ## ##:: ##::::: ##::::
: ##:: ##. ####:: ##::::: ##::::
: ##:: ##:. ###:: ##::::: ##::::
'####: ##::. ##:'####:::: ##::::
....::..::::..::....:::::..:::::
*/
  void init(Graph& full_graph, int* Part_offsets, int num_parts,
    int* labels, bool label,
    bool calc_auto, bool do_gdd, bool do_vert)
  {
    g = &full_graph;
    labels_g = labels;
    part_offsets = Part_offsets;
    labeled = label;    
    do_graphlet_freq = do_gdd;
    do_vert_output = do_vert;

    begin_vert = part_offsets[rank];
    end_vert = part_offsets[rank+1];
    num_verts_part = end_vert - begin_vert;

    sizes_verts = new int[end_vert-begin_vert];
    comm_sizes_send = new unsigned long[num_parts];
    comm_sizes_rec = new unsigned long[num_parts];

    comm_num_rec = new int[num_parts];
    comm_num_send = new int[num_parts];
    comm_verts_rec = new int*[num_parts];
    comm_verts_send = new int*[num_parts];
    vector<int>* temp_send = new vector<int>[num_parts];
    vector<int>* temp_rec = new vector<int>[num_parts];

    for (int v = begin_vert; v < end_vert; ++v)
    {
      int* outs = g->adjacent_vertices(v);
      int out_degree = g->out_degree(v);
      for (int j = 0; j < out_degree; ++j)
      {
        int out = outs[j];
        int out_part = get_part(out, part_offsets, num_parts);
        temp_send[out_part].push_back(v);
        temp_rec[out_part].push_back(out);
      }
    }

    for (int i = 0; i < num_parts; ++i)
    {
      int* temp_verts = temp_send[i].data();
      int temp_size = temp_send[i].size();
      if (temp_size)
      {
        quicksort_inc(temp_verts, 0, temp_size-1);
        int num_multiple = 0;
        for (int j = 1; j < temp_size; ++j)
          if (temp_verts[j] == temp_verts[j-1])
            ++num_multiple;

        int new_size = temp_size - num_multiple;
        comm_num_send[i] = new_size;
        comm_verts_send[i] = new int[new_size];
        comm_verts_send[i][0] = temp_verts[0];
        int counter = 1;
        for (int j = 1; j < temp_size; ++j)
          if (temp_verts[j] != temp_verts[j-1])
            comm_verts_send[i][counter++] = temp_verts[j];
        assert(counter == new_size);
      }
      else
      {
        comm_num_send[i] = 0;
        comm_verts_send[i] = NULL;
      }
    }

    for (int i = 0; i < num_parts; ++i)
    {
      int* temp_verts = temp_rec[i].data();
      int temp_size = temp_rec[i].size();
      if (temp_size)
      {
        quicksort_inc(temp_verts, 0, temp_size-1);
        int num_multiple = 0;
        for (int j = 1; j < temp_size; ++j)
          if (temp_verts[j] == temp_verts[j-1])
            ++num_multiple;

        int new_size = temp_size - num_multiple;
        comm_num_rec[i] = new_size;
        comm_verts_rec[i] = new int[new_size];
        comm_verts_rec[i][0] = temp_verts[0];
        int counter = 1;
        for (int j = 1; j < temp_size; ++j)
          if (temp_verts[j] != temp_verts[j-1])
            comm_verts_rec[i][counter++] = temp_verts[j];
        assert(counter == new_size);
      }
      else
      {
        comm_num_rec[i] = 0;
        comm_verts_rec[i] = NULL;
      }
    }

    delete [] temp_send;
    delete [] temp_rec;

    if (do_graphlet_freq || do_vert_output)
    {
      final_vert_counts = new double[num_verts_part];
      for (int i = 0; i < num_verts_part; ++i)
        final_vert_counts[i] = 0.0;
    }

    if (verbose)
      printf("rank: %d  begin: %d  end: %d\n", rank, begin_vert, end_vert);
  }
  
/*
'########:'##::::'##:'##:::::::'##:::::::
 ##.....:: ##:::: ##: ##::::::: ##:::::::
 ##::::::: ##:::: ##: ##::::::: ##:::::::
 ######::: ##:::: ##: ##::::::: ##:::::::
 ##...:::: ##:::: ##: ##::::::: ##:::::::
 ##::::::: ##:::: ##: ##::::::: ##:::::::
 ##:::::::. #######:: ########: ########:
..:::::::::.......:::........::........::
*/    
  double do_full_count(Graph* sub_graph, int* labels, int N)
  {  
    if (verbose && rank == 0) printf("%d Beginning partitioning\n", rank);

    num_iter = N;
    t = sub_graph;
    labels_t = labels;  
    
    part = partitioner(*t, labeled, labels_t);
    part.sort_subtemplates(); 
    
    num_colors = t->num_vertices();
    num_verts_graph = g->num_vertices();
    max_degree = g->max_degree();
    subtemplates = part.get_subtemplates();
    subtemplate_count = part.get_subtemplate_count();    

    if (verbose && rank == 0) 
      printf("%d Done partitioning %d %d\n", rank, num_colors, subtemplate_count);   

    create_tables();
    if (verbose && rank == 0) 
      printf("%d Done creating tables\n", rank);
  
    dt.init_part(subtemplates, subtemplate_count, g->num_vertices(), num_colors, begin_vert, end_vert, rank, nprocs);

    if (verbose && rank == 0) 
      printf("%d n %d, max degree %d\n", rank, num_verts_graph, max_degree); 
 
    double count = 0.0;
    transfer_time_sum = 0.0;
    transfer_size_sum = 0.0;

    for (int j = 0; j < num_iter; j++)
    {      
      double elt = timer();

      count += template_count();

      if (verbose && rank == 0) {     
        elt = timer() - elt;
        if (rank == 0)      
        {
          printf("Total time: %9.6lf\n", elt);
          printf("Computation time: %9.6lf\n", elt - transfer_time_sum);
          printf("Transfer time: %9.6lf\n", transfer_time_sum);
          printf("Transfer data: %9.6lf\n", transfer_size_sum/CHUNK_SIZE);
        }
      }
    }
    
    double final_count = count /= (double)num_iter;
    double prob_colorful = factorial(num_colors) / 
        ( factorial(num_colors - t->num_vertices()) * pow(num_colors, t->num_vertices()) );  
    int num_auto = t->num_vertices() <= 10 ? count_automorphisms(*t) : 1;
    final_count = floor(final_count / (prob_colorful * (double)num_auto) + 0.5);

    if (verbose && rank == 0) {
      printf("Probability colorful: %9.6lf\n", prob_colorful);
      printf("Num automorphisms: %d\n", num_auto);
    }

    if (do_graphlet_freq || do_vert_output)
    {
      for (int i = 0; i < num_verts_part; ++i)
        final_vert_counts[i] = 
            (double)floor( final_vert_counts[i] / ((double)num_auto * 
            (double)num_iter * prob_colorful) + 0.5);
    }

    delete_tables();
    part.clear_temparrays(); 

    return final_count;      
  }

  double* get_vert_counts()
  {
    return final_vert_counts;
  }

private:

/*
'########:'########:'##::::'##:'########:: ##:::::::
... ##..:: ##.....:: ###::'###: ##.... ##: ##:::::::
::: ##:::: ##::::::: ####'####: ##:::: ##: ##:::::::
::: ##:::: ######::: ## ### ##: ########:: ##:::::::
::: ##:::: ##...:::: ##. #: ##: ##.....::: ##:::::::
::: ##:::: ##::::::: ##:.:: ##: ##:::::::: ##:::::::
::: ##:::: ########: ##:::: ##: ##:::::::: ########:
:::..:::::........::..:::::..::..:::::::::........::
*/

  double template_count()
  {  
    double full_count;
    int num_verts = g->num_vertices();
    colors_g = new int[num_verts];

    if (rank == 0) {
#pragma omp parallel 
{
      /* thread-local RNG initialization */    
      xs1024star_t xs;
      xs1024star_seed((unsigned long)rand() + omp_get_thread_num(), &xs);

  #pragma omp for
      for (int v = 0; v < num_verts; ++v)
      {
        int color = (int)(xs1024star_next(&xs) % (unsigned long)num_colors);
        colors_g[v] = color;
      }
}
    }

    MPI_Bcast(colors_g, num_verts, MPI_INT, 0, MPI_COMM_WORLD);

    if (verbose && rank == 0) printf("Initing table node\n");

    int table_node = subtemplate_count - 1;
    dt.set_table_node(table_node);
    init_table_node(table_node);

    for (int s = subtemplate_count - 2; s >= 0 ; --s)
    {
      set_count = 0;
      total_count = 0;
      read_count = 0;

      int num_verts_sub = num_verts_table[s];
      int a = part.get_active_index(s);
      int p = part.get_passive_index(s);

      if (verbose && rank == 0) {
        printf("\nIniting with sub %d, verts: %d\n", s, num_verts_sub);
        printf("Active %d, Passive: %d\n", a, p);
      }    

      double elt = 0.0;
      float cc = 0.0;

      if (num_verts_sub == 1)
      {
        if (verbose && rank == 0) elt = timer();
      
        dt.set_sub_to_table_node(s);

        if (verbose && rank == 0) { 
          elt = timer() - elt;
          if (rank == 0)
            printf("%d s %d, it node %9.6lf s.\n", rank, s, elt);
        }
      }
      else
      {  
        if (verbose && rank == 0) elt = timer();

        dt.init_sub(s, a, p);
        cc = colorful_count_array_part(s);

        if (verbose && rank == 0) { 
          elt = timer() - elt;
          if (rank == 0)
            printf("%d s %d, array time %9.6lf s.\n", rank, s, elt);
        }
      }

      if (num_verts_sub > 1)
      {
        if (verbose && cc > 0.0)
          printf("Rank %d count: %9.6f\n", rank, cc);  

        MPI_Allreduce(MPI_IN_PLACE, &cc, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  
        if (verbose && cc > 0.0 && rank == 0)
          printf("Total count: %9.6f\n", cc);

        if (num_verts_table[a] > 1) {
          dt.clear_sub(a);
        }
        if (num_verts_table[p] > 1) {
          dt.clear_sub(p);
        }

        if (s == 0)
        {
          full_count = cc;
          break;
        }
        else if (nprocs > 1)
        {
          double trans_time = 0.0;          
          if (rank == 0 && verbose)
            trans_time = timer();
    
          for (int i = 0; i < nprocs; ++i) {     
            comm_sizes_send[i] = 0;
            for (int j = 0; j < comm_num_send[i]; ++j)
            {    
              int vert = comm_verts_send[i][j];
              comm_sizes_send[i] += (unsigned long)sizes_verts[vert-begin_vert];
            }            
          }

          MPI_Alltoall(comm_sizes_send, 1, MPI_UNSIGNED_LONG,
             comm_sizes_rec, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

          dt.init_comp_table(comm_sizes_rec);

          // Communicate counts tables among adjacent parts
          // TODO: replace this with alltoallv
          // --or option to use it over send/recv
          for (int i = 0; i < nprocs; ++i) 
          {
            for (int j = 0; j < nprocs; ++j) 
            {
              if (rank == i && rank == j) // compress own table
              {
                unsigned long count = comm_sizes_send[j];
                int num_vert_send = comm_num_send[j];
                int num_vert_rec = comm_num_rec[i];

                float* counts_comp = new float[count];
                short* colorsets_comp = new short[count];
                unsigned long* offsets_comp = new unsigned long[num_vert_send+1];
                dt.compress_for_send(count, 
                  num_vert_send, comm_verts_send[j], j,
                  counts_comp, colorsets_comp, offsets_comp);                 
                dt.append_to_table(count, num_vert_rec, 
                  part_offsets, i, comm_verts_rec[i],
                  counts_comp, colorsets_comp, offsets_comp);

                delete [] counts_comp;
                delete [] colorsets_comp;
                delete [] offsets_comp;
              }
              else if (rank == i)  //sending
              {
                unsigned long count = comm_sizes_send[j];
                if (count)
                {
                  unsigned long num_vert_send = (unsigned long)comm_num_send[j];
                  float* counts_comp = new float[count];
                  short* colorsets_comp = new short[count];
                  unsigned long* offsets_comp = new unsigned long[num_vert_send+1];
                  dt.compress_for_send(count, 
                    num_vert_send, comm_verts_send[j], j,
                    counts_comp, colorsets_comp, offsets_comp);

                  MPI_Send_chunk(counts_comp, count, j, rank);
                  MPI_Send_chunk(colorsets_comp, count, j, rank);
                  MPI_Send_chunk(offsets_comp, num_vert_send+1, j, rank);

                  delete [] counts_comp;
                  delete [] colorsets_comp;
                  delete [] offsets_comp;
                  transfer_size_sum += count*(4+2) + num_vert_send*8;
                }
              }
              else if (rank == j) // receiving
              { 
                unsigned long count = comm_sizes_rec[i];
                if (count)
                {
                  unsigned long num_vert_rec = (unsigned long)comm_num_rec[i];
                  float* counts_comp = new float[count];
                  short* colorsets_comp = new short[count];
                  unsigned long* offsets_comp = new unsigned long[num_vert_rec+1];

                  MPI_Recv_chunk(counts_comp, count, i, rank);
                  MPI_Recv_chunk(colorsets_comp, count, i, rank);
                  MPI_Recv_chunk(offsets_comp, num_vert_rec+1, i, rank);

                  dt.append_to_table(count, num_vert_rec, 
                    part_offsets, i, comm_verts_rec[i],
                    counts_comp, colorsets_comp, offsets_comp);
                  delete [] counts_comp;
                  delete [] colorsets_comp;
                  delete [] offsets_comp;
                }
              }
            }
          }

          MPI_Barrier(MPI_COMM_WORLD);
          dt.finalize();

          if (verbose && rank == 0)
          {
            trans_time = timer() - trans_time;
            printf("Transfer time: %9.6lf\n", trans_time);
            transfer_time_sum += trans_time;
          }
        }
        else
        {  
          unsigned long count = comm_sizes_send[0];
          int num_vert_send = comm_num_send[0];
          int num_vert_rec = comm_num_rec[0];
          float* counts_comp = new float[count];
          short* colorsets_comp = new short[count];
          unsigned long* offsets_comp = new unsigned long[num_vert_send+1];
          dt.compress_for_send(count, 
            num_vert_send, comm_verts_send[0], 0,
            counts_comp, colorsets_comp, offsets_comp);                 
          dt.append_to_table(count, num_vert_rec, 
            part_offsets, 0, comm_verts_rec[0],
            counts_comp, colorsets_comp, offsets_comp);
          delete [] counts_comp;
          delete [] colorsets_comp;
          delete [] offsets_comp;
        }
      }
    } // end for s in subtemplates

    return full_count;
  }
  
/*
'##::: ##::'#######::'########::'########:
 ###:: ##:'##.... ##: ##.... ##: ##.....::
 ####: ##: ##:::: ##: ##:::: ##: ##:::::::
 ## ## ##: ##:::: ##: ##:::: ##: ######:::
 ##. ####: ##:::: ##: ##:::: ##: ##...::::
 ##:. ###: ##:::: ##: ##:::: ##: ##:::::::
 ##::. ##:. #######:: ########:: ########:
..::::..:::.......:::........:::........::
*/
  void init_table_node(int s)
  {
    if (!labeled) {
#pragma omp parallel for schedule(static)
      for (int v = 0; v < num_verts_graph; ++v)
        dt.set_node(v, colors_g[v], 1.0);
    }
    else
    {
      int* labels_sub = part.get_labels(s);  
      int label_s = labels_sub[0];
#pragma omp parallel for schedule(static)
      for (int v = 0; v < num_verts_graph; ++v)
      {  
        int n = colors_g[v];
        int label_g = labels_g[v];
        if (label_g == label_s)
          dt.set(v, comb_num_indexes_set[s][n], 1.0);
      }
    }
  }
  

/*
:'######:::'#######::'##::::'##:'##::: ##:'########:
'##... ##:'##.... ##: ##:::: ##: ###:: ##:... ##..::
 ##:::..:: ##:::: ##: ##:::: ##: ####: ##:::: ##::::
 ##::::::: ##:::: ##: ##:::: ##: ## ## ##:::: ##::::
 ##::::::: ##:::: ##: ##:::: ##: ##. ####:::: ##::::
 ##::: ##: ##:::: ##: ##:::: ##: ##:. ###:::: ##::::
. ######::. #######::. #######:: ##::. ##:::: ##::::
:......::::.......::::.......:::..::::..:::::..:::::
*/ 
  float colorful_count_array_part(int s)
  {
    float cc = 0;
    int num_verts_sub = subtemplates[s].num_vertices();

    int active_index = part.get_active_index(s);
    int passive_index = part.get_passive_index(s);
    int num_verts_a = num_verts_table[active_index];  
    int num_verts_p = num_verts_table[passive_index];  
    int num_combinations = choose_table[num_verts_sub][num_verts_a];
    int num_colorsets_s = choose_table[num_colors][num_verts_sub];  
    int num_colorsets_a = choose_table[num_colors][num_verts_a];  
    int num_colorsets_p = choose_table[num_colors][num_verts_p];  
    unsigned long set_count_loop = 0;
    unsigned long total_count_loop = 0;
    unsigned long read_count_loop = 0; 

#pragma omp parallel
{   
    int *valid_nbrs = (int *) malloc(max_degree * sizeof(int));
        assert(valid_nbrs != NULL);
    int valid_nbrs_count = 0;

    float* counts_s = new float[num_colorsets_s];
    float* counts_v = new float[num_colorsets_a];
    float* counts_u = new float[num_colorsets_p];
    
#pragma omp for schedule(guided) reduction(+:cc) reduction(+:set_count_loop) \
        reduction(+:total_count_loop) reduction(+:read_count_loop)
    for (int v = begin_vert; v < end_vert; ++v)
    {
      int offset_vert = v-begin_vert;
      sizes_verts[offset_vert] = 0;
      valid_nbrs_count = 0;   
      int table_size_v = dt.table_size_a(v);
      if (table_size_v)
      {
        float* table_counts_v = dt.table_counts_a(v);
        short* table_colors_v = dt.colorsets_a(v);
        for (int i = 0; i < num_colorsets_a; ++i) counts_v[i] = 0.0;
        for (int i = 0; i < num_colorsets_s; ++i) counts_s[i] = 0.0;
        for (int i = 0; i < table_size_v; ++i)
        {     
          float count = table_counts_v[i];
          short color = table_colors_v[i];
          counts_v[color] = count;
        }        

        int* adjs = g->adjacent_vertices(v);
        int end = g->out_degree(v);
        ++read_count_loop;
       
        for (int i = 0; i < end; ++i) {
          int adj_i = adjs[i];
          if (dt.table_size_p(adj_i)) {
            valid_nbrs[valid_nbrs_count++] = adj_i;
          }
        }

        assert(valid_nbrs_count <= max_degree);        
        if (valid_nbrs_count)
        { 
          for (int vn = 0; vn < valid_nbrs_count; ++vn)
          {
            int u = valid_nbrs[vn];
            int table_size_u = dt.table_size_p(u);
            float* table_counts_u = dt.table_counts_p(u);
            short* table_colors_u = dt.colorsets_p(u);

            for (int i = 0; i < num_colorsets_p; ++i) counts_u[i] = 0.0;
            for (int i = 0; i < table_size_u; ++i)
            {
              float count = table_counts_u[i];
              short color = table_colors_u[i];
              counts_u[color] = count;
            }  

            for (int n = 0; n < num_colorsets_s; ++n)
            {
              float color_count = 0.0;                
              int* comb_indexes_a = comb_num_indexes[0][s][n];
              int* comb_indexes_p = comb_num_indexes[1][s][n];
              
              int p = num_combinations - 1;
              for (int a = 0; a < num_combinations; ++a, --p)
              {
                counts_s[n] += counts_v[comb_indexes_a[a]] * 
                      counts_u[comb_indexes_p[p]];
                ++read_count_loop;
              }
            }
          }

          for (int n = 0; n < num_colorsets_s; ++n)
          {
            float color_count = counts_s[n];
            if (color_count)
            {
              cc += color_count;
              ++set_count_loop;
              if (s != 0)
              {
                dt.set(v, comb_num_indexes_set[s][n], color_count);
                ++sizes_verts[offset_vert];  
              }
              else if (do_graphlet_freq || do_vert_output)
                final_vert_counts[offset_vert] += (double)color_count;
            }
            ++total_count_loop;  
          }
        }
      }
    }

    free(valid_nbrs);
    delete [] counts_s;
    delete [] counts_v;
    delete [] counts_u;
    
}    
#if DEBUG
    printf("%d done count fucks %d\n", rank, s);
#endif
    set_count = set_count_loop;
    total_count = total_count_loop;
    read_count = read_count_loop;

    dt.set_total_counts(set_count);

    return cc;
  }
  


/*
'########::::'###::::'########::'##:::::::'########::'######::
... ##..::::'## ##::: ##.... ##: ##::::::: ##.....::'##... ##:
::: ##:::::'##:. ##:: ##:::: ##: ##::::::: ##::::::: ##:::..::
::: ##::::'##:::. ##: ########:: ##::::::: ######:::. ######::
::: ##:::: #########: ##.... ##: ##::::::: ##...:::::..... ##:
::: ##:::: ##.... ##: ##:::: ##: ##::::::: ##:::::::'##::: ##:
::: ##:::: ##:::: ##: ########:: ########: ########:. ######::
:::..:::::..:::::..::........:::........::........:::......:::
*/
  void create_tables()
  {
    choose_table = init_choose_table(num_colors);
    create_num_verts_table();
    create_all_index_sets();
    create_all_color_sets();
    create_comb_num_system_indexes();
    delete_all_color_sets();
    delete_all_index_sets();
  }  
  
  void delete_tables()
  {
    for (int i = 0; i <= num_colors; ++i)
      delete [] choose_table[i];

    delete [] choose_table;
    delete_comb_num_system_indexes();
    delete [] num_verts_table;
  }
  
  void create_num_verts_table()
  {
    num_verts_table = new int[subtemplate_count];
    
    for (int s = 0; s < subtemplate_count; ++s)
      num_verts_table[s] = subtemplates[s].num_vertices();
  }
  
  void create_all_index_sets()
  {  
    index_sets = new int***[num_colors];
    
    for (int i = 0; i < (num_colors-1); ++i)
    {
      int num_vals = i + 2;
      index_sets[i] = new int**[(num_vals-1)];
      for (int j = 0; j < (num_vals-1); ++j)
      {  
        int set_size = j + 1;
        int num_combinations = choose(num_vals, set_size);
        index_sets[i][j] = new int*[num_combinations];
        
        int* set = init_permutation(set_size);
        
        for (int k = 0; k < num_combinations; ++k)
        {
          index_sets[i][j][k] = new int[set_size];
          for (int p = 0; p < set_size; ++p)
          {
            index_sets[i][j][k][p] = set[p] - 1;
          }

          next_set(set, set_size, num_vals);
        }
          
        delete [] set;
      }
    }
  }
  
  void delete_all_index_sets()
  {
    for (int i = 0; i < (num_colors-1); ++i)
    {
      int num_vals = i + 2;      
      for (int j = 0; j < (num_vals-1); ++j)
      {  
        int set_size = j + 1;
        int num_combinations = choose(num_vals, set_size);        
        for (int k = 0; k < num_combinations; ++k)
        {
          delete [] index_sets[i][j][k];
        }          
        delete [] index_sets[i][j];
      }      
      delete [] index_sets[i];
    }    
    delete [] index_sets;
  }  
  
  void create_all_color_sets()
  {
    color_sets = new int****[subtemplate_count];
    
    //printf("%d a %d\n",rank,subtemplate_count);
    for (int s = 0; s < subtemplate_count; ++s)
    {
      int num_verts_sub = subtemplates[s].num_vertices();
      //printf("%d b %d %d\n",rank,s,num_verts_sub);
      if (num_verts_sub > 1)
      {      
        int num_sets = choose(num_colors, num_verts_sub);
        color_sets[s] = new int***[num_sets];        
        //printf("%d c %d\n",rank,num_sets);
        int* colorset = init_permutation(num_verts_sub);
        for (int n = 0; n < num_sets; ++n)
        {
          int num_child_combs = num_verts_sub - 1;
          color_sets[s][n] = new int**[num_child_combs];
          //printf("%d d %d %d\n",rank,n, num_child_combs);
          for (int c = 0; c < num_child_combs; ++c)
          {
            int num_verts_1 = c+1;
            int num_verts_2 = num_verts_sub - num_verts_1;
            int** index_set_1 = index_sets[(num_verts_sub-2)][(num_verts_1-1)];
            int** index_set_2 = index_sets[(num_verts_sub-2)][(num_verts_2-1)];
        
            int num_child_sets = choose(num_verts_sub, (c + 1));
            color_sets[s][n][c] = new int*[num_child_sets];
            //printf("%d e %d %d %d %d\n",rank,c,num_verts_1,num_verts_2,num_child_sets);
            for (int i = 0; i < num_child_sets; ++i)
            {
              color_sets[s][n][c][i] = new int[num_verts_sub];
              
              for (int j = 0; j < num_verts_1; ++j)
                color_sets[s][n][c][i][j] = colorset[index_set_1[i][j]];
              for (int j = 0; j < num_verts_2; ++j)
                color_sets[s][n][c][i][j+num_verts_1] = colorset[index_set_2[i][j]];
            }
          }
          next_set(colorset, num_verts_sub, num_colors);
        }
        delete [] colorset;
      }      
    }
  }
  
  void delete_all_color_sets()
  {
    for (int s = 0; s < subtemplate_count; ++s)
    {
      int num_verts_sub = subtemplates[s].num_vertices();      
      if (num_verts_sub > 1)
      {      
        int num_sets = choose(num_colors, num_verts_sub);        
        for (int n = 0; n < num_sets; ++n)
        {
          int num_child_combs = num_verts_sub - 1;
          for (int c = 0; c < num_child_combs; ++c)
          {
            int num_child_sets = choose(num_verts_sub, (c + 1));
            for (int i = 0; i < num_child_sets; ++i)
            {
              delete [] color_sets[s][n][c][i];              
            }            
            delete [] color_sets[s][n][c];            
          }
          delete [] color_sets[s][n];
        }          
        delete [] color_sets[s];
      }      
    }    
    delete [] color_sets;
  }
  
  void create_comb_num_system_indexes()
  {
    comb_num_indexes = new int***[2];
    comb_num_indexes[0] = new int**[subtemplate_count];
    comb_num_indexes[1] = new int**[subtemplate_count];    
    comb_num_indexes_set = new int*[subtemplate_count];
    
    for (int s = 0; s < subtemplate_count; ++s)
    {
      int num_verts_sub = subtemplates[s].num_vertices();
      int num_combinations_s = choose(num_colors, num_verts_sub);
      
      if (num_verts_sub > 1)
      {  
        comb_num_indexes[0][s] = new int*[num_combinations_s];
        comb_num_indexes[1][s] = new int*[num_combinations_s];
      }
      comb_num_indexes_set[s] = new int[num_combinations_s];
      int* colorset_set = init_permutation(num_verts_sub);
      
      for (int n = 0; n < num_combinations_s; ++n)
      {      
        comb_num_indexes_set[s][n] = get_color_index(colorset_set, num_verts_sub);
      
        if (num_verts_sub > 1)
        {  
          int num_verts_a = part.get_num_verts_active(s);
          int num_verts_p = part.get_num_verts_passive(s);          
          int active_index = part.get_active_index(s);
          int passive_index = part.get_passive_index(s);
      
          int* colors_a;        
          int* colors_p;
          int** colorsets = color_sets[s][n][num_verts_a - 1];          
      
          int num_combinations_a = choose(num_verts_sub, num_verts_a);
          comb_num_indexes[0][s][n] = new int[num_combinations_a];        
          comb_num_indexes[1][s][n] = new int[num_combinations_a];        
          
          int p = num_combinations_a - 1;
          for (int a = 0; a < num_combinations_a; ++a, --p)
          {  
            colors_a = colorsets[a];          
            colors_p = colorsets[p] + num_verts_a;
            
            int color_index_a = get_color_index(colors_a, num_verts_a);
            int color_index_p = get_color_index(colors_p, num_verts_p);  

            comb_num_indexes[0][s][n][a] = color_index_a;
            comb_num_indexes[1][s][n][p] = color_index_p;
          }
        }        
        next_set(colorset_set, num_verts_sub, num_colors);
      }
      delete [] colorset_set;
    }
  }
  
  void delete_comb_num_system_indexes()
  {
    for (int s = 0; s < subtemplate_count; ++s)
    {
      int num_verts_sub = subtemplates[s].num_vertices();      
      int num_combinations_s = choose(num_colors, num_verts_sub);
      
      for (int n = 0; n < num_combinations_s; ++n)
      {  
        if (num_verts_sub > 1)
        {  
          delete [] comb_num_indexes[0][s][n];        
          delete [] comb_num_indexes[1][s][n];
        }
      }
      
      if (num_verts_sub > 1)
      {  
        delete [] comb_num_indexes[0][s];
        delete [] comb_num_indexes[1][s];
      }
      
      delete [] comb_num_indexes_set[s];
    }
    
    delete [] comb_num_indexes[0];
    delete [] comb_num_indexes[1];
    delete [] comb_num_indexes; 
    delete [] comb_num_indexes_set;
  }    
  
  Graph* g;
  Graph* t;
  int* colors_g;  
  int* labels_g;  
  int* labels_t;
  bool labeled;
    
  Graph* subtemplates;
  int subtemplate_count;
  int num_colors;
  int num_iter;
  
  dynamic_table_part dt;
  partitioner part;
  
  int** choose_table;
  int**** index_sets;
  int***** color_sets;
  int**** comb_num_indexes;
  int** comb_num_indexes_set;
  int* num_verts_table;  
  int num_verts_graph;
  int max_degree;
  
  unsigned long set_count;
  unsigned long total_count;
  unsigned long read_count;

  double* final_vert_counts;
  bool do_graphlet_freq;
  bool do_vert_output;
  bool calculate_automorphisms;

  double transfer_time_sum;
  double transfer_size_sum;

  int begin_vert;
  int end_vert;
  int num_verts_part;

  int* vert_offsets;

  unsigned long* comm_sizes_send;
  unsigned long* comm_sizes_rec;
  unsigned long* total_send;
  unsigned long* total_recv;
  int* sizes_verts;
  int* comm_num_rec;
  int* comm_num_send;
  int** comm_verts_rec;
  int** comm_verts_send;

  int* part_offsets;

};




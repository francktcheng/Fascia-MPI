
using namespace std;


class dynamic_table_part {
public:
  dynamic_table_part()
  { 
    table = NULL;
  }
  
  void init_part(Graph* subs, int num_subtemplates, 
    int num_vertices, int num_cols, int Begin_vert, int End_vert,
    int Rank, int Nprocs)
  {  
#if DEBUG
    printf("%d initing part %d %d\n", Rank, Begin_vert, End_vert);
#endif
    subtemplates = subs;
    num_subs = num_subtemplates;
    num_verts = num_vertices;
    begin_vert = Begin_vert;
    end_vert = End_vert;
    num_verts_part = end_vert - begin_vert;
    num_colors = num_cols;  
    rank = Rank;
    nprocs = Nprocs;
    init_choose_table();
    init_num_colorsets();
#if DEBUG
    printf("%d done initing tables %d\n", rank, num_verts_part);
#endif
    colorsets = (short**) malloc(num_subs * sizeof(short*));
    vert_offsets = (unsigned long**)malloc(num_subs* sizeof(unsigned long*));   
    table = (float **) malloc(num_subs * sizeof(float *));
    assert(colorsets != NULL);
    assert(vert_offsets != NULL);
    assert(table != NULL);
#if DEBUG
    printf("%d done malloc1\n", rank);
#endif
    is_sub_inited = (bool*) malloc(num_subs * sizeof(bool));
    sub_sizes = (unsigned long*) malloc(num_subs * sizeof(unsigned long));
    sub_sizes_part = (unsigned long*) malloc(num_subs * sizeof(unsigned long));
    assert(is_sub_inited != NULL);
    assert(sub_sizes != NULL);
    assert(sub_sizes_part != NULL);
#if DEBUG
    printf("%d done malloc2\n", rank);
#endif
    for (int s = 0; s < num_subs; ++s)
    {
      is_sub_inited[s] = false;
      sub_sizes[s] = 0;
      sub_sizes_part[s] = 0;
    }
#if DEBUG
    printf("%d done init\n", rank);
#endif
  }

  int init_sub(int s)
  {
    cur_sub = s;
    int num_verts_sub = subtemplates[s].num_vertices();
    assert(num_verts_sub == 1);
#if DEBUG
    printf("%d initing %d %d\n", rank, s, num_verts_part);
#endif
    int num_colorsets = choose_table[num_colors][num_verts_sub];
    int num_colorsets_a = 0;
    int num_colorsets_p = 0;

    unsigned long num_total = (unsigned long)num_verts * (unsigned long)num_colorsets;
#if DEBUG
    printf("%d initing %d, num tot %lu, cols %d\n", rank, s, num_total, num_colorsets);
#endif
    table[s] = (float*) malloc(num_total * sizeof(float));
    assert(table[s] != NULL);
#if DEBUG
    printf("%d malloc success\n", rank);
#endif
    cur_table = table[s];
    cur_table_a = NULL;
    cur_table_p = NULL;
    cur_num_colorsets = num_colorsets;
    cur_num_colorsets_a = num_colorsets_a;
    cur_num_colorsets_p = num_colorsets_p;
    sub_sizes[s] = num_total;
    sub_sizes_part[s] = num_total;

    for (unsigned i = 0; i < num_total; ++i)
      cur_table[i] = 0.0;

    is_sub_inited[s] = true;
#if DEBUG
    printf("%d done sub init\n", rank);
#endif
    return 0; 
  }
  
  int init_sub(int s, int a, int p)
  {
    assert(a != NULL_VAL);
    assert(p != NULL_VAL);
    cur_sub = s;
    cur_a = a;
    cur_p = p;

    int num_verts_a = subtemplates[a].num_vertices();
    int num_verts_p = subtemplates[p].num_vertices();
    if (num_verts_a == 1)
      a = table_node;
    if (num_verts_p == 1)
      p = table_node;

    cur_num_colorsets = choose_table[num_colors][subtemplates[s].num_vertices()];
    cur_num_colorsets_a = choose_table[num_colors][num_verts_a];
    cur_num_colorsets_p = choose_table[num_colors][num_verts_p];

    unsigned long num_total = 
      (unsigned long)num_verts * (unsigned long)cur_num_colorsets;
    sub_sizes[s] = num_total;

    unsigned long num_total_part = 
      (unsigned long)num_verts_part * (unsigned long)cur_num_colorsets;
    sub_sizes_part[s] = num_total_part;

    if (s)
    {
      table[s] = (float*) malloc(num_total_part * sizeof(float));
      assert(table[s] != NULL);
    }

    cur_table = table[s];
    cur_table_a = table[a];
    cur_table_p = table[p];

    if (s)
    {      
      for (unsigned i = 0; i < num_total_part; ++i)
        table[s][i] = 0.0; 
    }
  
    is_sub_inited[s] = true;
    cur_colorsets_a = colorsets[a];
    cur_colorsets_p = colorsets[p];
    cur_vert_offsets_a = vert_offsets[a];
    cur_vert_offsets_p = vert_offsets[p];

    return 0;   
  }

  void set_table_node(int s)
  {
    init_sub(s); 
    table_node = s;

#if DEBUG
    printf("%d set node\n", rank);
#endif
    unsigned long table_size = sub_sizes[s];
    colorsets[s] = (short*) malloc(table_size * sizeof(short));
    vert_offsets[s] = (unsigned long*)malloc((num_verts+1) * sizeof(unsigned long));
    assert(colorsets[s] != NULL);
    assert(vert_offsets[s] != NULL);
#if DEBUG
    printf("%d done other arrs\n", rank);
#endif
    unsigned long cur_offset = 0;
    for (int i = 0; i < num_verts; ++i)
    {
      vert_offsets[s][i] = cur_offset;
      for (short j = 0; j < num_colors; ++j)
        colorsets[s][cur_offset++] = j;
    }    
    vert_offsets[s][num_verts] = cur_offset;
    assert(table_size == cur_offset);
#if DEBUG
    printf("%d done node\n", rank);
#endif
  }

  void set_sub_to_table_node(int s)
  {
#if DEBUG
    printf("%d setting sub to table node %d %d\n", rank, s, table_node);
#endif
    table[s] = table[table_node];
    vert_offsets[s] = vert_offsets[table_node];
    colorsets[s] = colorsets[table_node];
    sub_sizes[s] = sub_sizes[table_node];
    sub_sizes_part[s] = sub_sizes_part[table_node];
  }

  void init_comp_table(unsigned long* part_counts)
  {
#if DEBUG
    printf("%d initing comp table\n", rank);
#endif
    unsigned long table_size = 0;
    for (int i = 0; i < nprocs; ++i)
    {
#if DEBUG
    printf("%d ts %lu, %d, %lu\n", rank, table_size, i, part_counts[i]);
#endif
      table_size += part_counts[i];
    }
#if DEBUG
    printf("%d TS %lu\n", rank, table_size);
#endif
    temp_colorsets = (short*)malloc(table_size * sizeof(short));
    temp_table = (float*)malloc(table_size * sizeof(float));
    temp_vert_offsets = (unsigned long*)malloc((num_verts+1)*sizeof(unsigned long));

    assert(temp_table != NULL);
    assert(temp_colorsets != NULL);
    assert(temp_vert_offsets != NULL);

    running_offset = 0;
#if DEBUG
    printf("%d done initing comp table\n", rank);
#endif
  }

  void compress_for_send(unsigned long tot_count, 
    int num_vert_send, int* comm_verts_send, int send_to,
    float* counts_comp, short* colorsets_comp, long unsigned* offsets_comp)
  {
#if DEBUG
    printf("%d compressing comp table %lu, %d\n", rank, tot_count, num_vert_send);
#if DEBUG1
    for (int i = 0; i < num_vert_send; ++i)
      printf("comm_verts_send %d\n", comm_verts_send[i]);
#endif
#endif
    long unsigned cur_offset = 0;
    for (int i = 0; i < num_vert_send; ++i)
    {      
      int vert = comm_verts_send[i];
      offsets_comp[i] = cur_offset;
      for (short j = 0; j < cur_num_colorsets; ++j)
      {
        float count = cur_table[(vert - begin_vert) * cur_num_colorsets + j];
        if (count > 0.0)
        {
          counts_comp[cur_offset] = count;
          colorsets_comp[cur_offset] = j;
          ++cur_offset;
        }
      }
    }
#if DEBUG
    printf("%d off %lu, %lu\n", rank, cur_offset, tot_count);
#endif   
    offsets_comp[num_vert_send] = cur_offset;
    assert(cur_offset == tot_count);
#if DEBUG
    printf("%d done compressing comp table\n", rank);
#endif
  }

  void append_to_table(unsigned long count, int num_vert_rec, 
    int* part_offsets, int from_rank, int* comm_verts_rec,
    float* counts_comp, short* colorsets_comp, unsigned long* offsets_comp)
  {
#if DEBUG
    printf("%d appended comp table\n", rank);
#endif
    copy(counts_comp, counts_comp+count, &temp_table[running_offset]);
    copy(colorsets_comp, colorsets_comp+count, &temp_colorsets[running_offset]);

    int cur_vert = 0;
    for (int v = part_offsets[from_rank]; v < part_offsets[from_rank+1]; ++v)
    {
      temp_vert_offsets[v] = offsets_comp[cur_vert] + running_offset;      
      if (cur_vert < num_vert_rec-1)
        if (comm_verts_rec[cur_vert] == v)
          ++cur_vert;
      //printf("%d %d %d %d - %d %d - %d\n", v, cur_vert, num_vert_rec, comm_verts_rec[cur_vert], rank, from_rank, offsets_comp[cur_vert]);
    }
    running_offset += count;
#if DEBUG
    printf("%lu\n", running_offset);
    printf("%d done appending comp table\n", rank);
#endif
  }

  void finalize()
  {
#if DEBUG
    printf("%d finalizing comp table, clearing\n", rank);
    printf("%lu\n", running_offset);
#endif
    temp_vert_offsets[num_verts] = running_offset;
    free(table[cur_sub]);
    table[cur_sub] = temp_table;
    colorsets[cur_sub] = temp_colorsets;
    vert_offsets[cur_sub] = temp_vert_offsets;
#if DEBUG
    printf("%d done finalizing comp table\n", rank);
#endif
  }

  void clear_sub(int s)
  {
#if DEBUG
    printf("%d clearing sub\n", rank);
#endif
    free(table[s]);
    free(colorsets[s]);
    free(vert_offsets[s]);
#if DEBUG
    printf("%d done clearing sub\n", rank);
#endif
  }  
   
  void clear_table()
  {
#if DEBUG
    printf("%d clearing table\n", rank);
#endif
    for (int s = 0; s < num_subs; s++)
    {
      if (is_sub_inited[s]) 
        free(table[s]);
    }

    free(table);
    free(is_sub_inited);
    free(is_vert_inited);    
  } 

  int table_size_a(int v)
  {
    assert(v < num_verts);
    assert(v >= 0);
    return (cur_vert_offsets_a[v+1] - cur_vert_offsets_a[v]); 
  }
  int table_size_p(int v)
  {
#if DEBUG
    //printf("%d ts %d %d %d\n", rank, v, 
    //  cur_vert_offsets_p[v+1], cur_vert_offsets_p[v]);
#endif
    assert(v < num_verts);
    assert(v >= 0);
    return (cur_vert_offsets_p[v+1] - cur_vert_offsets_p[v]); 
  }

  float* table_counts_a(int v)
  {    
    assert(v < num_verts);
    assert(v >= 0);
    unsigned long index = cur_vert_offsets_a[v];
    assert(index < sub_sizes[cur_a]);
    return &cur_table_a[cur_vert_offsets_a[v]];
  }
  float* table_counts_p(int v)
  {
    assert(v < num_verts);
    assert(v >= 0);
    unsigned long index = cur_vert_offsets_p[v];
    assert(index < sub_sizes[cur_p]);
    return &cur_table_p[cur_vert_offsets_p[v]];
  }

  short* colorsets_a(int v)
  {
    return &cur_colorsets_a[cur_vert_offsets_a[v]];
  }
  short* colorsets_p(int v)
  {
    return &cur_colorsets_p[cur_vert_offsets_p[v]];
  }


  void set(int vertex, int comb_num_index, float count)
  {
    unsigned long index = (unsigned long)(vertex-begin_vert);
    index *= (unsigned long)cur_num_colorsets;
    index += (unsigned long)comb_num_index;
    assert(index < sub_sizes[cur_sub]);
    cur_table[index] = count;
  }

  void set_node(int vertex, int comb_num_index, float count)
  {
    unsigned long index = (unsigned long)vertex;
    index *= (unsigned long)cur_num_colorsets;
    index += (unsigned long)comb_num_index;
    assert(index < sub_sizes[cur_sub]);
    cur_table[index] = count;
  }

  void set_total_counts(int total_count)
  { cur_total_count = total_count; }

  unsigned long get_total_counts()
  { return cur_total_count; }

  void print_table(int outnum)
  {
    char* outfile = new char[20];
    sprintf(outfile,"%d",outnum);

    ofstream out;
    out.open(outfile);

    int num_colorsets = choose_table[num_colors][subtemplates[cur_sub].num_vertices()];
    unsigned long cur_offset = 0;
    for (int v = begin_vert; v < end_vert; ++v)
    {
      out << v << ": ";
      for (int i = 0; i < num_colorsets; ++i)
        out << cur_table[cur_offset++] << " ";
      out << endl;
    }
    out.close();
    delete [] outfile;
  }
   
private:
  void init_choose_table()
  {
    choose_table = new int*[num_colors + 1];
  
    for (int i = 0; i <= num_colors; ++i)
      choose_table[i] = new int[num_colors + 1];
    
    for (int i = 0; i <= num_colors; ++i)
      for (int j = 0; j <= num_colors; ++j)
        choose_table[i][j] = choose(i, j);    
  }
  
  void init_num_colorsets()
  {
    num_colorsets = new int[num_subs];    
    for (int s = 0; s < num_subs; ++s)
      num_colorsets[s] = choose(num_colors, subtemplates[s].num_vertices());
  }

  int num_subs;
  int num_verts;
  int end_vert;
  int begin_vert;
  int num_colors;


  int** choose_table;
  int* num_colorsets;
  
  Graph* subtemplates;
  float** table;
  short** colorsets;
  unsigned long** vert_offsets;

  bool* is_vert_inited;
  bool* is_sub_inited;

  int cur_sub;
  int cur_a;
  int cur_p;
  float* cur_table;
  float* cur_table_a;
  float* cur_table_p;
  unsigned long* cur_vert_offsets_a;
  unsigned long* cur_vert_offsets_p;
  short* cur_colorsets_a;
  short* cur_colorsets_p;
  int cur_num_colorsets;
  int cur_num_colorsets_a;
  int cur_num_colorsets_p;
  unsigned long* sub_sizes;
  unsigned long* sub_sizes_part;

  float* temp_table;
  short* temp_colorsets;
  long unsigned* temp_vert_offsets;

  int rank;
  int nprocs;
  int num_verts_part;
  unsigned long cur_total_count;
  int table_node;

  unsigned long running_offset;

};

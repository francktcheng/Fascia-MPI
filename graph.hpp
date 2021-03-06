// Copyright (c) 2013, The Pennsylvania State University.
// All rights reserved.
// 
// See COPYING for license.

class Graph{
public:
  Graph() {};
  ~Graph() {};

  void init(int n, unsigned m, int* srcs, int* dsts)
  {
    num_verts = n;
    num_edgs = 2*m;
    max_deg = 0;
    adjacency_array = new int[2*m];
    degree_list = new unsigned[n+1];
    degree_list[0] = 0;

    unsigned* temp_deg_list = new unsigned[n];
    for (int i = 0; i < n; ++i)
      temp_deg_list[i] = 0;
    for (unsigned i = 0; i < m; ++i)
    {
      temp_deg_list[srcs[i]]++;
      temp_deg_list[dsts[i]]++;
    }
    for (int i = 0; i < n; ++i)
      max_deg = temp_deg_list[i] > max_deg ? temp_deg_list[i] : max_deg;
    for (int i = 0; i < n; ++i)
      degree_list[i+1] = degree_list[i] + temp_deg_list[i];
    copy(degree_list, degree_list+n, temp_deg_list);
    for (unsigned i = 0; i < m; ++i)
    {
      adjacency_array[temp_deg_list[srcs[i]]++] = dsts[i];
      adjacency_array[temp_deg_list[dsts[i]]++] = srcs[i];
    }

    delete [] temp_deg_list;
  }
  
  int* adjacent_vertices(int v)
  {
    return (&adjacency_array[degree_list[v]]);
  }
  
  unsigned out_degree(int v)
  {
    return degree_list[v+1] - degree_list[v];
  }
  
  int* adjacencies() const
  {
    return adjacency_array;
  }
  
  unsigned* degrees() const
  {
    return degree_list;
  }
  
  int num_vertices() const
  {
    return num_verts;
  }
  
  unsigned num_edges() const
  {
    return num_edgs;
  }

  unsigned max_degree() const
  {
    return max_deg;
  }
  
  Graph& operator= (const Graph& param)
  {
    num_verts = param.num_vertices();
    num_edgs = param.num_edges();
    max_deg = param.max_degree();
    
    adjacency_array = new int[2*num_edgs];
    degree_list = new unsigned[num_verts+1];    
    
    copy(param.adjacencies(), param.adjacencies() + num_edgs, adjacency_array);
    copy(param.degrees(), param.degrees() + (num_verts+1), degree_list);
      
    return *this;
  }
  
  void clear()
  {
    delete [] adjacency_array;
    delete [] degree_list;
  }
  
private:
  int num_verts;
  unsigned num_edgs;
  int max_deg;

  int* adjacency_array;
  unsigned* degree_list;

};

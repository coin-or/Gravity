#include <stdio.h>
#ifdef SYSV
#include <string.h>
#else
#include <string.h>
#endif
#undef min
typedef union {
  struct vertex_struct *V;
  struct arc_struct *A;
  struct graph_struct *G;
  char *S;
  long I;
} util;
typedef struct vertex_struct {
  struct arc_struct *arcs;
  char *name;
  util u, v, w, x, y, z;
} Vertex;
typedef struct arc_struct {
  struct vertex_struct *tip;
  struct arc_struct *next;
  long len;
  util a, b;
} Arc;
#define init_area(s)  *s= NULL
struct area_pointers {
  char *first;
  struct area_pointers *next;
};
typedef struct area_pointers *Area[1];
#define ID_FIELD_SIZE 161
typedef struct graph_struct {
  Vertex *vertices;
  long n;
  long m;
  char id[ID_FIELD_SIZE];
  char util_types[15];
  Area data;
  Area aux_data;
  util uu, vv, ww, xx, yy, zz;
} Graph;
typedef unsigned long siz_t;
extern long verbose;
extern long panic_code;
#define alloc_fault (-1)
#define no_room 1
#define early_data_fault 10
#define late_data_fault 11
#define syntax_error 20
#define bad_specs 30
#define very_bad_specs 40
#define missing_operand 50
#define invalid_operand 60
#define impossible 90
extern long gb_trouble_code;
extern char *gb_alloc ();
#define gb_typed_alloc(n,t,s) \
               (t*)gb_alloc((long)((n)*sizeof(t)),s)
extern void gb_free ();
#define n_1  uu.I
#define mark_bipartite(g,n1) g->n_1= n1,g->util_types[8]= 'I'
extern long extra_n;
extern char null_string[];
extern void make_compound_id ();
extern void make_double_compound_id ();
extern siz_t edge_trick;
#define gb_new_graph gb_nugraph
#define gb_new_arc gb_nuarc
#define gb_new_edge gb_nuedge
extern Graph *gb_new_graph ();
extern void gb_new_arc ();
extern Arc *gb_virgin_arc ();
extern void gb_new_edge ();
extern char *gb_save_string ();
extern void switch_to_graph ();
extern void gb_recycle ();
extern void hash_in ();
extern Vertex *hash_out ();
extern void hash_setup ();
extern Vertex *hash_lookup ();
extern Graph *board ();
extern Graph *simplex ();
extern Graph *subsets ();
extern Graph *perms ();
extern Graph *parts ();
extern Graph *binary ();
extern Graph *complement ();
extern Graph *gunion ();
extern Graph *intersection ();
extern Graph *lines ();
extern Graph *product ();
extern Graph *induced ();
#define complete(n) board((long)(n),0L,0L,0L,-1L,0L,0L)
#define transitive(n) board((long)(n),0L,0L,0L,-1L,0L,1L)
#define empty(n) board((long)(n),0L,0L,0L,2L,0L,0L)
#define circuit(n) board((long)(n),0L,0L,0L,1L,1L,0L)
#define cycle(n) board((long)(n),0L,0L,0L,1L,1L,1L)
#define disjoint_subsets(n,k) subsets((long)(k),1L,(long)(1-(n)),0L,0L,0L,1L,0L)
#define petersen() disjoint_subsets(5,2)
#define all_perms(n,directed) perms((long)(1-(n)),0L,0L,0L,0L,0L,\
   (long)(directed))
#define all_parts(n,directed) parts((long)(n),0L,0L,(long)(directed))
#define all_trees(n,directed) binary((long)(n),0L,(long)(directed))
#define cartesian 0
#define direct 1
#define strong 2
#define ind z.I
#define IND_GRAPH 1000000000
#define subst y.G
extern Graph *bi_complete ();
extern Graph *wheel ();
#define random_graph r_graph
#define random_bigraph r_bigraph
#define random_lengths r_lengths
extern Graph *random_graph ();
extern Graph *random_bigraph ();
extern long random_lengths ();
extern long *gb_fptr;
extern long gb_flip_cycle ();
#define gb_next_rand() (*gb_fptr>=0?*gb_fptr--:gb_flip_cycle())
#define max_stack 10
#define is(a) !strcmp(argv[argp],a)

typedef struct {
  long            i, j, k;
}               Face;
typedef struct {
  long            u, v;
}               Edge;
Graph          *
planar(n, density, seed)
  long            n;
  double          density;
  long            seed;
{
  Graph          *new_graph;
  Face           *faces;
  Edge           *edges;
  long            w, f, a, b, c, d, e;
  long            n_edges, max_edges;
  if (n <= 0)
    return NULL;
  n_edges = 0;
  if (n < 3) {
    if (n == 2) {
      max_edges = 1;
      edges = (Edge *) malloc(max_edges * sizeof(Edge));
      edges[n_edges].u = 0;
      edges[n_edges++].v = 1;
    }
  } else {
    gb_init_rand(seed);
    max_edges = (3 * (n - 2));
    edges = (Edge *) malloc(max_edges * sizeof(Edge));
    faces = (Face *) malloc((2 * n - 5) * sizeof(Face));
    faces[0].i = 0;
    faces[0].j = 1;
    faces[0].k = 2;
    edges[n_edges].u = 0;
    edges[n_edges++].v = 1;
    edges[n_edges].u = 0;
    edges[n_edges++].v = 2;
    edges[n_edges].u = 1;
    edges[n_edges++].v = 2;
    for (w = 1; w <= n - 3; w++) {
      f = gb_unif_rand(2 * w - 1);
      a = faces[f].i;
      b = faces[f].j;
      c = faces[f].k;
      d = w + 2;
      edges[n_edges].u = a;
      edges[n_edges++].v = d;
      edges[n_edges].u = b;
      edges[n_edges++].v = d;
      edges[n_edges].u = c;
      edges[n_edges++].v = d;
      faces[f].k = d;
      faces[2 * w - 1].i = a;
      faces[2 * w - 1].j = c;
      faces[2 * w - 1].k = d;
      faces[2 * w].i = b;
      faces[2 * w].j = c;
      faces[2 * w].k = d;
    }
    free((char *) faces);
  }
  while (n_edges > (density * max_edges) / 100.0) {
    e = gb_unif_rand(n_edges);
    edges[e] = edges[--n_edges];
  }
  new_graph = gb_new_graph(n);
  if (new_graph != NULL)
    for (e = 0; e < n_edges; e++)
      gb_new_edge(new_graph->vertices + edges[e].u,
		  new_graph->vertices + edges[e].v, 1L);
  return new_graph;
}
#include <math.h>
double 
gauss()
{
  double          u1, u2, v1, v2, s;
  static double   max_int = (double) 0x7fffffff;
  static double   x = 0.0;
  static          left = 0;
  if (left) {
    left = 0;
    return x;
  } else {
    do {
      u1 = (double) gb_next_rand() / max_int;
      u2 = (double) gb_next_rand() / max_int;
      v1 = 2.0 * u1 - 1.0;
      v2 = 2.0 * u2 - 1.0;
      s = v1 * v1 + v2 * v2;
    } while (s >= 1.0);
    x = v2 * sqrt((-2 * log(s)) / s);
    left = 1;
    return v1 * sqrt((-2 * log(s)) / s);
  }
}
int 
main(argc, argv)
  int             argc;
  char          **argv;
{
  Graph          *g[max_stack], *gtmp, *gtmp1;
  long            stack = -1;
  long            lower_weight, upper_weight;
  long            seed;
  long            length, size;
  long            board_height, board_width, board_size, space_dimension,
                  move_type;
  long            dimension, sum, bound;
  long            scalar;
  long            argp = 0;
  double          density;
  if (argc == 1) {
    printf("\n\
              Rudy: a rudimental graph generator by JRT\n\
              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
use :=   rudy <graph_expression>  [outputs the graph represented by\n\
                                   <graph_expression> to the standard\n\
                                   output.]\n");
    printf("\n\
<graph_expression> [is an expression in Reverse Polish Notation that uses\n\
        <simple_graph>'s as operands and has <postfix_unary_operator>'s and\n\
        <postfix_binary_operator>'s. The maximum stack size is 10.\n\
        A scalar parameter has to be integer or real depending on if its\n\
        name has a lower or an upper case initial, respectively.]\n");
    printf("\n\
<postfix_unary_operator> := -line |\n\
                            -random <lower_weight> <upper_weight> <seed> |\n\
                            -complement | \n\
                            -times <scalar> |\n\
                            -plus <scalar>\n");
    printf("\n\
-line   [produces the line graph of the operand. Weights are set to 1.]\n\
-random <lower_weight> <upper_weight> <seed> [replaces the edge weights of\n\
        the operand with integer weights drawn from a uniform distribution\n\
        with range [<low>..<high>]. The random generator is initialized\n\
        by <seed>.]\n\
-complement  [complements the graph. Weights are set to 1.]\n\
-times <scalar> [Multiplies all weights by <scalar>.]\n\
-plus <scalar> [Adds <scalar> to all weights.]\n");
    printf("\n\
<postfix_binary_operator> := + | x | :\n");
    printf("\n\
        [denote the two operands by G1=(V1,E1) and G2=(V2,E2) and the\n\
         result by G=(V,E).]\n\
x       [V is the cartesian product of V1 and V2, i.e., is the set of\n\
         all the ordered pairs (v1,v2) with v1 in V1 and v2 in V2. The edge\n\
         from (u1,u2) to (v1,u2) is in E whenever there's an edge from u1 \n\
         to v1 in E1; the edge from (u1,u2) to (u1,v2) is in E whenever \n\
         there's an edge from u2 to v2 in E2 (cartesian graph product).]\n\
+       [it is assumed that V1=V2; then V=V1=V2, E is the union of E1 and E2\n\
         (graph union).]\n\
:       [V1 and V2 are cosidered as disjoint sets; then V is the union of V1\n\
         and V2 while E is the union of E1, E2, and the edge set of the \n\
         complete bipartite graph K(V1,V2)\n");
    printf("\n\
<simple_graph> := -leap_2D <board_height> <board_width> <move_type> |\n\
                  -leap <board_size> <space_dimension> <move_type> |\n\
                  -wrapped_leap_2D <board_height> <board_width> <move_type> |\n\
                  -wrapped_leap <board_size> <space_dimension> <move_type> |\n\
                  -grid_2D <board_height> <board_width> |\n\
                  -grid <board_size> <space_dimension> |\n\
                  -toroidal_grid_2D <board_height> <board_width> |\n\
                  -toroidal_grid <board_size> <space_dimension> |\n\
                  -simplex <sum> <dimension> |\n\
                  -bounded_simplex <sum> <dimension> <bound> |\n\
                  -circuit <length> |\n\
                  -clique <size> |\n\
                  -planar <size> <Density> <seed>\n");
    printf("\
                  -spinglass2pm <n_rows> <n_columns> <%% of - bonds> <seed>\n\
                  -spinglass3pm <n_rows> <n_columns> <n_layers> \n\
                                <%% of - bonds> <seed>\n\
                  -spinglass2g <n_rows> <n_columns> <seed>\n\
                  -spinglass3g <n_rows> <n_columns> <n_layers> <seed>\n\
                  -rnd_graph <size> <Density> <seed>\n");
    printf("\n\
-leap_2D [generates a leap graph on a bidimensional board of size\n\
         <board_height> x <board_width>. <move_type> is the sum of\n\
         the square of the increments of each of the coordinates defining\n\
         a legal minimal move. Example from a chessboard: 1 for a root, 2\n\
         for a bishop, for a knight move (all in 2 dimensions). A node of\n\
         the graph is associated with a cell of the board. Two nodes are\n\
         adiacent if there exists a sequence of legal moves that brings\n\
         from one of the corresponding two cells to the other. The edge\n\
         weight is the minimal number of such moves.]\n");
    printf("\n\
-leap    [it is like -leap_2D, but it generates a leap graph on a \n\
         <space_dimension>-dimensional cubic board of size <board_size>.]\n");
    printf("\n\
-wrapped_leap_2D [Is like -leap_2D but the cell coordinates are computed\n\
         modulo the board dimensions. Weights are all 1].\n");
    printf("\n\
-wrapped_leap [Is like -leap but the cell coordinates are computed\n\
         modulo the board dimensions. Weights are all 1].\n");
    printf("\n\
-grid_2D [planar bidimensional grid of size <board_height> x <board_width>\n\
         Weights are all 1.]\n");
    printf("\n\
-grid    [<space_dimension>-dimensional cubic grid of size <board_size>\n\
         Weights are all 1.]\n");
    printf("\n\
-toroidal_grid_2D [toroidal bidimensional grid of size\n\
         <board_height> x <board_width>. Weights are all 1.]\n");
    printf("\n\
-toroidal_grid [<space_dimension>-dimensional cubic grid of size\n\
         <board_size> with periodic boundary wrapping. Weights are all 1.]\n");
    printf("\n\
-simplex [generates a simplex graph. A node of the graph is associated to a\n\
         <dimension>-dimensional vector of nonnegative integers whose sum is\n\
         <sum>. Two nodes are adjacent if the Euclidean distance of the\n\
         two associated vectors is sqrt(2). The edge weights are all 1.]\n");
    printf("\n\
-boundend_simplex [it is like -simplex, except that each component of\n\
         the vector has <bound> as an upper limit.]\n");
    printf("\n\
-circuit [generates a 2-regular connected graph of <length> nodes. The\n\
         edge weight are all 1.]\n");
    printf("\n\
-clique  [generates a complete graph of <size> nodes. The weight of the\n\
         edge (i,j) is |i-j|.]\n");
    printf("\n\
-planar  [generates a random planar graph of <size> nodes. The random number\n\
         generator is initialized by <seed>. The parameter <Density> is any\n\
         real in the interval [0,100]. The number of edges of the graph\n\
         is given by 3 * (<size> - 2) * <Density> / 100. Recall that the\n\
         maximum number of edges of a planar graph is 3 * (<size> - 2).\n\
         The weights are all set to 1.]\n");
    printf("\n\
-spinglass2pm [generates a toroidal 2D-grid for a spin glass model with\n\
         +/- J interactions. The grid has size <n_row> x <n_columns>.\n\
         The percentage of negative interactions is <%% of - bonds>.\n\
         This switch is for compatibility with a former generator. \n\
         Warning: the imput graph must be a <simple_graph>!]\n");
    printf("\n\
-spinglass3pm [generates a toroidal 3D-grid for a spin glass model with\n\
         +/-J interactions.\n\
         The grid has size <n_row> x <n_columns> x <n_layers>.\n\
         The percentage of negative interactions is <%% of - bonds>.\n\
         This switch is for compatibility with a former generator. \n\
         Warning: the imput graph must be a <simple_graph>!]\n");
    printf("\n\
-spinglass2g [generates a toroidal 2D grid for a spin glass model with\n\
         gaussian interactions. The grid has size <n_row> x <n_columns>.\n\
         This switch is for compatibility with a former generator. \n\
         Warning: the imput graph must be a <simple_graph>!]\n");
    printf("\n\
-spinglass3g [generates a toroidal 3D grid for a spin glass model with\n\
         +/- J interactions.\n\
         The grid has size <n_row> x <n_columns> x <n_layers>.\n\
         This switch is for compatibility with a former generator. \n\
         Warning: the imput graph must be a <simple_graph>!]\n");
    printf("\n\
-rnd_graph [generates a random graph of <size> nodes and density <Density>.\n\
         The parameter <Density> is any real in the interval [0,100].\n\
         The number of edges of the graph is given by the integer closest\n\
         to <size> * (<size> - 1) * <density> / 200. The random number\n\
         generator is initialized by <seed>. The edge weights are all 1.]\n");
    printf("\n");
    exit(1);
  }
  while (++argp < argc) {
    if (is("+")) {
      {
	if (stack <= 0)
	  goto two_op;
	if (g[stack - 1]->n != g[stack]->n)
	  goto diff_size;
	if (g[stack]->n <= 0)
	  goto empty_graph;
	g[stack - 1] = gunion(g[stack - 1], g[stack], 0L, 0L);
	gb_recycle(g[stack--]);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("x")) {
      {
	if (stack <= 0)
	  goto two_op;
	if (g[stack - 1]->n <= 0 || g[stack]->n <= 0)
	  goto empty_graph;
	g[stack - 1] = product(g[stack], g[stack - 1], cartesian, 0L);
	gb_recycle(g[stack--]);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is(":")) {
      {
	if (stack <= 0)
	  goto two_op;
	if (g[stack - 1]->n <= 0 || g[stack]->n <= 0)
	  goto empty_graph;
	gtmp = board(2L, 0L, 0L, 0L, -1L, 0L, 0L);
	(gtmp->vertices)->ind = IND_GRAPH;
	(gtmp->vertices)->subst = g[stack - 1];
	(gtmp->vertices + 1)->ind = IND_GRAPH;
	(gtmp->vertices + 1)->subst = g[stack];
	gtmp1 = induced(gtmp, NULL, 0L, 0L, 0L);
	gb_recycle(gtmp);
	gb_recycle(g[stack--]);
	gb_recycle(g[stack--]);
	g[++stack] = gtmp1;
	if (g[stack] == NULL)
	  goto panic;
	{
	  register Vertex *v;
	  register Arc   *a;
	  for (v = g[stack]->vertices; v < g[stack]->vertices + g[stack]->n; v++)
	    for (a = v->arcs; a; a = a->next)
	      a->len = 1;
	}
      }
      continue;
    }
    if (is("-complement")) {
      {
	if (stack < 0)
	  goto one_op;
	if (g[stack]->n <= 0)
	  goto empty_graph;
	g[stack] = complement(g[stack], 0L, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-random")) {
      {
	if (argp + 3 >= argc)
	  goto ill_num;
	lower_weight = atoi(argv[++argp]);
	upper_weight = atoi(argv[++argp]);
	seed = atoi(argv[++argp]);
	if (lower_weight < upper_weight)
	  upper_weight++;
	if (lower_weight > upper_weight || seed < 0)
	  goto ill_par;
	if (stack < 0)
	  goto one_op;
	if (g[stack]->n <= 0)
	  goto empty_graph;
	random_lengths(g[stack], 0L, lower_weight, upper_weight, NULL, seed);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-times")) {
      {
	if (argp + 1 >= argc)
	  goto ill_num;
	scalar = atoi(argv[++argp]);
	if (stack < 0)
	  goto one_op;
	if (g[stack]->n <= 0)
	  goto empty_graph;
	{
	  register Vertex *v;
	  register Arc   *a;
	  for (v = g[stack]->vertices; v < g[stack]->vertices + g[stack]->n; v++)
	    for (a = v->arcs; a; a = a->next)
	      a->len *= scalar;
	}
      }
      continue;
    }
    if (is("-plus")) {
      {
	if (argp + 1 >= argc)
	  goto ill_num;
	scalar = atoi(argv[++argp]);
	if (stack < 0)
	  goto one_op;
	if (g[stack]->n <= 0)
	  goto empty_graph;
	{
	  register Vertex *v;
	  register Arc   *a;
	  for (v = g[stack]->vertices; v < g[stack]->vertices + g[stack]->n; v++)
	    for (a = v->arcs; a; a = a->next)
	      a->len += scalar;
	}
      }
      continue;
    }
    if (is("-line")) {
      {
	if (stack < 0)
	  goto one_op;
	if (g[stack]->n <= 0)
	  goto empty_graph;
	g[stack] = lines(g[stack], 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-circuit")) {
      {
	if (argp + 1 >= argc)
	  goto ill_num;
	length = atoi(argv[++argp]);
	if (length <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = circuit(length);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-clique")) {
      {
	if (argp + 1 >= argc)
	  goto ill_num;
	size = atoi(argv[++argp]);
	if (size <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = complete(size);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-planar")) {
      {
	if (argp + 3 >= argc)
	  goto ill_num;
	size = atoi(argv[++argp]);
	sscanf(argv[++argp], "%lf", &density);
	seed = atoi(argv[++argp]);
	if (size <= 0 || density < 0.0 || density > 100.0 || seed < 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = planar(size, density, seed);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-leap_2D")) {
      {
	if (argp + 3 >= argc)
	  goto ill_num;
	board_height = atoi(argv[++argp]);
	board_width = atoi(argv[++argp]);
	move_type = atoi(argv[++argp]);
	if (board_height <= 0 || board_width <= 0 || move_type <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_height, board_width, 0L, 0L, -move_type, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-leap")) {
      {
	if (argp + 3 >= argc)
	  goto ill_num;
	board_size = atoi(argv[++argp]);
	space_dimension = atoi(argv[++argp]);
	move_type = atoi(argv[++argp]);
	if (space_dimension <= 0 || board_size <= 0 || move_type <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_size, -space_dimension, 0L, 0L, -move_type, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-wrapped_leap_2D")) {
      {
	if (argp + 3 >= argc)
	  goto ill_num;
	board_height = atoi(argv[++argp]);
	board_width = atoi(argv[++argp]);
	move_type = atoi(argv[++argp]);
	if (board_height <= 0 || board_width <= 0 || move_type <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_height, board_width, 0L, 0L, -move_type, -1L, 0L);
	if (g[stack] == NULL)
	  goto panic;
	g[stack] = complement(g[stack], 1L, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-wrapped_leap")) {
      {
	if (argp + 3 >= argc)
	  goto ill_num;
	board_size = atoi(argv[++argp]);
	space_dimension = atoi(argv[++argp]);
	move_type = atoi(argv[++argp]);
	if (space_dimension <= 0 || board_size <= 0 || move_type <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_size, -space_dimension, 0L, 0L, -move_type, -1L, 0L);
	if (g[stack] == NULL)
	  goto panic;
	g[stack] = complement(g[stack], 1L, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-grid_2D")) {
      {
	if (argp + 2 >= argc)
	  goto ill_num;
	board_height = atoi(argv[++argp]);
	board_width = atoi(argv[++argp]);
	if (board_height <= 0 || board_width <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_height, board_width, 0L, 0L, 1L, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-grid")) {
      {
	if (argp + 2 >= argc)
	  goto ill_num;
	board_size = atoi(argv[++argp]);
	space_dimension = atoi(argv[++argp]);
	if (space_dimension <= 0 || board_size <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_size, -space_dimension, 0L, 0L, 1L, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-toroidal_grid_2D")) {
      {
	if (argp + 2 >= argc)
	  goto ill_num;
	board_height = atoi(argv[++argp]);
	board_width = atoi(argv[++argp]);
	if (board_height <= 0 || board_width <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_height, board_width, 0L, 0L, 1L, -1L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-toroidal_grid")) {
      {
	if (argp + 2 >= argc)
	  goto ill_num;
	board_size = atoi(argv[++argp]);
	space_dimension = atoi(argv[++argp]);
	if (space_dimension <= 0 || board_size <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = board(board_size, -space_dimension, 0L, 0L, 1L, -1L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-simplex")) {
      {
	if (argp + 2 >= argc)
	  goto ill_num;
	sum = atoi(argv[++argp]);
	dimension = atoi(argv[++argp]);
	if (sum <= 0 || dimension <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = simplex(sum, -dimension, 0L, 0L, 0L, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-bounded_simplex")) {
      {
	if (argp + 3 >= argc)
	  goto ill_num;
	sum = atoi(argv[++argp]);
	dimension = atoi(argv[++argp]);
	bound = atoi(argv[++argp]);
	if (sum <= 0 || dimension <= 0 || bound <= 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	g[++stack] = simplex(sum, bound, -dimension, 0L, 0L, 0L, 0L);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-rnd_graph")) {
      {
	long            n_edges;
	double          rn_edges;
	if (argp + 3 >= argc)
	  goto ill_num;
	size = atoi(argv[++argp]);
	sscanf(argv[++argp], "%lf", &density);
	seed = atoi(argv[++argp]);
	if (size <= 0 || density < 0.0 || density > 100.0 || seed < 0)
	  goto ill_par;
	if (stack + 1 >= max_stack)
	  goto too_many;
	rn_edges = (size * (size - 1) * density) / 200.0;
	n_edges = rn_edges;
	if (rn_edges - n_edges >= 0.5)
	  n_edges++;
	g[++stack] = random_graph(size, n_edges, 0, 0, 0, NULL, NULL, 1, 1, seed);
	if (g[stack] == NULL)
	  goto panic;
      }
      continue;
    }
    if (is("-spinglass2pm")) {
      {
	int             r, c, perc, seed, nnodes, nedges;
	double         *cost;
	if (stack != -1)
	  goto simp_gr;
	if (argp + 4 >= argc)
	  goto ill_num;
	if ((argp + 5) != argc)
	  goto simp_gr;
	r = atoi(argv[++argp]);
	c = atoi(argv[++argp]);
	perc = atoi(argv[++argp]);
	seed = atoi(argv[++argp]);
	nnodes = r * c;
	nedges = 2 * nnodes;
	cost = (double *) malloc((nedges + 1) * sizeof(double));
	{
	  int             j, k, t, nummo;
	  gb_init_rand(seed);
	  nummo = (nedges * perc) / 100;
	  for (j = 1; j <= nummo; j++)
	    cost[j] = -1.0;
	  for (j = nummo + 1; j <= nedges; j++)
	    cost[j] = 1.0;
	  for (j = nedges; j > 1; j--) {
	    k = gb_unif_rand(j) + 1;
	    t = cost[j];
	    cost[j] = cost[k];
	    cost[k] = t;
	  }
	}
	{
	  int             k, i, j;
	  printf("%d %d \n", nnodes, nedges);
	  for (k = 1; k <= nedges; k++) {
	    if (k <= nnodes) {
	      i = k;
	      j = k + 1;
	      if (!(i % c))
		j -= c;
	    } else {
	      i = k - nnodes;
	      j = i + c;
	      if (j > nnodes)
		j = j - nnodes;
	    }
	    printf("%d %d %1.0f\n", i, j, cost[k]);
	  }
	}
	free((char *) cost);
	exit(0);
      }
      continue;
    }
    if (is("-spinglass3pm")) {
      {
	int             r, c, l, perc, seed, nnodes, nedges;
	double         *cost;
	if (stack != -1)
	  goto simp_gr;
	if (argp + 5 >= argc)
	  goto ill_num;
	if ((argp + 6) != argc)
	  goto simp_gr;
	r = atoi(argv[++argp]);
	c = atoi(argv[++argp]);
	l = atoi(argv[++argp]);
	perc = atoi(argv[++argp]);
	seed = atoi(argv[++argp]);
	nnodes = r * c * l;
	nedges = 3 * nnodes;
	cost = (double *) malloc((nedges + 1) * sizeof(double));
	{
	  int             j, k, t, nummo;
	  gb_init_rand(seed);
	  nummo = (nedges * perc) / 100;
	  for (j = 1; j <= nummo; j++)
	    cost[j] = -1.0;
	  for (j = nummo + 1; j <= nedges; j++)
	    cost[j] = 1.0;
	  for (j = nedges; j > 1; j--) {
	    k = gb_unif_rand(j) + 1;
	    t = cost[j];
	    cost[j] = cost[k];
	    cost[k] = t;
	  }
	}
	{
	  int             i, j, rc, layer, rcl1, a, b, e;
	  printf("%d %d \n", nnodes, nedges);
	  rc = r * c;
	  for (e = 1; e <= nedges; e++) {
	    layer = (e - 1) / (3 * rc) + 1;
	    rcl1 = rc * (layer - 1);
	    a = 3 * rcl1 + rc;
	    b = a + rc;
	    if (e <= a) {
	      i = e - 2 * rcl1;
	      j = i + 1;
	      if (!(i % c))
		j -= c;
	    } else if (e <= b) {
	      i = e - (2 * rcl1) - rc;
	      j = i + c;
	      if (j > layer * rc)
		j -= rc;
	    } else {
	      i = e - 2 * rcl1 - (2 * rc);
	      j = i + rc;
	      if (j > nnodes)
		j -= nnodes;
	    }
	    printf("%d %d %1.0f\n", i, j, cost[e]);
	  }
	}
	free((char *) cost);
	exit(0);
      }
      continue;
    }
    if (is("-spinglass2g")) {
      {
	int             r, c, seed, nnodes, nedges;
	double         *cost;
	if (stack != -1)
	  goto simp_gr;
	if (argp + 3 >= argc)
	  goto ill_num;
	if ((argp + 4) != argc)
	  goto simp_gr;
	r = atoi(argv[++argp]);
	c = atoi(argv[++argp]);
	seed = atoi(argv[++argp]);
	nnodes = r * c;
	nedges = 2 * nnodes;
	cost = (double *) malloc((nedges + 1) * sizeof(double));
	{
	  int             j, k, t;
	  int             scalefactor = 100000;
	  gb_init_rand(seed);
	  for (j = 1; j <= nedges; j++)
	    cost[j] = gauss() * scalefactor;
	  for (j = nedges; j > 1; j--) {
	    k = gb_unif_rand(j) + 1;
	    t = cost[j];
	    cost[j] = cost[k];
	    cost[k] = t;
	  }
	}
	{
	  int             k, i, j;
	  printf("%d %d \n", nnodes, nedges);
	  for (k = 1; k <= nedges; k++) {
	    if (k <= nnodes) {
	      i = k;
	      j = k + 1;
	      if (!(i % c))
		j -= c;
	    } else {
	      i = k - nnodes;
	      j = i + c;
	      if (j > nnodes)
		j = j - nnodes;
	    }
	    printf("%d %d %1.0f\n", i, j, cost[k]);
	  }
	}
	free((char *) cost);
	exit(0);
      }
      continue;
    }
    if (is("-spinglass3g")) {
      {
	int             r, c, l, seed, nnodes, nedges;
	double         *cost;
	if (stack != -1)
	  goto simp_gr;
	if (argp + 4 >= argc)
	  goto ill_num;
	if ((argp + 5) != argc)
	  goto simp_gr;
	r = atoi(argv[++argp]);
	c = atoi(argv[++argp]);
	l = atoi(argv[++argp]);
	seed = atoi(argv[++argp]);
	nnodes = r * c * l;
	nedges = 3 * nnodes;
	cost = (double *) malloc((nedges + 1) * sizeof(double));
	{
	  int             j, k, t;
	  int             scalefactor = 100000;
	  gb_init_rand(seed);
	  for (j = 1; j <= nedges; j++)
	    cost[j] = gauss() * scalefactor;
	  for (j = nedges; j > 1; j--) {
	    k = gb_unif_rand(j) + 1;
	    t = cost[j];
	    cost[j] = cost[k];
	    cost[k] = t;
	  }
	}
	{
	  int             i, j, rc, layer, rcl1, a, b, e;
	  printf("%d %d \n", nnodes, nedges);
	  rc = r * c;
	  for (e = 1; e <= nedges; e++) {
	    layer = (e - 1) / (3 * rc) + 1;
	    rcl1 = rc * (layer - 1);
	    a = 3 * rcl1 + rc;
	    b = a + rc;
	    if (e <= a) {
	      i = e - 2 * rcl1;
	      j = i + 1;
	      if (!(i % c))
		j -= c;
	    } else if (e <= b) {
	      i = e - (2 * rcl1) - rc;
	      j = i + c;
	      if (j > layer * rc)
		j -= rc;
	    } else {
	      i = e - 2 * rcl1 - (2 * rc);
	      j = i + rc;
	      if (j > nnodes)
		j -= nnodes;
	    }
	    printf("%d %d %1.0f\n", i, j, cost[e]);
	  }
	}
	free((char *) cost);
	exit(0);
      }
      continue;
    }
    goto error;
  }
  {
    register Vertex *v;
    register Arc   *a;
    long            i, j, free_util;
    if (stack != 0)
      goto not_empty;
    for (free_util = 0; free_util < 6; free_util++)
      if (g[0]->util_types[free_util] == 'Z')
	break;
    if (free_util >= 6)
      goto no_utility;
    i = 0;
    for (v = g[0]->vertices; v < g[0]->vertices + g[0]->n; v++) {
      switch (free_util) {
      case 0:
	v->u.I = ++i;
	break;
      case 1:
	v->v.I = ++i;
	break;
      case 2:
	v->w.I = ++i;
	break;
      case 3:
	v->x.I = ++i;
	break;
      case 4:
	v->y.I = ++i;
	break;
      case 5:
	v->z.I = ++i;
	break;
      }
    }
    printf("%ld %ld \n", g[0]->n, g[0]->m / 2);
    for (v = g[0]->vertices; v < g[0]->vertices + g[0]->n; v++)
      for (a = v->arcs; a; a = a->next) {
	switch (free_util) {
	case 0:
	  i = v->u.I;
	  j = a->tip->u.I;
	  break;
	case 1:
	  i = v->v.I;
	  j = a->tip->v.I;
	  break;
	case 2:
	  i = v->w.I;
	  j = a->tip->w.I;
	  break;
	case 3:
	  i = v->x.I;
	  j = a->tip->x.I;
	  break;
	case 4:
	  i = v->y.I;
	  j = a->tip->y.I;
	  break;
	case 5:
	  i = v->z.I;
	  j = a->tip->z.I;
	  break;
	}
	if (i < j)
	  printf("%ld %ld %ld\n", i, j, a->len);
      }
#ifdef CONN_COMP
    printf("%ld\n", g[0]->n);
    for (v = g[0]->vertices; v < g[0]->vertices + g[0]->n; v++) {
      switch (free_util) {
      case 0:
	i = v->u.I;
	break;
      case 1:
	i = v->v.I;
	break;
      case 2:
	i = v->w.I;
	break;
      case 3:
	i = v->x.I;
	break;
      case 4:
	i = v->y.I;
	break;
      case 5:
	i = v->z.I;
	break;
      }
      printf("%ld\n", i);
    }
#endif
  }
  exit(0);
error:
  fprintf(stderr, "Unknown rudy option\n");
  goto terminate;
empty_graph:
  fprintf(stderr, "An operand is empty\n");
  goto terminate;
one_op:
  fprintf(stderr, "The operand requires one argument\n");
  goto terminate;
two_op:
  fprintf(stderr, "The operand requires two arguments\n");
  goto terminate;
too_many:
  fprintf(stderr, "Too many operands in the stack\n");
  goto terminate;
diff_size:
  fprintf(stderr, "The two arguments have different sizes\n");
  goto terminate;
ill_par:
  fprintf(stderr, "Illegal parameters\n");
  goto terminate;
ill_num:
  fprintf(stderr, "Illegal number of parameters\n");
  goto terminate;
simp_gr:
  fprintf(stderr, "The input must be a <simple_graph>\n");
  goto terminate;
not_empty:
  fprintf(stderr, "The stack is not empty\n");
  goto terminate;
no_utility:
  fprintf(stderr, "No utility fields available\n");
  goto terminate;
panic:
  fprintf(stderr, "Something went wrong (panic code %ld)!\n", panic_code);
  goto terminate;
terminate:
  {
    long            i;
    for (i = 0; i <= argp; i++)
      fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, "\n");
    exit(1);
  }
}



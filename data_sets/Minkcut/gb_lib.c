/*  LAST EDIT: Fri Sep 29 17:58:46 1995 by Christoph Helmberg (kombi8!bzfhelmb)  */
/* -------------------------- gb_lib.c -------------------- */
#include <stdio.h>
#ifdef SYSV
#include <string.h>
#else
#include <string.h>
#endif
#undef min

#define gb_next_rand() (*gb_fptr>=0?*gb_fptr--:gb_flip_cycle())
extern long *gb_fptr;
extern long gb_flip_cycle ();

extern void gb_init_rand ();

extern long gb_unif_rand ();

typedef union {
  struct vertex_struct *V;
  struct arc_struct *A;
  struct graph_struct *G;
  char *S;
  long I;
}

util;

typedef struct vertex_struct {
  struct arc_struct *arcs;
  char *name;
  util u, v, w, x, y, z;
}

Vertex;

typedef struct arc_struct {
  struct vertex_struct *tip;
  struct arc_struct *next;
  long len;
  util a, b;
}

Arc;

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
}

Graph;

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

#define panic(c) \
{panic_code= c; \
gb_free(working_storage); \
gb_trouble_code= 0; \
return NULL; \
} \

#define BUF_SIZE 4096 \

#define MAX_D 91 \

#define MAX_NNN 1000000000.0 \

#define UL_BITS 8*sizeof(unsigned long) \

#define vert_offset(v,delta)((Vertex*)(((siz_t)v)+delta)) \
 \

#define tmp u.V \

#define tlen z.A \

#define mult v.I
#define minlen w.I \

#define map z.V \

#define ind z.I \

#define IND_GRAPH 1000000000
#define subst y.G \

static Area working_storage;

static char buffer[BUF_SIZE];

static long nn[MAX_D + 1];
static long wr[MAX_D + 1];
static long del[MAX_D + 1];
static long sig[MAX_D + 2];
static long xx[MAX_D + 1], yy[MAX_D + 1];

static char *short_imap = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ\
abcdefghijklmnopqrstuvwxyz_^~&@,;.:?!%#$+-*/|<=>()[]{}`'";

Graph *
board (n1, n2, n3, n4, piece, wrap, directed)
     long n1, n2, n3, n4;
     long piece;
     long wrap;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  long n;
  long p;
  long l;

  if (piece == 0)
    piece = 1;
  if (n1 <= 0) {
    n1 = n2 = 8;
    n3 = 0;
  }
  nn[1] = n1;
  if (n2 <= 0) {
    k = 2;
    d = -n2;
    n3 = n4 = 0;
  }
  else {
    nn[2] = n2;
    if (n3 <= 0) {
      k = 3;
      d = -n3;
      n4 = 0;
    }
    else {
      nn[3] = n3;
      if (n4 <= 0) {
	k = 4;
	d = -n4;
      }
      else {
	nn[4] = n4;
	d = 4;
	goto done;
      }
    }
  }
  if (d == 0) {
    d = k - 1;
    goto done;
  }

  if (d > MAX_D)
    panic (bad_specs);
  for (j = 1; k <= d; j++, k++)
    nn[k] = nn[j];

done:

  {
    float nnn;
    for (n = 1, nnn = 1.0, j = 1; j <= d; j++) {
      nnn *= (float) nn[j];
      if (nnn > MAX_NNN)
	panic (very_bad_specs);
      n *= nn[j];
    }
    new_graph = gb_new_graph (n);
    if (new_graph == NULL)
      panic (no_room);
    sprintf (new_graph->id, "board(%ld,%ld,%ld,%ld,%ld,%ld,%d)",
	     n1, n2, n3, n4, piece, wrap, directed ? 1 : 0);
    strcpy (new_graph->util_types, "ZZZIIIZZZZZZZZ");

    {
      register char *q;
      nn[0] = xx[0] = xx[1] = xx[2] = xx[3] = 0;
      for (k = 4; k <= d; k++)
	xx[k] = 0;
      for (v = new_graph->vertices;; v++) {
	q = buffer;
	for (k = 1; k <= d; k++) {
	  sprintf (q, ".%ld", xx[k]);
	  while (*q)
	    q++;
	}
	v->name = gb_save_string (&buffer[1]);
	v->x.I = xx[1];
	v->y.I = xx[2];
	v->z.I = xx[3];
	for (k = d; xx[k] + 1 == nn[k]; k--)
	  xx[k] = 0;
	if (k == 0)
	  break;
	xx[k]++;
      }
    }

  }

  {
    register long w = wrap;
    for (k = 1; k <= d; k++, w >>= 1) {
      wr[k] = w & 1;
      del[k] = sig[k] = 0;
    }
    sig[0] = del[0] = sig[d + 1] = 0;
  }

  p = piece;
  if (p < 0)
    p = -p;
  while (1) {

    for (k = d; sig[k] + (del[k] + 1) * (del[k] + 1) > p; k--)
      del[k] = 0;
    if (k == 0)
      break;
    del[k]++;
    sig[k + 1] = sig[k] + del[k] * del[k];
    for (k++; k <= d; k++)
      sig[k + 1] = sig[k];
    if (sig[d + 1] < p)
      continue;

    while (1) {

      for (k = 1; k <= d; k++)
	xx[k] = 0;
      for (v = new_graph->vertices;; v++) {

	for (k = 1; k <= d; k++)
	  yy[k] = xx[k] + del[k];
	for (l = 1;; l++) {

	  for (k = 1; k <= d; k++) {
	    if (yy[k] < 0) {
	      if (!wr[k])
		goto no_more;
	      do
		yy[k] += nn[k];
	      while (yy[k] < 0);
	    }
	    else if (yy[k] >= nn[k]) {
	      if (!wr[k])
		goto no_more;
	      do
		yy[k] -= nn[k];
	      while (yy[k] >= nn[k]);
	    }
	  }

	  if (piece < 0) {
	    for (k = 1; k <= d; k++)
	      if (yy[k] != xx[k])
		goto unequal;
	    goto no_more;
	  unequal:;
	  }

	  for (k = 2, j = yy[1]; k <= d; k++)
	    j = nn[k] * j + yy[k];
	  if (directed)
	    gb_new_arc (v, new_graph->vertices + j, l);
	  else
	    gb_new_edge (v, new_graph->vertices + j, l);

	  if (piece > 0)
	    goto no_more;
	  for (k = 1; k <= d; k++)
	    yy[k] += del[k];
	}
      no_more:

	;
	for (k = d; xx[k] + 1 == nn[k]; k--)
	  xx[k] = 0;
	if (k == 0)
	  break;
	xx[k]++;
      }

      for (k = d; del[k] <= 0; k--)
	del[k] = -del[k];
      if (sig[k] == 0)
	break;
      del[k] = -del[k];

    }
  }

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
simplex (n, n0, n1, n2, n3, n4, directed)
     unsigned long n;
     long n0, n1, n2, n3, n4;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  if (n0 == 0)
    n0 = -2;
  if (n0 < 0) {
    k = 2;
    nn[0] = n;
    d = -n0;
    n1 = n2 = n3 = n4 = 0;
  }
  else {
    if (n0 > n)
      n0 = n;
    nn[0] = n0;
    if (n1 <= 0) {
      k = 2;
      d = -n1;
      n2 = n3 = n4 = 0;
    }
    else {
      if (n1 > n)
	n1 = n;
      nn[1] = n1;
      if (n2 <= 0) {
	k = 3;
	d = -n2;
	n3 = n4 = 0;
      }
      else {
	if (n2 > n)
	  n2 = n;
	nn[2] = n2;
	if (n3 <= 0) {
	  k = 4;
	  d = -n3;
	  n4 = 0;
	}
	else {
	  if (n3 > n)
	    n3 = n;
	  nn[3] = n3;
	  if (n4 <= 0) {
	    k = 5;
	    d = -n4;
	  }
	  else {
	    if (n4 > n)
	      n4 = n;
	    nn[4] = n4;
	    d = 4;
	    goto done;
	  }
	}
      }
    }
  }
  if (d == 0) {
    d = k - 2;
    goto done;
  }
  nn[k - 1] = nn[0];

  if (d > MAX_D)
    panic (bad_specs);
  for (j = 1; k <= d; j++, k++)
    nn[k] = nn[j];

done:

  {
    long nverts;
    register long *coef = gb_typed_alloc (n + 1, long, working_storage);
    if (gb_trouble_code)
      panic (no_room + 1);
    for (k = 0; k <= nn[0]; k++)
      coef[k] = 1;

    for (j = 1; j <= d; j++) {
      for (k = n, i = n - nn[j] - 1; i >= 0; k--, i--)
	coef[k] -= coef[i];
      s = 1;
      for (k = 1; k <= n; k++) {
	s += coef[k];
	if (s > 1000000000)
	  panic (very_bad_specs);
	coef[k] = s;
      }
    }

    nverts = coef[n];
    gb_free (working_storage);
    new_graph = gb_new_graph (nverts);
    if (new_graph == NULL)
      panic (no_room);
  }

  sprintf (new_graph->id, "simplex(%lu,%ld,%ld,%ld,%ld,%ld,%d)",
	   n, n0, n1, n2, n3, n4, directed ? 1 : 0);
  strcpy (new_graph->util_types, "VVZIIIZZZZZZZZ");

  v = new_graph->vertices;
  yy[d + 1] = 0;
  sig[0] = n;
  for (k = d; k >= 0; k--)
    yy[k] = yy[k + 1] + nn[k];
  if (yy[0] >= n) {
    k = 0;
    xx[0] = (yy[1] >= n ? 0 : n - yy[1]);
    while (1) {

      for (s = sig[k] - xx[k], k++; k <= d; s -= xx[k], k++) {
	sig[k] = s;
	if (s <= yy[k + 1])
	  xx[k] = 0;
	else
	  xx[k] = s - yy[k + 1];
      }
      if (s != 0)
	panic (impossible + 1);

      {
	register char *p = buffer;
	for (k = 0; k <= d; k++) {
	  sprintf (p, ".%ld", xx[k]);
	  while (*p)
	    p++;
	}
	v->name = gb_save_string (&buffer[1]);
	v->x.I = xx[0];
	v->y.I = xx[1];
	v->z.I = xx[2];
      }

      hash_in (v);

      for (j = 0; j < d; j++)
	if (xx[j]) {
	  register Vertex *u;
	  xx[j]--;
	  for (k = j + 1; k <= d; k++)
	    if (xx[k] < nn[k]) {
	      register char *p = buffer;
	      xx[k]++;
	      for (i = 0; i <= d; i++) {
		sprintf (p, ".%ld", xx[i]);
		while (*p)
		  p++;
	      }
	      u = hash_out (&buffer[1]);
	      if (u == NULL)
		panic (impossible + 2);
	      if (directed)
		gb_new_arc (u, v, 1L);
	      else
		gb_new_edge (u, v, 1L);
	      xx[k]--;
	    }
	  xx[j]++;
	}

      v++;

      for (k = d - 1;; k--) {
	if (xx[k] < sig[k] && xx[k] < nn[k])
	  break;
	if (k == 0)
	  goto last;
      }
      xx[k]++;

    }
  }
last:if (v != new_graph->vertices + new_graph->n)
    panic (impossible);

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
subsets (n, n0, n1, n2, n3, n4, size_bits, directed)
     unsigned long n;
     long n0, n1, n2, n3, n4;
     unsigned long size_bits;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  if (n0 == 0)
    n0 = -2;
  if (n0 < 0) {
    k = 2;
    nn[0] = n;
    d = -n0;
    n1 = n2 = n3 = n4 = 0;
  }
  else {
    if (n0 > n)
      n0 = n;
    nn[0] = n0;
    if (n1 <= 0) {
      k = 2;
      d = -n1;
      n2 = n3 = n4 = 0;
    }
    else {
      if (n1 > n)
	n1 = n;
      nn[1] = n1;
      if (n2 <= 0) {
	k = 3;
	d = -n2;
	n3 = n4 = 0;
      }
      else {
	if (n2 > n)
	  n2 = n;
	nn[2] = n2;
	if (n3 <= 0) {
	  k = 4;
	  d = -n3;
	  n4 = 0;
	}
	else {
	  if (n3 > n)
	    n3 = n;
	  nn[3] = n3;
	  if (n4 <= 0) {
	    k = 5;
	    d = -n4;
	  }
	  else {
	    if (n4 > n)
	      n4 = n;
	    nn[4] = n4;
	    d = 4;
	    goto done;
	  }
	}
      }
    }
  }
  if (d == 0) {
    d = k - 2;
    goto done;
  }
  nn[k - 1] = nn[0];

  if (d > MAX_D)
    panic (bad_specs);
  for (j = 1; k <= d; j++, k++)
    nn[k] = nn[j];

done:

  {
    long nverts;
    register long *coef = gb_typed_alloc (n + 1, long, working_storage);
    if (gb_trouble_code)
      panic (no_room + 1);
    for (k = 0; k <= nn[0]; k++)
      coef[k] = 1;

    for (j = 1; j <= d; j++) {
      for (k = n, i = n - nn[j] - 1; i >= 0; k--, i--)
	coef[k] -= coef[i];
      s = 1;
      for (k = 1; k <= n; k++) {
	s += coef[k];
	if (s > 1000000000)
	  panic (very_bad_specs);
	coef[k] = s;
      }
    }

    nverts = coef[n];
    gb_free (working_storage);
    new_graph = gb_new_graph (nverts);
    if (new_graph == NULL)
      panic (no_room);
  }

  sprintf (new_graph->id, "subsets(%lu,%ld,%ld,%ld,%ld,%ld,0x%lx,%d)",
	   n, n0, n1, n2, n3, n4, size_bits, directed ? 1 : 0);
  strcpy (new_graph->util_types, "ZZZIIIZZZZZZZZ");

  v = new_graph->vertices;
  yy[d + 1] = 0;
  sig[0] = n;
  for (k = d; k >= 0; k--)
    yy[k] = yy[k + 1] + nn[k];
  if (yy[0] >= n) {
    k = 0;
    xx[0] = (yy[1] >= n ? 0 : n - yy[1]);
    while (1) {

      for (s = sig[k] - xx[k], k++; k <= d; s -= xx[k], k++) {
	sig[k] = s;
	if (s <= yy[k + 1])
	  xx[k] = 0;
	else
	  xx[k] = s - yy[k + 1];
      }
      if (s != 0)
	panic (impossible + 1);

      {
	register char *p = buffer;
	for (k = 0; k <= d; k++) {
	  sprintf (p, ".%ld", xx[k]);
	  while (*p)
	    p++;
	}
	v->name = gb_save_string (&buffer[1]);
	v->x.I = xx[0];
	v->y.I = xx[1];
	v->z.I = xx[2];
      }

      {
	register Vertex *u;
	for (u = new_graph->vertices; u <= v; u++) {
	  register char *p = u->name;
	  long ss = 0;
	  for (j = 0; j <= d; j++, p++) {
	    for (s = (*p++) - '0'; *p >= '0'; p++)
	      s = 10 * s + *p - '0';

	    if (xx[j] < s)
	      ss += xx[j];
	    else
	      ss += s;
	  }
	  if ((size_bits & (((unsigned long) 1) << ss)) && ss < UL_BITS) {
	    if (directed)
	      gb_new_arc (u, v, 1L);
	    else
	      gb_new_edge (u, v, 1L);
	  }
	}
      }

      v++;

      for (k = d - 1;; k--) {
	if (xx[k] < sig[k] && xx[k] < nn[k])
	  break;
	if (k == 0)
	  goto last;
      }
      xx[k]++;

    }
  }
last:if (v != new_graph->vertices + new_graph->n)
    panic (impossible);

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
perms (n0, n1, n2, n3, n4, max_inv, directed)
     long n0, n1, n2, n3, n4;
     unsigned long max_inv;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  register long n;

  if (n0 == 0) {
    n0 = 1;
    n1 = 0;
  }
  else if (n0 < 0) {
    n1 = n0;
    n0 = 1;
  }
  n = BUF_SIZE;

  if (n0 == 0)
    n0 = -2;
  if (n0 < 0) {
    k = 2;
    nn[0] = n;
    d = -n0;
    n1 = n2 = n3 = n4 = 0;
  }
  else {
    if (n0 > n)
      n0 = n;
    nn[0] = n0;
    if (n1 <= 0) {
      k = 2;
      d = -n1;
      n2 = n3 = n4 = 0;
    }
    else {
      if (n1 > n)
	n1 = n;
      nn[1] = n1;
      if (n2 <= 0) {
	k = 3;
	d = -n2;
	n3 = n4 = 0;
      }
      else {
	if (n2 > n)
	  n2 = n;
	nn[2] = n2;
	if (n3 <= 0) {
	  k = 4;
	  d = -n3;
	  n4 = 0;
	}
	else {
	  if (n3 > n)
	    n3 = n;
	  nn[3] = n3;
	  if (n4 <= 0) {
	    k = 5;
	    d = -n4;
	  }
	  else {
	    if (n4 > n)
	      n4 = n;
	    nn[4] = n4;
	    d = 4;
	    goto done;
	  }
	}
      }
    }
  }
  if (d == 0) {
    d = k - 2;
    goto done;
  }
  nn[k - 1] = nn[0];

  if (d > MAX_D)
    panic (bad_specs);
  for (j = 1; k <= d; j++, k++)
    nn[k] = nn[j];

done:

  {
    register long ss;
    for (k = 0, s = ss = 0; k <= d; ss += s * nn[k], s += nn[k], k++)
      if (nn[k] >= BUF_SIZE)
	panic (bad_specs);

    if (s >= BUF_SIZE)
      panic (bad_specs + 1);
    n = s;
    if (max_inv == 0 || max_inv > ss)
      max_inv = ss;
  }

  {
    long nverts;
    register long *coef = gb_typed_alloc (max_inv + 1, long, working_storage);
    if (gb_trouble_code)
      panic (no_room + 1);
    coef[0] = 1;
    for (j = 1, s = nn[0]; j <= d; s += nn[j], j++)
      for (k = 1; k <= nn[j]; k++) {
	register long ii;
	for (i = max_inv, ii = i - k - s; ii >= 0; ii--, i--)
	  coef[i] -= coef[ii];
	for (i = k, ii = 0; i <= max_inv; i++, ii++) {
	  coef[i] += coef[ii];
	  if (coef[i] > 1000000000)
	    panic (very_bad_specs + 1);
	}
      }

    for (k = 1, nverts = 1; k <= max_inv; k++) {
      nverts += coef[k];
      if (nverts > 1000000000)
	panic (very_bad_specs);
    }
    gb_free (working_storage);
    new_graph = gb_new_graph (nverts);
    if (new_graph == NULL)
      panic (no_room);
    sprintf (new_graph->id, "perms(%ld,%ld,%ld,%ld,%ld,%lu,%d)",
	     n0, n1, n2, n3, n4, max_inv, directed ? 1 : 0);
    strcpy (new_graph->util_types, "VVZZZZZZZZZZZZ");
  }

  {
    register long *xtab, *ytab, *ztab;
    long m = 0;

    xtab = gb_typed_alloc (3 * n + 3, long, working_storage);
    if (gb_trouble_code) {
      gb_recycle (new_graph);
      panic (no_room + 2);
    }
    ytab = xtab + (n + 1);
    ztab = ytab + (n + 1);
    for (j = 0, k = 1, s = nn[0];; k++) {
      xtab[k] = ztab[k] = j;
      if (k == s) {
	if (++j > d)
	  break;
	else
	  s += nn[j];
      }
    }

    v = new_graph->vertices;
    while (1) {

      {
	register char *p;
	register long *q;
	for (p = &buffer[n - 1], q = &xtab[n]; q > xtab; p--, q--)
	  *p = short_imap[*q];
	v->name = gb_save_string (buffer);
	hash_in (v);

      }

      for (j = 1; j < n; j++)
	if (xtab[j] > xtab[j + 1]) {
	  register Vertex *u;

	  buffer[j - 1] = short_imap[xtab[j + 1]];
	  buffer[j] = short_imap[xtab[j]];
	  u = hash_out (buffer);
	  if (u == NULL)
	    panic (impossible + 2);
	  if (directed)
	    gb_new_arc (u, v, 1L);
	  else
	    gb_new_edge (u, v, 1L);
	  buffer[j - 1] = short_imap[xtab[j]];
	  buffer[j] = short_imap[xtab[j + 1]];
	}

      v++;

      for (k = n; k; k--) {
	if (m < max_inv && ytab[k] < k - 1)
	  if (ytab[k] < ytab[k - 1] || ztab[k] > ztab[k - 1])
	    goto move;
	if (ytab[k]) {
	  for (j = k - ytab[k]; j < k; j++)
	    xtab[j] = xtab[j + 1];
	  m -= ytab[k];
	  ytab[k] = 0;
	  xtab[k] = ztab[k];
	}
      }
      goto last;
    move:j = k - ytab[k];
      xtab[j] = xtab[j - 1];
      xtab[j - 1] = ztab[k];
      ytab[k]++;
      m++;

    }
  last:if (v != new_graph->vertices + new_graph->n)
      panic (impossible);
    gb_free (working_storage);
  }

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
parts (n, max_parts, max_size, directed)
     unsigned long n;
     unsigned long max_parts;
     unsigned long max_size;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  if (max_parts == 0 || max_parts > n)
    max_parts = n;
  if (max_size == 0 || max_size > n)
    max_size = n;
  if (max_parts > MAX_D)
    panic (bad_specs);

  {
    long nverts;
    register long *coef = gb_typed_alloc (n + 1, long, working_storage);
    if (gb_trouble_code)
      panic (no_room + 1);
    coef[0] = 1;
    for (k = 1; k <= max_parts; k++) {
      for (j = n, i = n - k - max_size; i >= 0; i--, j--)
	coef[j] -= coef[i];
      for (j = k, i = 0; j <= n; i++, j++) {
	coef[j] += coef[i];
	if (coef[j] > 1000000000)
	  panic (very_bad_specs);
      }
    }
    nverts = coef[n];
    gb_free (working_storage);
    new_graph = gb_new_graph (nverts);
    if (new_graph == NULL)
      panic (no_room);
    sprintf (new_graph->id, "parts(%lu,%lu,%lu,%d)",
	     n, max_parts, max_size, directed ? 1 : 0);
    strcpy (new_graph->util_types, "VVZZZZZZZZZZZZ");
  }

  v = new_graph->vertices;
  xx[0] = max_size;
  sig[1] = n;
  for (k = max_parts, s = 1; k > 0; k--, s++)
    yy[k] = s;
  if (max_size * max_parts >= n) {
    k = 1;
    xx[1] = (n - 1) / max_parts + 1;
    while (1) {

      for (s = sig[k] - xx[k], k++; s; k++) {
	sig[k] = s;
	xx[k] = (s - 1) / yy[k] + 1;
	s -= xx[k];
      }
      d = k - 1;

      {
	register char *p = buffer;
	for (k = 1; k <= d; k++) {
	  sprintf (p, "+%ld", xx[k]);
	  while (*p)
	    p++;
	}
	v->name = gb_save_string (&buffer[1]);
	hash_in (v);

      }

      if (d < max_parts) {
	xx[d + 1] = 0;
	for (j = 1; j <= d; j++) {
	  if (xx[j] != xx[j + 1]) {
	    long a, b;
	    for (b = xx[j] / 2, a = xx[j] - b; b; a++, b--) {
	      register Vertex *u;
	      register char *p = buffer;
	      for (k = j + 1; xx[k] > a; k++)
		nn[k - 1] = xx[k];
	      nn[k - 1] = a;
	      for (; xx[k] > b; k++)
		nn[k] = xx[k];
	      nn[k] = b;
	      for (; k <= d; k++)
		nn[k + 1] = xx[k];
	      for (k = 1; k <= d + 1; k++) {
		sprintf (p, "+%ld", nn[k]);
		while (*p)
		  p++;
	      }
	      u = hash_out (&buffer[1]);
	      if (u == NULL)
		panic (impossible + 2);
	      if (directed)
		gb_new_arc (v, u, 1L);
	      else
		gb_new_edge (v, u, 1L);
	    }

	  }
	  nn[j] = xx[j];
	}
      }

      v++;

      if (d == 1)
	goto last;
      for (k = d - 1;; k--) {
	if (xx[k] < sig[k] && xx[k] < xx[k - 1])
	  break;
	if (k == 1)
	  goto last;
      }
      xx[k]++;

    }
  }
last:if (v != new_graph->vertices + new_graph->n)
    panic (impossible);

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);

  }
  return new_graph;
}

Graph *
binary (n, max_height, directed)
     unsigned long n;
     unsigned long max_height;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  if (2 * n + 2 > BUF_SIZE)
    panic (bad_specs);
  if (max_height == 0 || max_height > n)
    max_height = n;
  if (max_height > 30)
    panic (very_bad_specs);

  {
    long nverts;
    if (n >= 20 && max_height >= 6) {
      register float ss;
      d = (1L << max_height) - 1 - n;
      if (d > 8)
	panic (bad_specs + 1);
      if (d < 0)
	nverts = 0;
      else {
	nn[0] = nn[1] = 1;
	for (k = 2; k <= d; k++)
	  nn[k] = 0;
	for (j = 2; j <= max_height; j++) {
	  for (k = d; k; k--) {
	    for (ss = 0.0, i = k; i >= 0; i--)
	      ss += ((float) nn[i]) * ((float) nn[k - i]);
	    if (ss > MAX_NNN)
	      panic (very_bad_specs + 1);
	    for (s = 0, i = k; i >= 0; i--)
	      s += nn[i] * nn[k - i];
	    nn[k] = s;
	  }
	  i = (1L << j) - 1;
	  if (i <= d)
	    nn[i]++;
	}
	nverts = nn[d];
      }
    }

    else {
      nn[0] = nn[1] = 1;
      for (k = 2; k <= n; k++)
	nn[k] = 0;
      for (j = 2; j <= max_height; j++)
	for (k = n - 1; k; k--) {
	  for (s = 0, i = k; i >= 0; i--)
	    s += nn[i] * nn[k - i];
	  nn[k + 1] = s;
	}
      nverts = nn[n];
    }
    new_graph = gb_new_graph (nverts);
    if (new_graph == NULL)
      panic (no_room);
    sprintf (new_graph->id, "binary(%lu,%lu,%d)",
	     n, max_height, directed ? 1 : 0);
    strcpy (new_graph->util_types, "VVZZZZZZZZZZZZ");
  }

  {
    register long *xtab, *ytab, *ltab, *stab;

    xtab = gb_typed_alloc (8 * n + 4, long, working_storage);
    if (gb_trouble_code) {
      gb_recycle (new_graph);
      panic (no_room + 2);
    }
    d = n + n;
    ytab = xtab + (d + 1);
    ltab = ytab + (d + 1);
    stab = ltab + (d + 1);
    ltab[0] = 1L << max_height;
    stab[0] = n;

    v = new_graph->vertices;
    if (ltab[0] > n) {
      k = 0;
      xtab[0] = n ? 1 : 0;
      while (1) {

	for (j = k + 1; j <= d; j++) {
	  if (xtab[j - 1]) {
	    ltab[j] = ltab[j - 1] >> 1;
	    ytab[j] = ytab[j - 1] + ltab[j];
	    stab[j] = stab[j - 1];
	  }
	  else {
	    ytab[j] = ytab[j - 1] & (ytab[j - 1] - 1);
	    ltab[j] = ytab[j - 1] - ytab[j];
	    stab[j] = stab[j - 1] - 1;
	  }
	  if (stab[j] <= ytab[j])
	    xtab[j] = 0;
	  else
	    xtab[j] = 1;
	}

	;

	{
	  register char *p = buffer;
	  for (k = 0; k <= d; k++, p++)
	    *p = (xtab[k] ? '.' : 'x');
	  v->name = gb_save_string (buffer);
	  hash_in (v);

	}

	;

	for (j = 0; j < d; j++)
	  if (xtab[j] == 1 && xtab[j + 1] == 1) {
	    for (i = j + 1, s = 0; s >= 0; s += (xtab[i + 1] << 1) - 1, i++)
	      xtab[i] = xtab[i + 1];
	    xtab[i] = 1;
	    {
	      register char *p = buffer;
	      register Vertex *u;
	      for (k = 0; k <= d; k++, p++)
		*p = (xtab[k] ? '.' : 'x');
	      u = hash_out (buffer);
	      if (u) {
		if (directed)
		  gb_new_arc (v, u, 1L);
		else
		  gb_new_edge (v, u, 1L);
	      }
	    }
	    for (i--; i > j; i--)
	      xtab[i + 1] = xtab[i];
	    xtab[i + 1] = 1;
	  }

	;
	v++;

	for (k = d - 1;; k--) {
	  if (k <= 0)
	    goto last;
	  if (xtab[k])
	    break;
	}
	for (k--;; k--) {
	  if (xtab[k] == 0 && ltab[k] > 1)
	    break;
	  if (k == 0)
	    goto last;
	}
	xtab[k]++;

	;
      }
    }
  }
last:if (v != new_graph->vertices + new_graph->n)
    panic (impossible);
  gb_free (working_storage);

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
complement (g, copy, self, directed)
     Graph *g;
     long copy;
     long self;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  register long n;
  register Vertex *u;
  register siz_t delta;
  if (g == NULL)
    panic (missing_operand);

  n = g->n;
  new_graph = gb_new_graph (n);
  if (new_graph == NULL)
    panic (no_room);
  delta = ((siz_t) (new_graph->vertices)) - ((siz_t) (g->vertices));
  for (u = new_graph->vertices, v = g->vertices; v < g->vertices + n; u++, v++)
    u->name = gb_save_string (v->name);

  sprintf (buffer, ",%d,%d,%d)", copy ? 1 : 0, self ? 1 : 0, directed ? 1 : 0);
  make_compound_id (new_graph, "complement(", g, buffer);

  for (v = g->vertices; v < g->vertices + n; v++) {
    register Vertex *vv;
    u = vert_offset (v, delta);

    {
      register Arc *a;
      for (a = v->arcs; a; a = a->next)
	vert_offset (a->tip, delta)->tmp = u;
    }
    if (directed) {
      for (vv = new_graph->vertices; vv < new_graph->vertices + n; vv++)
	if ((vv->tmp == u && copy) || (vv->tmp != u && !copy))
	  if (vv != u || self)
	    gb_new_arc (u, vv, 1L);
    }
    else {
      for (vv = (self ? u : u + 1); vv < new_graph->vertices + n; vv++)
	if ((vv->tmp == u && copy) || (vv->tmp != u && !copy))
	  gb_new_edge (u, vv, 1L);
    }
  }
  for (v = new_graph->vertices; v < new_graph->vertices + n; v++)
    v->tmp = NULL;

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);

  }
  return new_graph;
}

Graph *
gunion (g, gg, multi, directed)
     Graph *g, *gg;
     long multi;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  register long n;
  register Vertex *u;
  register siz_t delta, ddelta;
  if (g == NULL || gg == NULL)
    panic (missing_operand);

  n = g->n;
  new_graph = gb_new_graph (n);
  if (new_graph == NULL)
    panic (no_room);
  delta = ((siz_t) (new_graph->vertices)) - ((siz_t) (g->vertices));
  for (u = new_graph->vertices, v = g->vertices; v < g->vertices + n; u++, v++)
    u->name = gb_save_string (v->name);

  sprintf (buffer, ",%d,%d)", multi ? 1 : 0, directed ? 1 : 0);
  make_double_compound_id (new_graph, "gunion(", g, ",", gg, buffer);
  ddelta = ((siz_t) (new_graph->vertices)) -
    ((siz_t) (gg->vertices));

  for (v = g->vertices; v < g->vertices + n; v++) {
    register Arc *a;
    register Vertex *vv = vert_offset (v, delta);

    register Vertex *vvv = vert_offset (vv, -ddelta);

    for (a = v->arcs; a; a = a->next) {
      u = vert_offset (a->tip, delta);

      {
	register Arc *b;
	if (directed) {
	  if (multi || u->tmp != vv)
	    gb_new_arc (vv, u, a->len);
	  else {
	    b = u->tlen;
	    if (a->len < b->len)
	      b->len = a->len;
	  }
	  u->tmp = vv;
	  u->tlen = vv->arcs;
	}
	else if (u >= vv) {
	  if (multi || u->tmp != vv)
	    gb_new_edge (vv, u, a->len);
	  else {
	    b = u->tlen;
	    if (a->len < b->len)
	      b->len = (b + 1)->len = a->len;
	  }
	  u->tmp = vv;
	  u->tlen = vv->arcs;
	  if (u == vv && a->next == a + 1)
	    a++;
	}
      }

    }
    if (vvv < gg->vertices + gg->n)
      for (a = vvv->arcs; a; a = a->next) {
	u = vert_offset (a->tip, ddelta);
	if (u < new_graph->vertices + n) {
	  register Arc *b;
	  if (directed) {
	    if (multi || u->tmp != vv)
	      gb_new_arc (vv, u, a->len);
	    else {
	      b = u->tlen;
	      if (a->len < b->len)
		b->len = a->len;
	    }
	    u->tmp = vv;
	    u->tlen = vv->arcs;
	  }
	  else if (u >= vv) {
	    if (multi || u->tmp != vv)
	      gb_new_edge (vv, u, a->len);
	    else {
	      b = u->tlen;
	      if (a->len < b->len)
		b->len = (b + 1)->len = a->len;
	    }
	    u->tmp = vv;
	    u->tlen = vv->arcs;
	    if (u == vv && a->next == a + 1)
	      a++;
	  }
	}

	;
      }
  }
  for (v = new_graph->vertices; v < new_graph->vertices + n; v++)
    v->tmp = NULL, v->tlen = NULL;

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
intersection (g, gg, multi, directed)
     Graph *g, *gg;
     long multi;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  register long n;
  register Vertex *u;
  register siz_t delta, ddelta;
  if (g == NULL || gg == NULL)
    panic (no_room + 1);

  n = g->n;
  new_graph = gb_new_graph (n);
  if (new_graph == NULL)
    panic (no_room);
  delta = ((siz_t) (new_graph->vertices)) - ((siz_t) (g->vertices));
  for (u = new_graph->vertices, v = g->vertices; v < g->vertices + n; u++, v++)
    u->name = gb_save_string (v->name);

  sprintf (buffer, ",%d,%d)", multi ? 1 : 0, directed ? 1 : 0);
  make_double_compound_id (new_graph, "intersection(", g, ",", gg, buffer);
  ddelta = ((siz_t) (new_graph->vertices)) -
    ((siz_t) (gg->vertices));

  for (v = g->vertices; v < g->vertices + n; v++) {
    register Arc *a;
    register Vertex *vv = vert_offset (v, delta);

    register Vertex *vvv = vert_offset (vv, -ddelta);

    if (vvv >= gg->vertices + gg->n)
      continue;

    for (a = v->arcs; a; a = a->next) {
      u = vert_offset (a->tip, delta);
      if (u->tmp == vv) {
	u->mult++;
	if (a->len < u->minlen)
	  u->minlen = a->len;
      }
      else
	u->tmp = vv, u->mult = 0, u->minlen = a->len;
      if (u == vv && !directed && a->next == a + 1)
	a++;

    }

    for (a = vvv->arcs; a; a = a->next) {
      u = vert_offset (a->tip, ddelta);
      if (u >= new_graph->vertices + n)
	continue;
      if (u->tmp == vv) {
	long l = u->minlen;
	if (a->len > l)
	  l = a->len;
	if (u->mult < 0) {
	  register Arc *b = u->tlen;
	  if (l < b->len) {
	    b->len = l;
	    if (!directed)
	      (b + 1)->len = l;
	  }
	}

	else {
	  if (directed)
	    gb_new_arc (vv, u, l);
	  else {
	    if (vv <= u)
	      gb_new_edge (vv, u, l);
	    if (vv == u && a->next == a + 1)
	      a++;
	  }
	  if (!multi) {
	    u->tlen = vv->arcs;
	    u->mult = -1;
	  }
	  else if (u->mult == 0)
	    u->tmp = NULL;
	  else
	    u->mult--;
	}

	;
      }
    }
  }

  for (v = new_graph->vertices; v < new_graph->vertices + n; v++) {
    v->tmp = NULL;
    v->tlen = NULL;
    v->mult = 0;
    v->minlen = 0;
  }

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
lines (g, directed)
     Graph *g;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  register long m;
  register Vertex *u;
  if (g == NULL)
    panic (no_room + 1);

  m = (directed ? g->m : (g->m) / 2);
  new_graph = gb_new_graph (m);
  if (new_graph == NULL)
    panic (no_room);
  make_compound_id (new_graph, "lines(", g, directed ? ",1)" : ",0)");
  u = new_graph->vertices;
  for (v = g->vertices + g->n - 1; v >= g->vertices; v--) {
    register Arc *a;
    register long mapped = 0;
    for (a = v->arcs; a; a = a->next) {
      register Vertex *vv = a->tip;
      if (!directed) {
	if (vv < v)
	  continue;
	if (vv >= g->vertices + g->n)
	  goto near_panic;

      }

      u->u.V = v;
      u->v.V = vv;
      u->w.A = a;
      if (!directed) {
	if (u >= new_graph->vertices + m || (a + 1)->tip != v)
	  goto near_panic;
	if (v == vv && a->next == a + 1)
	  a++;
	else
	  (a + 1)->tip = u;
      }
      sprintf (buffer, "%.*s-%c%.*s", (BUF_SIZE - 3) / 2, v->name,
	       directed ? '>' : '-', BUF_SIZE / 2 - 1, vv->name);
      u->name = gb_save_string (buffer);

      if (!mapped) {
	u->map = v->map;

	v->map = u;
	mapped = 1;
      }
      u++;
    }
  }
  if (u != new_graph->vertices + m)
    goto near_panic;

  if (directed)
    for (u = new_graph->vertices; u < new_graph->vertices + m; u++) {
      v = u->v.V;
      if (v->arcs) {
	v = v->map;
	do {
	  gb_new_arc (u, v, 1L);
	  v++;
	} while (v->u.V == u->v.V);
      }
    }

  else
    for (u = new_graph->vertices; u < new_graph->vertices + m; u++) {
      register Vertex *vv;
      register Arc *a;
      register long mapped = 0;
      v = u->u.V;
      for (vv = v->map; vv < u; vv++)
	gb_new_edge (u, vv, 1L);
      v = u->v.V;
      for (a = v->arcs; a; a = a->next) {
	vv = a->tip;
	if (vv < u && vv >= new_graph->vertices)
	  gb_new_edge (u, vv, 1L);
	else if (vv >= v && vv < g->vertices + g->n)
	  mapped = 1;
      }
      if (mapped && v > u->u.V)
	for (vv = v->map; vv->u.V == v; vv++)
	  gb_new_edge (u, vv, 1L);
    }

  for (u = new_graph->vertices, v = NULL; u < new_graph->vertices + m; u++) {
    if (u->u.V != v) {
      v = u->u.V;
      v->map = u->map;
      u->map = NULL;
    }
    if (!directed)
      ((u->w.A) + 1)->tip = v;
  }

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
near_panic:

  m = u - new_graph->vertices;

  for (u = new_graph->vertices, v = NULL; u < new_graph->vertices + m; u++) {
    if (u->u.V != v) {
      v = u->u.V;
      v->map = u->map;
      u->map = NULL;
    }
    if (!directed)
      ((u->w.A) + 1)->tip = v;
  }

  gb_recycle (new_graph);
  panic (invalid_operand);

}

Graph *
product (g, gg, type, directed)
     Graph *g, *gg;
     long type;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  register Vertex *u, *vv;
  register long n;
  if (g == NULL || gg == NULL)
    panic (no_room + 1);

  {
    float test_product = ((float) (g->n)) * ((float) (gg->n));
    if (test_product > MAX_NNN)
      panic (very_bad_specs);
  }
  n = (g->n) * (gg->n);
  new_graph = gb_new_graph (n);
  if (new_graph == NULL)
    panic (no_room);
  for (u = new_graph->vertices, v = g->vertices, vv = gg->vertices;
       u < new_graph->vertices + n; u++) {
    sprintf (buffer, "%.*s,%.*s", BUF_SIZE / 2 - 1, v->name, (BUF_SIZE - 1) / 2, vv->name);
    u->name = gb_save_string (buffer);
    if (++vv == gg->vertices + gg->n)
      vv = gg->vertices, v++;
  }
  sprintf (buffer, ",%d,%d)", (type ? 2 : 0) - (int) (type & 1), directed ? 1 : 0);
  make_double_compound_id (new_graph, "product(", g, ",", gg, buffer);

  if ((type & 1) == 0) {
    register Vertex *uu, *uuu;
    register Arc *a;
    register siz_t delta;
    delta = ((siz_t) (new_graph->vertices)) - ((siz_t) (gg->vertices));
    for (u = gg->vertices; u < gg->vertices + gg->n; u++)
      for (a = u->arcs; a; a = a->next) {
	v = a->tip;
	if (!directed) {
	  if (u > v)
	    continue;
	  if (u == v && a->next == a + 1)
	    a++;
	}
	for (uu = vert_offset (u, delta), vv = vert_offset (v, delta);
	     uu < new_graph->vertices + n; uu += gg->n, vv += gg->n)
	  if (directed)
	    gb_new_arc (uu, vv, a->len);
	  else
	    gb_new_edge (uu, vv, a->len);
      }

    for (u = g->vertices, uu = new_graph->vertices; uu < new_graph->vertices + n;
	 u++, uu += gg->n)
      for (a = u->arcs; a; a = a->next) {
	v = a->tip;
	if (!directed) {
	  if (u > v)
	    continue;
	  if (u == v && a->next == a + 1)
	    a++;
	}
	vv = new_graph->vertices + ((gg->n) * (v - g->vertices));
	for (uuu = uu; uuu < uu + gg->n; uuu++, vv++)
	  if (directed)
	    gb_new_arc (uuu, vv, a->len);
	  else
	    gb_new_edge (uuu, vv, a->len);
      }

  }

  if (type) {
    Vertex *uu;
    Arc *a;
    siz_t delta0 =
    ((siz_t) (new_graph->vertices)) - ((siz_t) (gg->vertices));
    siz_t del = (gg->n) * sizeof (Vertex);
    register siz_t delta, ddelta;
    for (uu = g->vertices, delta = delta0; uu < g->vertices + g->n; uu++, delta += del)
      for (a = uu->arcs; a; a = a->next) {
	vv = a->tip;
	if (!directed) {
	  if (uu > vv)
	    continue;
	  if (uu == vv && a->next == a + 1)
	    a++;
	}
	ddelta = delta0 + del * (vv - g->vertices);
	for (u = gg->vertices; u < gg->vertices + gg->n; u++) {
	  register Arc *aa;
	  for (aa = u->arcs; aa; aa = aa->next) {
	    long length = a->len;
	    if (length > aa->len)
	      length = aa->len;
	    v = aa->tip;
	    if (directed)
	      gb_new_arc (vert_offset (u, delta), vert_offset (v, ddelta), length);
	    else
	      gb_new_edge (vert_offset (u, delta),
			   vert_offset (v, ddelta), length);
	  }
	}
      }
  }

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);

  }
  return new_graph;
}

Graph *
induced (g, description, self, multi, directed)
     Graph *g;
     char *description;
     long self;
     long multi;
     long directed;
{

  Graph *new_graph;
  register long i, j, k;
  register long d;
  register Vertex *v;
  register long s;

  register Vertex *u;
  register long n = 0;
  register long nn = 0;
  if (g == NULL)
    panic (missing_operand);

  for (v = g->vertices; v < g->vertices + g->n; v++)
    if (v->ind > 0) {
      if (n > IND_GRAPH)
	panic (very_bad_specs);
      if (v->ind >= IND_GRAPH) {
	if (v->subst == NULL)
	  panic (missing_operand + 1);

	n += v->subst->n;
      }
      else
	n += v->ind;
    }
    else if (v->ind < -nn)
      nn = -(v->ind);
  if (n > IND_GRAPH || nn > IND_GRAPH)
    panic (very_bad_specs + 1);
  n += nn;

  new_graph = gb_new_graph (n);
  if (new_graph == NULL)
    panic (no_room);

  for (k = 1, u = new_graph->vertices; k <= nn; k++, u++) {
    u->mult = -k;
    sprintf (buffer, "%ld", -k);
    u->name = gb_save_string (buffer);
  }
  for (v = g->vertices; v < g->vertices + g->n; v++)
    if ((k = v->ind) < 0)
      v->map = (new_graph->vertices) - (k + 1);
    else if (k > 0) {
      u->mult = k;
      v->map = u;
      if (k <= 2) {
	u->name = gb_save_string (v->name);
	u++;
	if (k == 2) {
	  sprintf (buffer, "%s'", v->name);
	  u->name = gb_save_string (buffer);
	  u++;
	}
      }
      else if (k >= IND_GRAPH) {
	register Graph *gg = v->subst;
	register Vertex *vv = gg->vertices;
	register Arc *a;
	siz_t delta = ((siz_t) u) - ((siz_t) vv);
	for (j = 0; j < v->subst->n; j++, u++, vv++) {
	  sprintf (buffer, "%.*s:%.*s", BUF_SIZE / 2 - 1, v->name, (BUF_SIZE - 1) / 2, vv->name);
	  u->name = gb_save_string (buffer);
	  for (a = vv->arcs; a; a = a->next) {
	    register Vertex *vvv = a->tip;
	    Vertex *uu = vert_offset (vvv, delta);
	    if (vvv == vv && !self)
	      continue;
	    if (uu->tmp == u && !multi) {
	      register Arc *b = uu->tlen;
	      if (a->len < b->len) {
		b->len = a->len;
		if (!directed)
		  (b + 1)->len = a->len;
	      }
	      continue;
	    }

	    if (!directed) {
	      if (vvv < vv)
		continue;
	      if (vvv == vv && a->next == a + 1)
		a++;
	      gb_new_edge (u, uu, a->len);
	    }
	    else
	      gb_new_arc (u, uu, a->len);
	    uu->tmp = u;
	    uu->tlen = ((directed || u <= uu) ? u->arcs : uu->arcs);
	  }
	}
      }

      else
	for (j = 0; j < k; j++, u++) {
	  sprintf (buffer, "%.*s:%ld", BUF_SIZE - 12, v->name, j);
	  u->name = gb_save_string (buffer);
	}
    }

  sprintf (buffer, ",%s,%d,%d,%d)", description ? description : null_string,
	   self ? 1 : 0, multi ? 1 : 0, directed ? 1 : 0);
  make_compound_id (new_graph, "induced(", g, buffer);

  for (v = g->vertices; v < g->vertices + g->n; v++) {
    u = v->map;
    if (u) {
      register Arc *a;
      register Vertex *uu, *vv;
      k = u->mult;
      if (k < 0)
	k = 1;
      else if (k >= IND_GRAPH)
	k = v->subst->n;
      for (; k; k--, u++) {
	if (!multi)
	  for (a = u->arcs; a; a = a->next) {
	    a->tip->tmp = u;
	    if (directed || a->tip > u || a->next == a + 1)
	      a->tip->tlen = a;
	    else
	      a->tip->tlen = a + 1;
	  }

	;
	for (a = v->arcs; a; a = a->next) {
	  vv = a->tip;
	  uu = vv->map;
	  if (uu == NULL)
	    continue;
	  j = uu->mult;
	  if (j < 0)
	    j = 1;
	  else if (j >= IND_GRAPH)
	    j = vv->subst->n;
	  if (!directed) {
	    if (vv < v)
	      continue;
	    if (vv == v) {
	      if (a->next == a + 1)
		a++;
	      j = k, uu = u;
	    }
	  }

	  for (; j; j--, uu++) {
	    if (u == uu && !self)
	      continue;
	    if (uu->tmp == u && !multi) {
	      register Arc *b = uu->tlen;
	      if (a->len < b->len) {
		b->len = a->len;
		if (!directed)
		  (b + 1)->len = a->len;
	      }
	      continue;
	    }

	    if (directed)
	      gb_new_arc (u, uu, a->len);
	    else
	      gb_new_edge (u, uu, a->len);
	    uu->tmp = u;
	    uu->tlen = ((directed || u <= uu) ? u->arcs : uu->arcs);
	  }

	}
      }
    }
  }

  for (v = g->vertices; v < g->vertices + g->n; v++)
    if (v->map)
      v->ind = v->map->mult;
  for (v = new_graph->vertices; v < new_graph->vertices + n; v++)
    v->u.I = v->v.I = v->z.I = 0;

  if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  return new_graph;
}

Graph *
bi_complete (n1, n2, directed)
     unsigned long n1;
     unsigned long n2;
     long directed;
{
  Graph *new_graph = board (2L, 0L, 0L, 0L, 1L, 0L, directed);
  if (new_graph) {
    new_graph->vertices->ind = n1;
    (new_graph->vertices + 1)->ind = n2;
    new_graph = induced (new_graph, NULL, 0L, 0L, directed);
    if (new_graph) {
      sprintf (new_graph->id, "bi_complete(%lu,%lu,%d)",
	       n1, n2, directed ? 1 : 0);
      mark_bipartite (new_graph, n1);
    }
  }
  return new_graph;
}

Graph *
wheel (n, n1, directed)
     unsigned long n;
     unsigned long n1;
     long directed;
{
  Graph *new_graph = board (2L, 0L, 0L, 0L, 1L, 0L, directed);

  if (new_graph) {
    new_graph->vertices->ind = n1;
    (new_graph->vertices + 1)->ind = IND_GRAPH;
    (new_graph->vertices + 1)->subst = board (n, 0L, 0L, 0L, 1L, 1L, directed);

    new_graph = induced (new_graph, NULL, 0L, 0L, directed);
    if (new_graph) {
      sprintf (new_graph->id, "wheel(%lu,%lu,%d)",
	       n, n1, directed ? 1 : 0);
    }
  }
  return new_graph;
}

#define gb_typed_alloc(n,t,s)(t*)gb_alloc((long)((n)*sizeof(t)),s) \

#define n_1 uu.I \

#define arcs_per_block 102 \

#define gb_new_graph gb_nugraph
#define gb_new_arc gb_nuarc
#define gb_new_edge gb_nuedge \

#define string_block_size 1016 \

#define hash_link u.V
#define hash_head v.V \

#define HASH_MULT 314159
#define HASH_PRIME 516595003 \

#undef tmp

static Arc *next_arc;
static Arc *bad_arc;
static char *next_string;
static char *bad_string;
static Arc dummy_arc[2];
static Graph dummy_graph;
static Graph *cur_graph = &dummy_graph;

long verbose = 0;
long panic_code = 0;

long gb_trouble_code = 0;

long extra_n = 4;
char null_string[1];

siz_t edge_trick = sizeof (Arc) - (sizeof (Arc) & (sizeof (Arc) - 1));

char *
gb_alloc (n, s)
     long n;
     Area s;
{
  long m = sizeof (char *);
  Area t;
  char *loc;
  if (n <= 0 || n > 0xffff00 - 2 * m) {
    gb_trouble_code |= 2;
    return NULL;
  }
  n = ((n + m - 1) / m) * m;
  loc = (char *) calloc ((unsigned) ((n + 2 * m + 255) / 256), 256);
  if (loc) {
    *t = (struct area_pointers *) (loc + n);
    (*t)->first = loc;
    (*t)->next = *s;
    *s = *t;
  }
  else
    gb_trouble_code |= 1;
  return loc;
}

void 
gb_free (s)
     Area s;
{
  Area t;
  while (*s) {
    *t = (*s)->next;
    free ((*s)->first);
    *s = *t;
  }
}

Graph *
gb_new_graph (n)
     long n;
{
  cur_graph = (Graph *) calloc (1, sizeof (Graph));
  if (cur_graph) {
    cur_graph->vertices = gb_typed_alloc (n + extra_n, Vertex, cur_graph->data);
    if (cur_graph->vertices) {
      Vertex *p;
      cur_graph->n = n;
      for (p = cur_graph->vertices + n + extra_n - 1; p >= cur_graph->vertices; p--)
	p->name = null_string;
      sprintf (cur_graph->id, "gb_new_graph(%ld)", n);
      strcpy (cur_graph->util_types, "ZZZZZZZZZZZZZZ");
    }
    else {
      free ((char *) cur_graph);
      cur_graph = NULL;
    }
  }
  next_arc = bad_arc = NULL;
  next_string = bad_string = NULL;
  gb_trouble_code = 0;
  return cur_graph;
}

void 
make_compound_id (g, s1, gg, s2)
     Graph *g;
     char *s1;
     Graph *gg;
     char *s2;
{
  int avail = ID_FIELD_SIZE - strlen (s1) - strlen (s2);
  char tmp[ID_FIELD_SIZE];
  strcpy (tmp, gg->id);
  if (strlen (tmp) < avail)
    sprintf (g->id, "%s%s%s", s1, tmp, s2);
  else
    sprintf (g->id, "%s%.*s...)%s", s1, avail - 5, tmp, s2);
}

void 
make_double_compound_id (g, s1, gg, s2, ggg, s3)
     Graph *g;
     char *s1;
     Graph *gg;
     char *s2;
     Graph *ggg;
     char *s3;
{
  int avail = ID_FIELD_SIZE - strlen (s1) - strlen (s2) - strlen (s3);
  if (strlen (gg->id) + strlen (ggg->id) < avail)
    sprintf (g->id, "%s%s%s%s%s", s1, gg->id, s2, ggg->id, s3);
  else
    sprintf (g->id, "%s%.*s...)%s%.*s...)%s", s1, avail / 2 - 5, gg->id,
	     s2, (avail - 9) / 2, ggg->id, s3);
}

Arc *
gb_virgin_arc ()
{
  register Arc *cur_arc = next_arc;
  if (cur_arc == bad_arc) {
    cur_arc = gb_typed_alloc (arcs_per_block, Arc, cur_graph->data);
    if (cur_arc == NULL)
      cur_arc = dummy_arc;
    else {
      next_arc = cur_arc + 1;
      bad_arc = cur_arc + arcs_per_block;
    }
  }
  else
    next_arc++;
  return cur_arc;
}

void 
gb_new_arc (u, v, len)
     Vertex *u, *v;
     long len;
{
  register Arc *cur_arc = gb_virgin_arc ();
  cur_arc->tip = v;
  cur_arc->next = u->arcs;
  cur_arc->len = len;
  u->arcs = cur_arc;
  cur_graph->m++;
}

void 
gb_new_edge (u, v, len)
     Vertex *u, *v;
     long len;
{
  register Arc *cur_arc = gb_virgin_arc ();
  if (cur_arc != dummy_arc)
    next_arc++;
  if (u < v) {
    cur_arc->tip = v;
    cur_arc->next = u->arcs;
    (cur_arc + 1)->tip = u;
    (cur_arc + 1)->next = v->arcs;
    u->arcs = cur_arc;
    v->arcs = cur_arc + 1;
  }
  else {
    (cur_arc + 1)->tip = v;
    (cur_arc + 1)->next = u->arcs;
    u->arcs = cur_arc + 1;
    cur_arc->tip = u;
    cur_arc->next = v->arcs;
    v->arcs = cur_arc;
  }
  cur_arc->len = (cur_arc + 1)->len = len;
  cur_graph->m += 2;
}

char *
gb_save_string (s)
     register char *s;
{
  register char *p = s;
  register long len;

  while (*p++);
  len = p - s;
  p = next_string;
  if (p + len > bad_string) {
    long size = string_block_size;
    if (len > size)
      size = len;
    p = gb_alloc (size, cur_graph->data);
    if (p == NULL)
      return null_string;
    bad_string = p + size;
  }
  while (*s)
    *p++ = *s++;
  *p++ = '\0';
  next_string = p;
  return p - len;
}

void 
switch_to_graph (g)
     Graph *g;
{
  cur_graph->ww.A = next_arc;
  cur_graph->xx.A = bad_arc;
  cur_graph->yy.S = next_string;
  cur_graph->zz.S = bad_string;
  cur_graph = (g ? g : &dummy_graph);
  next_arc = cur_graph->ww.A;
  bad_arc = cur_graph->xx.A;
  next_string = cur_graph->yy.S;
  bad_string = cur_graph->zz.S;
  cur_graph->ww.A = NULL;
  cur_graph->xx.A = NULL;
  cur_graph->yy.S = NULL;
  cur_graph->zz.S = NULL;
}

void 
gb_recycle (g)
     Graph *g;
{
  if (g) {
    gb_free (g->data);
    gb_free (g->aux_data);
    free ((char *) g);
  }
}

void 
hash_in (v)
     Vertex *v;
{
  register char *t = v->name;
  register Vertex *u;

  {
    register long h;
    for (h = 0; *t; t++) {
      h += (h ^ (h >> 1)) + HASH_MULT * (unsigned char) *t;
      while (h >= HASH_PRIME)
	h -= HASH_PRIME;
    }
    u = cur_graph->vertices + (h % cur_graph->n);
  }

  v->hash_link = u->hash_head;
  u->hash_head = v;
}

Vertex *
hash_out (s)
     char *s;
{
  register char *t = s;
  register Vertex *u;

  {
    register long h;
    for (h = 0; *t; t++) {
      h += (h ^ (h >> 1)) + HASH_MULT * (unsigned char) *t;
      while (h >= HASH_PRIME)
	h -= HASH_PRIME;
    }
    u = cur_graph->vertices + (h % cur_graph->n);
  }

  for (u = u->hash_head; u; u = u->hash_link)
    if (strcmp (s, u->name) == 0)
      return u;
  return NULL;
}

void 
hash_setup (g)
     Graph *g;
{
  Graph *save_cur_graph;
  if (g && g->n > 0) {
    register Vertex *v;
    save_cur_graph = cur_graph;
    cur_graph = g;
    for (v = g->vertices; v < g->vertices + g->n; v++)
      v->hash_head = NULL;
    for (v = g->vertices; v < g->vertices + g->n; v++)
      hash_in (v);
    g->util_types[0] = g->util_types[1] = 'V';

    cur_graph = save_cur_graph;
  }
}

Vertex *
hash_lookup (s, g)
     char *s;
     Graph *g;
{
  Graph *save_cur_graph;
  if (g && g->n > 0) {
    register Vertex *v;
    save_cur_graph = cur_graph;
    cur_graph = g;
    v = hash_out (s);
    cur_graph = save_cur_graph;
    return v;
  }
  else
    return NULL;
}

#define random_graph r_graph
#define random_bigraph r_bigraph
#define random_lengths r_lengths \

#undef panic
#define panic(c){panic_code= c;gb_trouble_code= 0;return NULL;} \

#define dist_code(x)(x?"dist":"0") \

#define rand_len (min_len==max_len?min_len:min_len+gb_unif_rand(max_len-min_len)) \

static char name_buffer[] = "9999999999";

typedef struct {
  long prob;
  long inx;
}

magic_entry;

typedef struct node_struct {
  long key;
  struct node_struct *link;
  long j;
}

node;
static Area temp_nodes;
static node *base_node;

static char buffer[] = "1,-1000000001,-1000000000,dist,1000000000)";

static magic_entry *
walker (n, nn, dist, g)
     long n;
     long nn;
     register long *dist;

     Graph *g;
{
  magic_entry *table;
  long t;
  node *hi = NULL, *lo = NULL;
  register node *p, *q;
  base_node = gb_typed_alloc (nn, node, temp_nodes);
  table = gb_typed_alloc (nn, magic_entry, g->aux_data);
  if (!gb_trouble_code) {

    t = 0x40000000 / nn;
    p = base_node;
    while (nn > n) {
      p->key = 0;
      p->link = lo;
      p->j = --nn;
      lo = p++;
    }
    for (dist = dist + n - 1; n > 0; dist--, p++) {
      p->key = *dist;
      p->j = --n;
      if (*dist > t)
	p->link = hi, hi = p;
      else
	p->link = lo, lo = p;
    }

    while (hi) {
      register magic_entry *r;
      register long x;
      p = hi, hi = p->link;
      q = lo, lo = q->link;
      r = table + q->j;
      x = t * q->j + q->key - 1;
      r->prob = x + x + 1;
      r->inx = p->j;

      if ((p->key -= t - q->key) > t)
	p->link = hi, hi = p;
      else
	p->link = lo, lo = p;
    }

    while (lo) {
      register magic_entry *r;
      register long x;
      q = lo, lo = q->link;
      r = table + q->j;
      x = t * q->j + t - 1;
      r->prob = x + x + 1;

    }

  }
  gb_free (temp_nodes);
  return table;
}

Graph *
random_graph (n, m, multi, self, directed, dist_from, dist_to, min_len, max_len,
	      seed)
     unsigned long n;
     unsigned long m;
     long multi;
     long self;
     long directed;
     long *dist_from;
     long *dist_to;
     long min_len, max_len;
     long seed;
{

  Graph *new_graph;
  long mm;
  register long k;

  long nn = 1;
  long kk = 31;
  magic_entry *from_table, *to_table;

  if (n == 0)
    panic (bad_specs);
  if (min_len > max_len)
    panic (very_bad_specs);
  if (((unsigned long) (max_len)) - ((unsigned long) (min_len)) >=
      ((unsigned long) 0x80000000))
    panic (bad_specs + 1);

  {
    register long acc;
    register long *p;
    if (dist_from) {
      for (acc = 0, p = dist_from; p < dist_from + n; p++) {
	if (*p < 0)
	  panic (invalid_operand);

	if (*p > 0x40000000 - acc)
	  panic (invalid_operand + 1);

	acc += *p;
      }
      if (acc != 0x40000000)
	panic (invalid_operand + 2);
    }
    if (dist_to) {
      for (acc = 0, p = dist_to; p < dist_to + n; p++) {
	if (*p < 0)
	  panic (invalid_operand + 5);

	if (*p > 0x40000000 - acc)
	  panic (invalid_operand + 6);

	acc += *p;
      }
      if (acc != 0x40000000)
	panic (invalid_operand + 7);
    }
  }

  gb_init_rand (seed);

  new_graph = gb_new_graph (n);
  if (new_graph == NULL)
    panic (no_room);
  for (k = 0; k < n; k++) {
    sprintf (name_buffer, "%ld", k);
    (new_graph->vertices + k)->name = gb_save_string (name_buffer);
  }
  sprintf (new_graph->id, "random_graph(%lu,%lu,%d,%d,%d,%s,%s,%ld,%ld,%ld)",
   n, m, multi > 0 ? 1 : multi < 0 ? -1 : 0, self ? 1 : 0, directed ? 1 : 0,
	dist_code (dist_from), dist_code (dist_to), min_len, max_len, seed);

  {
    if (dist_from) {
      while (nn < n)
	nn += nn, kk--;
      from_table = walker (n, nn, dist_from, new_graph);
    }
    if (dist_to) {
      while (nn < n)
	nn += nn, kk--;
      to_table = walker (n, nn, dist_to, new_graph);
    }
    if (gb_trouble_code) {
      gb_recycle (new_graph);
      panic (alloc_fault);
    }
  }

  for (mm = m; mm; mm--) {
    register Vertex *u, *v;
  repeat:
    if (dist_from) {
      register magic_entry *magic;
      register long uu = gb_next_rand ();
      k = uu >> kk;
      magic = from_table + k;
      if (uu <= magic->prob)
	u = new_graph->vertices + k;
      else
	u = new_graph->vertices + magic->inx;
    }

    else
      u = new_graph->vertices + gb_unif_rand (n);
    if (dist_to) {
      register magic_entry *magic;
      register long uu = gb_next_rand ();
      k = uu >> kk;
      magic = to_table + k;
      if (uu <= magic->prob)
	v = new_graph->vertices + k;
      else
	v = new_graph->vertices + magic->inx;
    }

    else
      v = new_graph->vertices + gb_unif_rand (n);
    if (u == v && !self)
      goto repeat;
    if (multi <= 0)
      if (gb_trouble_code)
	goto trouble;
      else {
	register Arc *a;
	long len;
	for (a = u->arcs; a; a = a->next)
	  if (a->tip == v)
	    if (multi == 0)
	      goto repeat;
	    else {
	      len = rand_len;
	      if (len < a->len) {
		a->len = len;
		if (!directed) {
		  if (u <= v)
		    (a + 1)->len = len;
		  else
		    (a - 1)->len = len;
		}
	      }
	      goto done;
	    }
      }

    if (directed)
      gb_new_arc (u, v, rand_len);
    else
      gb_new_edge (u, v, rand_len);
  done:;
  }

trouble:if (gb_trouble_code) {
    gb_recycle (new_graph);
    panic (alloc_fault);
  }
  gb_free (new_graph->aux_data);
  return new_graph;
}

Graph *
random_bigraph (n1, n2, m, multi, dist1, dist2, min_len, max_len, seed)
     unsigned long n1, n2;
     unsigned long m;
     long multi;
     long *dist1, *dist2;
     long min_len, max_len;
     long seed;
{
  unsigned long n = n1 + n2;
  Area new_dists;
  long *dist_from, *dist_to;
  Graph *new_graph;
  init_area (new_dists);
  if (n1 == 0 || n2 == 0)
    panic (bad_specs);
  if (min_len > max_len)
    panic (very_bad_specs);
  if (((unsigned long) (max_len)) - ((unsigned long) (min_len)) >=
      ((unsigned long) 0x80000000))
    panic (bad_specs + 1);
  dist_from = gb_typed_alloc (n, long, new_dists);
  dist_to = gb_typed_alloc (n, long, new_dists);
  if (gb_trouble_code) {
    gb_free (new_dists);
    panic (no_room + 2);
  }

  {
    register long *p, *q;
    register long k;
    p = dist1;
    q = dist_from;
    if (p)
      while (p < dist1 + n1)
	*q++ = *p++;
    else
      for (k = 0; k < n1; k++)
	*q++ = (0x40000000 + k) / n1;
    p = dist2;
    q = dist_to + n1;
    if (p)
      while (p < dist2 + n2)
	*q++ = *p++;
    else
      for (k = 0; k < n2; k++)
	*q++ = (0x40000000 + k) / n2;
  }

  new_graph = random_graph (n, m, multi, 0L, 0L,
			    dist_from, dist_to, min_len, max_len, seed);
  sprintf (new_graph->id, "random_bigraph(%lu,%lu,%lu,%d,%s,%s,%ld,%ld,%ld)",
	   n1, n2, m, multi > 0 ? 1 : multi < 0 ? -1 : 0, dist_code (dist1), dist_code (dist2),
	   min_len, max_len, seed);
  mark_bipartite (new_graph, n1);
  gb_free (new_dists);
  return new_graph;
}

long 
random_lengths (g, directed, min_len, max_len, dist, seed)
     Graph *g;
     long directed;
     long min_len, max_len;
     long *dist;
     long seed;
{
  register Vertex *u, *v;
  register Arc *a;
  long nn = 1, kk = 31;
  magic_entry *dist_table;
  if (g == NULL)
    return missing_operand;
  gb_init_rand (seed);
  if (min_len > max_len)
    return very_bad_specs;
  if (((unsigned long) (max_len)) - ((unsigned long) (min_len)) >=
      ((unsigned long) 0x80000000))
    return bad_specs;

  if (dist) {
    register long acc;
    register long *p;
    register long n = max_len - min_len + 1;
    for (acc = 0, p = dist; p < dist + n; p++) {
      if (*p < 0)
	return -1;
      if (*p > 0x40000000 - acc)
	return 1;
      acc += *p;
    }
    if (acc != 0x40000000)
      return 2;
    while (nn < n)
      nn += nn, kk--;
    dist_table = walker (n, nn, dist, g);
    if (gb_trouble_code) {
      gb_trouble_code = 0;
      return alloc_fault;
    }
  }

  sprintf (buffer, ",%d,%ld,%ld,%s,%ld)", directed ? 1 : 0,
	   min_len, max_len, dist_code (dist), seed);
  make_compound_id (g, "random_lengths(", g, buffer);

  for (u = g->vertices; u < g->vertices + g->n; u++)
    for (a = u->arcs; a; a = a->next) {
      v = a->tip;
      if (directed == 0 && u > v)
	a->len = (a - 1)->len;
      else {
	register long len;
	if (dist == 0)
	  len = rand_len;
	else {
	  long uu = gb_next_rand ();
	  long k = uu >> kk;
	  magic_entry *magic = dist_table + k;
	  if (uu <= magic->prob)
	    len = min_len + k;
	  else
	    len = min_len + magic->inx;
	}
	a->len = len;
	if (directed == 0 && u == v && a->next == a + 1)
	  (++a)->len = len;
      }
    }

  return 0;
}

#define gb_next_rand()(*gb_fptr>=0?*gb_fptr--:gb_flip_cycle()) \

#define mod_diff(x,y)(((x)-(y))&0x7fffffff) \

#define two_to_the_31 ((unsigned long)0x80000000) \

static long A[56] =
{-1};

long *gb_fptr = A;

long 
gb_flip_cycle ()
{
  register long *ii, *jj;
  for (ii = &A[1], jj = &A[32]; jj <= &A[55]; ii++, jj++)
    *ii = mod_diff (*ii, *jj);
  for (jj = &A[1]; ii <= &A[55]; ii++, jj++)
    *ii = mod_diff (*ii, *jj);
  gb_fptr = &A[54];
  return A[55];
}

void 
gb_init_rand (seed)
     long seed;
{
  register long i;
  register long prev = seed, next = 1;
  seed = prev = mod_diff (prev, 0);
  A[55] = prev;
  for (i = 21; i; i = (i + 21) % 55) {
    A[i] = next;

    next = mod_diff (prev, next);
    if (seed & 1)
      seed = 0x40000000 + (seed >> 1);
    else
      seed >>= 1;
    next = mod_diff (next, seed);

    prev = A[i];
  }

  (void) gb_flip_cycle ();
  (void) gb_flip_cycle ();
  (void) gb_flip_cycle ();
  (void) gb_flip_cycle ();
  (void) gb_flip_cycle ();

}

long 
gb_unif_rand (m)
     long m;
{
  register unsigned long t = two_to_the_31 - (two_to_the_31 % m);
  register long r;
  do {
    r = gb_next_rand ();
  } while (t <= (unsigned long) r);
  return r % m;
}
/* --------------------end of gb_lib.c -------------------- */


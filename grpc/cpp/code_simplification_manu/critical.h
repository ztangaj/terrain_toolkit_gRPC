#ifndef CRITICAL_H
#define CRITICAL_H

// cross platflam defining
#ifndef defined ( WIN32 )
#define __func__ __FUNCTION__
#endif
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif


#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <stdexcept>
#include <map>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <list>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <exception>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <array>
#include <stack>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>    /* For command line args */
#include <boost/config.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/astar_search.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace simplify{

#define N 10000000

/*********************************************************************************
** Declarations
**********************************************************************************/
typedef struct Pt {                                                   /* 3D point */
    double x,y,z;
} Point3;
typedef std::vector<Point3> p3list;                                   /* Store 3D points in a vector */

typedef struct Dist {                                                 /* 3D point distances: For Range Query */
    double d;                                                         /* Distance (q,v_id) */
    int v_id;                                                         /* Vertex ID of point */
} Distance;

typedef std::vector<Distance> D;                                      /* Store points, distances in a vector */

struct mycomp {                                                       /* Comparison Functor for sorting */
  bool operator() (Distance i, Distance j) { return (i.d  < j.d); }
} dcompare;


p3list vertex_props;                                                  /* External Vertex Property array */

typedef double cost;

struct VertexProperties {
  std::size_t index;
};

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, VertexProperties,         /* Graph */
                    boost::property<boost::edge_weight_t, cost> > mygraph_t;


typedef boost::property_map<mygraph_t, boost::edge_weight_t>::type WeightMap;       /* Edge Weights */
typedef mygraph_t::vertex_descriptor vert;
typedef mygraph_t::edge_descriptor edge_descriptor;
typedef mygraph_t::vertex_iterator vertex_iterator;
typedef mygraph_t::adjacency_iterator adj_iterator;
typedef std::pair<int, int> edge;


struct triangle {
  std::array<int, 3> v;
  bool operator<(triangle const & o) const {
    return v[0] == o.v[0] ? v[1] == o.v[1] ? v[2] < o.v[2] : v[1] < o.v[1] : v[0] < o.v[0];
  }
  triangle sorted() const {
    auto t = *this;
    std::sort(t.v.begin(), t.v.end());
    return t;
  }
};

std::ostream & operator<<(std::ostream & o, triangle const & t) {
o << " Tri --> (" << t.v[0] << ", " << t.v[1] << ", " << t.v[2] << ")";
return o;
}



/* When Propoerty(2) is satisfied, we must propagate the delta
 * values stored (i.e. u's) in v_bar (vertex to remove) to all
 * it's parents (v's).
 */
typedef struct deltas {
  vert v,u;                                                             /* vertex IDs */
  double delta;                                                         /* delta value */
} delta_pairs;


typedef std::vector<delta_pairs> DELTAS_TO_MOVE;                                /* Holds (v,u,delta_prime) values to update L */

typedef double coord;

std::vector<int> hull_vertices;
coord points[N][2], *P[N+1]; /* an extra position is used */

/*********************************************************************************
** Prototypes
**********************************************************************************/
bool readgraph(std::istream & , mygraph_t &, std::size_t, std::size_t,
               boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type &);

bool net_distance(mygraph_t &g, double & dist, std::size_t s, std::size_t t );

void dump_graph(mygraph_t &g, boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type &id,
                 WeightMap &weightmap);

void dump_graph_visual( std::string fname,
                        std::vector<vert> &v,
                        std::vector< std::pair<vert,vert> > &e );
void dump_graph_visual( std::string fname, mygraph_t &G,
                        boost::property_map< mygraph_t,
                        std::size_t VertexProperties::*>::type &id );


bool remove_vertices( mygraph_t &g, mygraph_t &g_orig,
                      boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type& id,
                      WeightMap &weightmap);

bool triangulate(mygraph_t &g, mygraph_t &small_g, vert v_bar, int degree,
                 boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type& id);

bool check_property_1(mygraph_t &big_g, mygraph_t &tiny_g, vert v_bar, bool triangulated );
bool check_property_2(mygraph_t &g_orig, mygraph_t &g, vert v_bar, DELTAS_TO_MOVE &dmove );

void update_delta_in_L( vert, vert, double );

double get_delta( vert i1, vert i2, mygraph_t &g);
double get_euclidean( vert i1, vert i2);

bool remove_single_vertex( vert v_bar, mygraph_t &g, mygraph_t &g_orig);

void dump_dotfile_graph(mygraph_t &G, boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type &id,
                        WeightMap &weightmap, int idx);
int vertices_left(mygraph_t &g);
bool compute_bounds( double &LB, double &UB, vert &v_i, vert &v_j, mygraph_t &g);

void archive_graph(mygraph_t &G, boost::property_map< mygraph_t,
                   std::size_t VertexProperties::*>::type &id,
                   WeightMap &weightmap, std::ofstream &of);
void get_stats();
bool is_friend( vert u, vert v, vert v_bar);
double estimated_distance(vert v_bar, vert u);

void graph_to_off(mygraph_t &G, boost::property_map< mygraph_t,
                  std::size_t VertexProperties::*>::type &id, std::string outf);
void canonical_triangle( triangle &t );
void DLS( mygraph_t &G, vert start, vert goal, int depth_limit,
          std::vector<vert> &v, std::vector< std::pair<vert,vert> > &e );

/*********************************************************************************
** Class Definitions
**********************************************************************************/

/* euclidean distance heuristic */
template <class Graph, class CostType, class LocMap>
class distance_heuristic : public boost::astar_heuristic<Graph, CostType>
{
public:
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  distance_heuristic(LocMap l, Vertex goal)
    : m_location(l), m_goal(goal) {}
  CostType operator()(Vertex u)
  {
    CostType dx = m_location[m_goal].x - m_location[u].x;
    CostType dy = m_location[m_goal].y - m_location[u].y;
    CostType dz = m_location[m_goal].z - m_location[u].z;

    return ::sqrt(dx*dx + dy*dy + dz*dz );  /* 3D Euclidean */
    // return 0; /* Dijkstra */
  }
private:
  LocMap m_location;
  Vertex m_goal;
};


struct found_goal {}; // exception for termination

/* visitor that terminates when we find the goal */
template <class Vertex>
class astar_goal_visitor : public boost::default_astar_visitor
{
public:
  astar_goal_visitor(Vertex goal) : m_goal(goal) {}
  template <class Graph>
  void examine_vertex(Vertex u, Graph& g) {
    if(u == m_goal)
      throw found_goal();
  }
private:
  Vertex m_goal;
};


int read_points(p3list &vertices)
{
  int n=0;
  for(size_t i=0; i < vertices.size(); i++){
    points[n][0] = vertices[i].x;
    points[n][1] = vertices[i].y;
    P[n] = points[n];
                assert(++n <= N);
  }
  return n;
}

void print_hull(coord **P, int m) {
  int i;
  for (i=0; i<m; i++)
    hull_vertices.push_back((P[i]-points[0])/2);
                // printf("%d ", (P[i]-points[0])/2);
  printf("\n");
}


int ccw(coord **P, int i, int j, int k) {
	coord	a = P[i][0] - P[j][0],
		b = P[i][1] - P[j][1],
		c = P[k][0] - P[j][0],
		d = P[k][1] - P[j][1];
	return a*d - b*c <= 0;	   /* true if points i, j, k counterclockwise */
}


#define CMPM(c,A,B) \
	v = (*(coord**)A)[c] - (*(coord**)B)[c];\
	if (v>0) return 1;\
	if (v<0) return -1;

int cmpl(const void *a, const void *b) {
	double v;
	CMPM(0,a,b);
	CMPM(1,b,a);
	return 0;
}

int cmph(const void *a, const void *b) {return cmpl(b,a);}


int make_chain(coord** V, int n, int (*cmp)(const void*, const void*)) {
	int i, j, s = 1;
	coord* t;

	qsort(V, n, sizeof(coord*), cmp);
	for (i=2; i<n; i++) {
		for (j=s; j>=1 && ccw(V, i, j, j-1); j--){}
		s = j+1;
		t = V[s]; V[s] = V[i]; V[i] = t;
	}
	return s;
}

int ch2d(coord **P, int n)  {
	int u = make_chain(P, n, cmpl);		/* make lower hull */
	if (!n) return 0;
	P[n] = P[0];
	return u+make_chain(P+u, n-u+1, cmph);	/* make upper hull */
}

/*
 * Compute the 3D Euclidean Distance

double euclidean3d(Node *a, Node *b)
{
  double X = (b->x - a->x)*(b->x - a->x);
  double Y = (b->y - a->y)*(b->y - a->y);
  double Z = (b->z - a->z)*(b->z - a->z);
  return sqrt(X+Y+Z);
}
*/

double FindDistanceToSegment( double x1, double y1,
                              double x2, double y2,
                              double pointX, double pointY)
{
    double diffX = x2 - x1;
    float diffY = y2 - y1;
    if ((diffX == 0) && (diffY == 0))
    {
        diffX = pointX - x1;
        diffY = pointY - y1;
        return sqrt(diffX * diffX + diffY * diffY);
    }

    float t = ( (pointX - x1) * diffX + (pointY - y1) * diffY) /
                (diffX * diffX + diffY * diffY);

    if (t < 0)
    {
        //point is nearest to the first point i.e x1 and y1
        diffX = pointX - x1;
        diffY = pointY - y1;
    }
    else if (t > 1)
    {
        //point is nearest to the end point i.e x2 and y2
        diffX = pointX - x2;
        diffY = pointY - y2;
    }
    else
    {
        //if perpendicular line intersect the line segment.
        diffX = pointX - (x1 + t * diffX);
        diffY = pointY - (y1 + t * diffY);
    }

    //returning shortest distance
    return sqrt(diffX * diffX + diffY * diffY);
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             DelPoint;

#include "critical.h"

#define TIMING

#ifdef TIMING
#define INIT_TIMER auto st = std::chrono::high_resolution_clock::now();
#define START_TIMER             st = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << \
    std::chrono::duration_cast<std::chrono::milliseconds>( \
            std::chrono::high_resolution_clock::now()-st \
    ).count() << " ms " << std::endl;
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif


#define PI 3.14159265
#define SIN(x)            (sin ((x)* PI/180))
#define COS(x)            (cos ((x)* PI/180))
#define MIN(a,b)          ((a < b) ?  (a) : (b))
#define MAX(a,b)          ((a < b) ?  (b) : (a))
#define LAMBDA(theta_m)   (MIN( SIN((theta_m))/2, ( SIN((theta_m))*COS((theta_m)) ) )) 

/* --- Error codes --- */
#define E_NONE         0        /* no error */
#define E_NOMEM      (-1)       /* not enough memory */
#define E_FOPEN      (-2)       /* cannot open file */
#define E_FREAD      (-3)       /* read error on file */
#define E_FWRITE     (-4)       /* write error on file */
#define E_OPTION     (-5)       /* unknown option */
#define E_OPTARG     (-6)       /* missing option argument */
#define E_ARGCNT     (-7)       /* too few/many arguments */
#define E_UNKNOWN    (-8)       /* unknown error */

#define E_OUT_OF_BOUNDS (-10)   /* Point is out of extent */
#define SEC_SINCE(t)  ((clock()-(t)) /(double)CLOCKS_PER_SEC)



/*----------------------------------------------------------------------
 Global Variables
 ----------------------------------------------------------------------*/
double beta     = 1.0;
double theta_m  = 0.0;
double lambda   = 0.0;
int REMOVED_VERTICES = 0;

/*
 * L(v) -->  {v-bar1 -> delta1}, {v-bar2 -> delta2} ...
 * A map of map for faster retrieval
 */
typedef std::map<vert, double> inMap;                                                     /* Entries are vertex -> delta */
typedef std::map<vert, inMap>  mainMap;
typedef std::map<vert, std::vector<vert> > hostMap;                                       /* Host entries type */
typedef std::map<vert, int> remMap;

mainMap L;                                                                                /* Main map = L(v) */
remMap REMOVED;                                                                           /* Removed Vertices */
hostMap HST;
std::map<int,int> boundary_points_map;                                                    /* Boundary Vertices which must stay */

void dumpL( mainMap &LL);
void dumpHST(hostMap &HH);
void dumpREMOVED(remMap &REM);

typedef std::pair<vert,vert> ed;
std::vector<ed> ADDBACK_EDGES;

/*----------------------------------------------------------------------
 Debug Logs
 ----------------------------------------------------------------------*/

#define DEBUG_FLAG 0                                                               /* Set to false to stop debugging messages */
#define DEBUG_MORE 0                                                               /* To give much more detailed debug messages */

struct dbglog {
    std::ostream &os_;
    mutable bool has_endl_;
    dbglog (std::ostream &os = std::cerr) : os_(os), has_endl_(false) {}
    ~dbglog () { if (!has_endl_) os_ << std::endl; }
    template <typename T> static bool has_endl (const T &) { return false; }
    static bool has_endl (char c) { return (c == '\n'); }
    static bool has_endl (std::string s) { return has_endl(*s.rbegin()); }
    static bool has_endl (const char *s) { return has_endl(std::string(s)); }
    template <typename T>
    static bool same_manip (T & (*m)(T &), T & (*e)(T &)) { return (m == e); }
    const dbglog & operator << (std::ostream & (*m)(std::ostream &)) const {
        has_endl_ = same_manip(m, std::endl);
        os_ << m;
        return *this;
    }
    template <typename T>
    const dbglog & operator << (const T &v) const {
        has_endl_ = has_endl(v);
        os_ << v;
        return *this;
    }
};


#define D_(msg) if (!DEBUG_ENABLED) {} \
                   else dbglog() << __FILE__ << ":" << __LINE__ << " " << msg

#define DM_(msg) if (!DEBUG_MORE_ENABLED) {} \
                   else dbglog() << __FILE__ << ":" << __LINE__ << " " << msg

#if DEBUG_FLAG
  #define DEBUG_ENABLED 1
#else
  #define DEBUG_ENABLED 0
#endif

#if DEBUG_MORE
  #define DEBUG_MORE_ENABLED 1
#else
  #define DEBUG_MORE_ENABLED 0
#endif

/*----------------------------------------------------------------------
  Functions
 ----------------------------------------------------------------------*/



/*----------------------------------------------------------------------
 * Canonical form of a triangular face is in clockwise manner.
 * Important for the generation of triangles in the OFF file.
 ----------------------------------------------------------------------*/
void canonical_triangle( triangle &t )
{
  DM_("Entering: ") << __func__ << std::endl ;

  vert a = t.v[0];
  vert b = t.v[1];
  vert c = t.v[2];

  double x1 = vertex_props[ a ].x;
  double y1 = vertex_props[ a ].y;

  double x2 = vertex_props[ b ].x;
  double y2 = vertex_props[ b ].y;

  double x3 = vertex_props[ c ].x;
  double y3 = vertex_props[ c ].y;

  /* After translation by (x3,y3) : subtraction from (x1,y1) and (x2,y2)
    and compute determinant of

    | 1         1          1 |
    | x1-x3     y1-y3      0 |
    | x2-x3     y2-y3      0 |

  */
  /* Note: Not caring for case 0 */
  if( ((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3)) < 0 ){ /* Clockwise */
    /* Make this counter clockwise */
    t.v[0] = c;
    t.v[1] = b;
    t.v[2] = a;
  }
  else { /* Do Nothing */}
  DM_("Leaving: ") << __func__ << std::endl ;
}

/*----------------------------------------------------------------------
 * Find the vertex that
 * 1) Find the vertex pair which has a common edge between them
 * 2) Return the vertex which is FURTHER away from our edge (U,V)
 ----------------------------------------------------------------------*/
vert find_further_vertex(std::vector<vert> A_int_B, mygraph_t &G, vert u, vert v)
{
  vert f1=-1,f2=-1;

  for( vert x : A_int_B ) {
    for( vert y : A_int_B ) {
        if(x==y) continue;
        bool edge_exists = boost::edge(x,y,G).second;
        if(edge_exists){
          D_("Found Connecting Edge !! (") << x << "," << y <<")";
          f1 = x; f2 = y;
          break;
        }
    } /* Inner loop */
  }   /* Outer loop */

  /* Get the end points of the line segment */
  double x1 = vertex_props[u].x;
  double y1 = vertex_props[u].y;

  double x2 = vertex_props[v].x;
  double y2 = vertex_props[v].y;

  double f1x = vertex_props[f1].x;
  double f1y = vertex_props[f1].y;

  double f2x = vertex_props[f2].x;
  double f2y = vertex_props[f2].y;

  double d1 = FindDistanceToSegment(x1,y1,x2,y2, f1x,f1y);
  double d2 = FindDistanceToSegment(x1,y1,x2,y2, f2x,f2y);

  return (d1 > d2)? f1 : f2;
}

/*----------------------------------------------------------------------
 * Converts a graph to a .OFF file. We figure out what the faces are
 * in our reduced graph G' in order to write out an OFF file.
 ----------------------------------------------------------------------*/
void graph_to_off(mygraph_t &G, boost::property_map< mygraph_t,
                  std::size_t VertexProperties::*>::type &id, std::string outf)
{
  DM_("Entering: ") << __func__ << std::endl ;

  boost::graph_traits<mygraph_t>::vertex_iterator i, end;
  boost::graph_traits<mygraph_t>::edge_iterator ei, edge_end;
  std::set< vert > A,B;
  triangle t;
  std::set<triangle> s;

  /* Loop through the edges to find triangles in graph */
  for (boost::tie(ei, edge_end) = edges(G); ei != edge_end; ++ei){
    vert u = id[source(*ei, G)];
    vert v = id[target(*ei, G)];
    D_(" **** Processing Edge:--- (") << u << "," << v << ")" << std::endl;

    adj_iterator ai, aend;
    /* Loop through the adjacent vertices of u and form set A */
    for(boost::tie(ai,aend) = boost::adjacent_vertices(u, G);
            ai != aend; ++ai)
    {
      if(*ai != v){
        A.insert(*ai);
      }
    }
    /* Loop through the adjacent vertices of v and form set B */
    for(boost::tie(ai,aend) = boost::adjacent_vertices(v, G);
            ai != aend; ++ai)
    {
      if(*ai != u){
        //  std::cout << "adj(v) ---> "<< *ai << std::endl;
        B.insert(*ai);
      }
    }
    std::vector<vert> A_intersect_B;

    std::set_intersection(A.begin(),A.end(),B.begin(),B.end(),
                          std::back_inserter(A_intersect_B));

    /* Check to see if we have a bad set of triangles forming in G */
    /*********************** ERROR CHECK ****************************************************/
    if(A_intersect_B.size() > 2){
        D_( "[ERROR] ****** Found ") << A_intersect_B.size() << " common vertices!";

        #if DEBUG_FLAG
          std::vector<vert> bad_vertices;
          std::vector<std::pair<vert,vert>> bad_edges;
          bad_vertices.push_back(u);
          bad_vertices.push_back(v);
          bad_edges.push_back(std::make_pair(u,v));
        #endif

        for( vert w : A_intersect_B ) {

          #if DEBUG_FLAG
            bad_vertices.push_back(w);
            bad_edges.push_back(std::make_pair(u,w));
            bad_edges.push_back(std::make_pair(w,v));
          #endif

          t.v[0] = u; t.v[1] = v; t.v[2] = w;
          canonical_triangle(t);
          D_( "BAD TRIANGLE (canonicalized) -- : " ) << t << std::endl;
          /* Remove the triangles that are subsuming other smaller triangles */
          D_( "---------->>> Removing Bad Vertex = ") << w << "\n";
          vert x = find_further_vertex(A_intersect_B, G, u, v);

          /* Remove x from the intersection vector */
          A_intersect_B.erase(
                    std::remove(A_intersect_B.begin(), A_intersect_B.end(), x),
                    A_intersect_B.end());

        } /* End of vert w loop */
        #if DEBUG_FLAG
          dump_graph_visual( "error_triangles", bad_vertices, bad_edges );
        #endif
    }
    /*********************** ERROR CHECK [DONE] ****************************************************/


    for( vert w : A_intersect_B ) {
      t.v[0] = u; t.v[1] = v; t.v[2] = w;
      D_( "TRIANGLE -- : " ) << t << std::endl;
      s.insert(t.sorted()); /* Inserted the triangle with sorted vertices into the set */
    }
    A_intersect_B.clear();
    A.clear(); B.clear();
  } /* End of edge iteration */

  /* Lets dump to a .OFF file */
  std::ofstream off_file;
  off_file.open (outf);
  off_file << std::fixed << std::setprecision(2);  
  int num_of_verts = boost::num_vertices(G) - REMOVED_VERTICES;
  int num_of_faces = s.size();
  int num_of_edges = boost::num_edges(G);

  D_( "(#V, #F, #E) = (") << num_of_verts << ","
                          << num_of_faces << ","
                          << num_of_edges << ")" << std::endl;

  /* Write out header of OFF file */
  off_file << "OFF" << std::endl;
  off_file << num_of_verts << " " << num_of_faces << " " << num_of_edges << std::endl;

  std::map<vert,int> VertexMap;
  int k=0;
  /* Loop through the vertices and fill up the gaps in IDs */
  for(boost::tie(i,end) = boost::vertices(G); i != end; ++i) {
    int degree = boost::degree( *i, G);
    if(degree > 0){
      VertexMap[*i] = k;
      /* Writing out Vertex (x,y,z) to file ... */
      off_file  << vertex_props[*i].x << " "
                << vertex_props[*i].y << " "
                << vertex_props[*i].z << std::endl;
      k++;
    }
  }

  /* Loop through the set of unique triangles */
  for (auto t : s) {
    DM_( "Old: " ) << t << std::endl;
    canonical_triangle(t);
    DM_( "New: " ) << t << std::endl;

    /* Writing out Faces to file ... */
    off_file << "3 "<< VertexMap[ t.v[0] ] << " "
                     << VertexMap[ t.v[1] ] << " "
                     << VertexMap[ t.v[2] ] << std::endl;
  }
  off_file.close();
  DM_("Leaving: ") << __func__ << std::endl ;
} /* End of fn */


/*----------------------------------------------------------------------
 * Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
 * intersect the intersection point may be stored in the doubles i_x and i_y.
 ----------------------------------------------------------------------*/
char get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
    double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y)
{
  DM_("Entering: ") << __func__ << std::endl ;

  double s1_x, s1_y, s2_x, s2_y;
  s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

  double s, t;
  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
  t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1){
    // Collision detected
    if (i_x != NULL)
      *i_x = p0_x + (t * s1_x);
    if (i_y != NULL)
      *i_y = p0_y + (t * s1_y);
    DM_("Leaving: ") << __func__ << std::endl ;
    return 1;
  }
  DM_("Leaving: ") << __func__ << std::endl ;
  return 0; // No collision
}

/*----------------------------------------------------------------------
 * Check if incoming edge is intersecting with any other edge in graph G
 *
 ----------------------------------------------------------------------*/
bool doesIntersect( mygraph_t &G,
                    boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type& id,
                    vert a, vert b )
{
  DM_("Entering: ") << __func__ << std::endl ;

  double A0 = vertex_props[a].x ;
  double B0 = vertex_props[a].y ;

  double A1 = vertex_props[b].x ;
  double B1 = vertex_props[b].y ;

  /* Now we need a list of edges from the graph G */
  boost::graph_traits<mygraph_t>::vertex_iterator i, end;
  boost::graph_traits<mygraph_t>::edge_iterator ei, edge_end;

  /* Loop through the edges to find triangles in graph */
  for (boost::tie(ei, edge_end) = edges(G); ei != edge_end; ++ei){

    vert u = id[source(*ei, G)];
    vert v = id[target(*ei, G)];

    double A2 = vertex_props[u].x;
    double B2 = vertex_props[u].y;
    double A3 = vertex_props[v].x;
    double B3 = vertex_props[v].y;

    /*
     * Check if line ((A0,B0),(A1,B1)) intersects with line ((A2,B2),(A3,B3))
     */
    double ix=0,iy=0;
    if( get_line_intersection( A0,B0,A1,B1, A2,B2,A3,B3, &ix, &iy ) == 1){

        /* Skip end point intersections */
        if( (ix == A0 && iy == B0 ) ||
            (ix == A1 && iy == B1 ) ||
            (ix == A2 && iy == B2 ) ||
            (ix == A3 && iy == B3 ))
          continue;

        D_( "Incoming Edge (")
            << std::fixed << std::setprecision (5) << a << ","<< b << ") intersected with existing edge ("
            << u << "," << v <<")\n"
            << "Incoming Line: ("<< A0 <<"," << B0 <<") ---> ("<< A1 <<"," << B1 <<") intersects with "
            << "Graph Line: ("<< A2 <<"," << B2 <<") ---> ("<< A3 <<"," << B3 <<") at point ("
            << ix << "," <<  iy << ")";

      DM_("Leaving: ") << __func__ << std::endl ;
      return 1;
    }

  } /* End of for edges */

DM_("Leaving: ") << __func__ << std::endl ;
return 0;
}



/*----------------------------------------------------------------------
 * Remove vertices by looping through the graph's vertices and 
 * checking property(*) satisfication
 *
 * Loops through the graph trying to remove vertices one by one.
 * Algo:
 * foreach vertex v-bar in original graph
 *  - Get the adjacent nodes
 *  - If > 3 --> Triangulate the points by pushing to a vector, keep mapping to
 *               original vertex IDs though. Test the triangulation/graph for Property(*)
 *               If successful -> Clear vertex and add the other edges in 
 *               triangulation to g
 *  - Else <= 3 --> Check for property(*) and clear or leave vertex.
 ----------------------------------------------------------------------*/
bool remove_vertices( mygraph_t &g, mygraph_t &g_orig,
                      boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type& id,
                      WeightMap &weightmap)
{
  DM_("Entering: ") << __func__ << std::endl ;

  vertex_iterator i,end;
  vert v_bar;


  /* vertex to remove */
  for(boost::tie(i,end) = boost::vertices(g); i != end; ++i) {                            /* Iterate through vertices in graph */

    v_bar = id[*i];
    int degree = boost::out_degree( v_bar, g);
    bool one =0, two =0;
    ADDBACK_EDGES.clear();                                                                /* Clear up the global caches */

    D_( "\n\n ********************************************************************* ");
    D_( " ****************************************** Processing Vertex: ")
        << v_bar << "(" << degree << ")" << std::endl;
    D_( "\n\n ********************************************************************* ");

    if ( boundary_points_map.find(v_bar) == boundary_points_map.end() ) {/* not found */
    } else { /* found */
        D_("Skipping Boundary Vertex");
        continue;
    }

    #if DEBUG_FLAG // && DEBUG_MORE
      /**********************************************************************************/
      /* Check a DLS of depth = 2 */
      /* Lets do a depth limited exploration */
      std::vector<vert> before_DLS_vertices;
      std::vector<std::pair<vert,vert>> before_DLS_edges;
      DLS( g, v_bar, -1, 2, before_DLS_vertices, before_DLS_edges);
      dump_graph_visual( "before_error_graph", before_DLS_vertices, before_DLS_edges );
      before_DLS_vertices.clear();
      before_DLS_edges.clear();
      /**********************************************************************************/
    #endif

    /***********************************************************************************
      * CHECK PROPERTY 1 (depends on triangulation due to net_dist on new graph )
      *********************************************************************************/
    if( degree >= 3 ){  /* Time to triangulate and check property(*)  */
      mygraph_t small_g( degree );

      D_("--- Triangulating this --- ");
      triangulate(g, small_g, v_bar, degree, id);                                         /* Triangulate */

      #if DEBUG_FLAG && DEBUG_MORE
        D_( "--- Small Graph looks like this... --- ");
        dump_graph( small_g, id, weightmap);
      #endif

      /* Note: We must pass in the original (unchanged graph) and 
       * the small graph (gap triangulation) 
       */
      one = check_property_1( g_orig, small_g, v_bar, 1 );
      D_( "-- [tri=1] Was Property 1 satisfied? ==> ") << one ;
    }
    else{
        /* Note: We must pass in the original (unchanged graph) and the 
         * changing graph (remaining vertices) 
         */
      D_("--- NOT  Triangulating this --- ");
      one = check_property_1( g_orig, g, v_bar, 0 );
      D_("-- [tri=0] Was Property 1 satisfied? ==> ") << one <<"\n\n";
    }

    /***********************************************************************************
      * CHECK PROPERTY 2 (independant of triangulation )
      *********************************************************************************/
    DELTAS_TO_MOVE dmove;
    two = check_property_2( g_orig, g, v_bar, dmove );
    D_( "------ Was Property 2 satisfied? ==> ") << two <<"\n";

    if(two){
      DELTAS_TO_MOVE::iterator dpos;
      for( dpos = dmove.begin(); dpos != dmove.end(); ++dpos){
        DM_( "Updating L(v_bar) with (") << (*dpos).v << "," << (*dpos).u << "," << (*dpos).delta << ") \n";
        /* Updating L(v) */
        update_delta_in_L((*dpos).v, (*dpos).u, (*dpos).delta );                          /* We passed in (v,u,delta_prime) */
      }
    }
    dmove.clear(); 

    /***********************************************************************************
     * CHECK PROPERTY (*) to finally remove the vertex and introduce new edges 
     * if needed to our reduced graph g (G')
     *********************************************************************************/

    if(one && two && vertices_left(g_orig) > 2){                                          /* We cannot remove the very LAST vertex in g' */

        remove_single_vertex( v_bar, g, g_orig );                                           /* Remove single vertex */

      /* ADD THE NEW EDGES HERE NOW BECAUSE THE EDGES ARE NOW VALID FOR ADDITION ! */
      std::vector<ed>::iterator pos;
      int stop=0;
      for( pos = ADDBACK_EDGES.begin(); pos != ADDBACK_EDGES.end(); ++pos){

          if( !boost::lookup_edge( pos->first, pos->second, g).second &&
              !boost::lookup_edge( pos->second, pos->first, g).second
              ){                  /* This edge isn't there, so add it */

              /* Check for Intersection : If intersection --> Don't add this new edge! */
              if( !doesIntersect(g, id, boost::vertex(pos->first, g),
                                        boost::vertex(pos->second,g)) )
              {

                doesIntersect(g, id,  boost::vertex(pos->first, g),
                                      boost::vertex(pos->second,g));

                D_( "Adding Edge (+) = (" ) << pos->first
                      << " ---> "<< pos->second <<")";

                edge_descriptor e; bool inserted;
                boost::tie(e, inserted) =
                boost::add_edge(  boost::vertex(pos->first,  g),
                                  boost::vertex(pos->second, g), g);                        /* Add edge to graph */
                weightmap[e] = get_euclidean(boost::vertex(pos->first, g),                  /* Add edge weight */
                                           boost::vertex(pos->second, g));

             }
             else {
                D_( "******* Edge : (") << pos->first
                              << " ---> "<< pos->second
                              << ") needs to be discarded, intersects with EXISTING edge!";
                stop = 1;
             }


          } /* End of if lookup is ok */
          else
            D_( "******* Edge : (") << pos->first << " ---> "<< pos->second <<") already exists!\n";

      } /* End of iterating through edges */

      //if(stop){
      //    exit(-1);
      //}

      REMOVED_VERTICES ++;
    } /* End of if */
    ADDBACK_EDGES.clear();
    dmove.clear();
 } /* End of foreach vertex */
  DM_("Leaving: ") << __func__ << std::endl ;
  return 1;
}

/*----------------------------------------------------------------------
 * remove a single vertex operation
 * Algo:
 * for vertex v-bar to be removed
 *  - Get the adjacent nodes
 *  - foreach v \in adj nodes
 *      - Check and update delta in L(v) for v_bar
 *  - end foreach
 *  - remove the v_bar from graph g and add v_bar to removed_list
 ----------------------------------------------------------------------*/
bool remove_single_vertex( vert v_bar, mygraph_t &g, mygraph_t &g_orig)
{
  DM_("Entering: ") << __func__ << std::endl ;

  adj_iterator ai, aend;
  /* Iterate through adjacent vertices of v_bar */
  for(boost::tie(ai,aend) = boost::adjacent_vertices(v_bar, g);
          ai != aend; ++ai)
  {
    //double delta = get_delta(v_bar, *ai, g_orig);                                         /* Compute Delta on Original Graph */
    //std::cout << "In remove_single_vertex: Delta = "<< delta << std::endl;
    //update_delta_in_L(  *ai, v_bar, delta );                                              /* Update Delta in L(v) */
    /* Add to HST(v_bar) adjacent vertices */
    hostMap::iterator h_iter;
    h_iter = HST.find(v_bar);                                                             /* key = v_bar */
    if(h_iter != HST.end()){                                                              /* Found the key !! */
        h_iter->second.push_back(*ai);
    }
    else
      HST.insert( std::make_pair(v_bar, std::vector<vert>{*ai}  ) );
  } /* End of for adj vertices (v_bar) */
  REMOVED.insert(std::make_pair(v_bar,1));                                                /* Add to Removed vertices */
  boost::clear_vertex(v_bar, g);                                                          /* Remove all edges to this vertex */

  DM_("Leaving: ") << __func__ << std::endl ;

return 1;
}

/*----------------------------------------------------------------------
 * Returns the estimated distance on a graph
 ----------------------------------------------------------------------*/
double estimated_distance(vert v_bar, vert u)
{
  DM_("Entering: ") << __func__ << std::endl ;

  mainMap::iterator l_iter;
  double delta=0.0;
  
  l_iter = L.find(v_bar);                                                               /* Should find inMap */
  if( l_iter != L.end()){                                                               /* Found the key !! */
    
    /* Iterate through the (v,delta) pairs */
    inMap::iterator pos;
    pos = (l_iter->second).find(u);
    if( pos != (l_iter->second).end()){ /* We found our delta */
      delta = pos->second;
    }
  }
  DM_("Leaving: ") << __func__ << std::endl ;
  return delta + get_euclidean(v_bar,u);
}

/*----------------------------------------------------------------------
 * Checks Property(2) for v_bar in graph g.
 * Foreach v in adjacent(v_bar,g)
 * {
 *   Foreach w in L(v_bar)
 *   {
 *      compute delta_prime and check property
 *   }
 *  }
 ----------------------------------------------------------------------*/
bool check_property_2(mygraph_t &g_orig, mygraph_t &g, vert v_bar, DELTAS_TO_MOVE &dmove )
{
  DM_("Entering: ") << __func__ << std::endl ;

  adj_iterator ai, aend;
  for(boost::tie(ai,aend) = boost::adjacent_vertices(v_bar, g); ai != aend; ++ai){
    vert v = *ai;
    DM_("\n H: (v,v_bar) = ") << v <<","<<v_bar << std::endl;
    mainMap::iterator l_iter;
    l_iter = L.find(v_bar);                                                               /* Should find inMap */
    if( l_iter != L.end()){                                                               /* Found the key !! */

      //std::cout << "L("<< v_bar <<") entry was found! \n";
      /* Iterate through the (v,delta) pairs */
      inMap::iterator pos;
      for(pos= (l_iter->second).begin(); pos!= (l_iter->second).end(); ++pos){
          vert u              = pos->first;

            double delta        = pos->second;                                              /* delta = 2nd part(u,delta) in L(v-bar) */
            double uv_bar       = get_euclidean(u,v_bar);
            double vv_bar       = get_euclidean(v,v_bar);
            double uv           = get_euclidean(u,v);
            double delta_prime  = uv_bar + vv_bar - uv + delta;                             /* new delta_prime computed */
        
            dmove.push_back({v,u,delta_prime});
            double mid          = fabs(delta_prime + uv);
            double d = 0.0;
            net_distance( g_orig, d, u,v );                                                 /* net dist (u,v) on original graph */

            DM_( "(v,u,delta,uv_bar,vv_bar,uv,delta_prime,net(u,v)) = (")
             << v << "," << u << "," << delta << ","
             << uv_bar << "," << vv_bar << ","
             << uv << "," << delta_prime << "," << d <<  ")\n";
            double LB = d;                                                                  /* Lower Bound */
            double UB = beta * d;                                                           /* Upper Bound */

            DM_( "\n[prop 2 fn] (lb,mid,ub) = (")
              << std::fixed << std::setprecision (4)
              << LB << ", " << mid <<", " << UB << ")\n";
        
            if( !(( (abs(LB - mid) <= 0.001) || (LB < mid) ) &&
                ( (abs(UB - mid) <= 0.001) || (UB > mid) )))
            {
                dmove.clear();
                DM_("Leaving: ") << __func__ << std::endl ;
                return 0;
            }
      } /* End of for */
    }
    else{ 
      DM_( "L("<< v_bar <<") entry was NOT found, so empty set satisfies (2)! \n");
      double vv_bar       = get_euclidean(v,v_bar);
      dmove.push_back({v,v_bar, vv_bar});
      DM_("Leaving: ") << __func__ << std::endl ;
      return 1;
    }
  } /* End of v loop */

DM_("Leaving: ") << __func__ << std::endl ;
return 1;
}

/*----------------------------------------------------------------------
 * Checks Property(1) for v_bar in graph g.
 ----------------------------------------------------------------------*/
bool check_property_1(mygraph_t &big_g, mygraph_t &tiny_g, vert v_bar, bool triangulated )
{
  DM_("Entering: ") << __func__ << std::endl ;

  boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type                  /* Get vertex property map */
      tid = get(&VertexProperties::index, tiny_g);
  std::vector<vert> vlist;

  /* If "triangulated" is true --> the small graph rep'n of a triangulation
   * is being passed through.
   */
  if( triangulated ){
      boost::graph_traits<mygraph_t>::vertex_iterator i, end;
      for(boost::tie(i,end) = boost::vertices(tiny_g); i != end; ++i)
          vlist.push_back(*i);  /* Add vertex to list */
   }
  else {  /* v_bar with degree < 4 was passed in, so DIRECT check of property 1! */
    adj_iterator ai, aend;
    for(boost::tie(ai,aend) = boost::adjacent_vertices(v_bar, tiny_g); 
        ai != aend; ++ai)                                                                 /* Iterate through adjacent vertices of v_bar */
      vlist.push_back( *ai );
  }
  /* Vertex Pair Generation */
  for(std::size_t i =0; i< vlist.size(); i++){
      for( std::size_t j=i+1; j< vlist.size(); j++){
          std::size_t s =  vlist[i] ;                                                     /* start vertex */
          std::size_t t =  vlist[j] ;                                                     /* end vertex */
          double d_small = 0.0, d_big = 0.0;

          net_distance( tiny_g, d_small, s,t );                                           /* Compute A* distance(s,t) on reduced triangulation graph */
          net_distance( big_g,  d_big,   tid[s], tid[t] );                                /* Compute A* distance(s,t) on large original graph */

          DM_( "A* on small:(" ) << s <<"," << t << ") --> "<< d_small <<"\n";
          DM_( "A* on big:  (" ) << tid[s] <<"," << tid[t] << ") --> "<< d_big   <<"\n";
          double LB = d_big, UB = beta*beta*d_big;                                                  /* Bounds */

          /* IMPORTANT! : Add edges to be transported back to the 
           * main graph when no triangulation 
           */
          if(!triangulated)
            ADDBACK_EDGES.push_back( std::make_pair( tid[s], tid[t] ));

          DM_( "[prop 1] -------- >> ( lb, small_distance, ub) = (" )
               << std::fixed << std::setprecision (4)
               << LB << ", " << d_small <<", " << UB << ")\n";

          if( !(( (abs(LB - d_small) <= 0.001) || (LB < d_small) ) &&
              ( (abs(UB - d_small) <= 0.001) || (UB > d_small) )) )
          {
            ADDBACK_EDGES.clear();
            DM_("Leaving: ") << __func__ << std::endl ;
            return 0;
          }
       }
  }
  vlist.clear();                                                                          /* Remove all the vertices */
DM_("Leaving: ") << __func__ << std::endl ;
return 1;
}


/*----------------------------------------------------------------------
 * Triangulate the vertices and generate a smaller graph to 
 * test property(*) on
 ----------------------------------------------------------------------*/
bool triangulate(mygraph_t &g, mygraph_t &small_g, vert v_bar, int degree,
                 boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type& id)
{
  DM_("Entering: ") << __func__ << std::endl ;

  /* Local Map to track mapping in small_g to g */
#ifndef WIN32
  vert vertex_mapping[degree + 1];
#else
  vert* vertex_mapping = new vert[degree + 1];
#endif
  std::vector< std::pair<DelPoint,unsigned> > points;

  /* For the smaller graph of triangulation */
  boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type                   /* small: Get vertex property map */
      sid = get(&VertexProperties::index, small_g);
  WeightMap sweightmap = get(boost::edge_weight, small_g);                                /* small: Get edge weights */

  adj_iterator ai, aend;
  int v_small = 0;
  for(boost::tie(ai,aend) = boost::adjacent_vertices(v_bar, g);
          ai != aend; ++ai)                                                               /* Iterate through adjacent vertices of v_bar */
  {
    vert adj_vertex = id[*ai];
    // std::cout << "*** Adjacent Vertex: "<< adj_vertex << std::endl;
    vertex_mapping[ v_small ] = adj_vertex;                                               /* Store the mapping of vertex to orig graph */

    //std::cout << "*** Vertex_mapping ["<< v_small << "] = "<< adj_vertex << std::endl;

    double x = vertex_props[adj_vertex].x;
    double y = vertex_props[adj_vertex].y;
    points.push_back( std::make_pair(DelPoint(x,y),v_small) );                            /* Make points vector for triangulation */
    ++ v_small;
  } /* End of adj loop */

  Delaunay T;                                                                           /* Triangulate !! */
  T.insert( points.begin(),points.end() );
  DM_( "# Of vertices in T = ") << T.number_of_vertices() << " degree = "
             << degree << " points size: "<< points.size() << std::endl;
  CGAL_assertion( T.number_of_vertices() == degree );

  Delaunay::Finite_vertices_iterator vit;                                                 /* Iterate through vertices of T to make graph small_g */
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {
    if( points[ vit->info() ].first != vit->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
    vert i = points[ vit->info() ].second;
    sid[i] = vertex_mapping[i];                                                           /* Ensure to set vertex descriptor to the vertex_mapping[i];
                                                                                          * original graph's vertex ID */
   }	/* End of for */

  Delaunay::Finite_edges_iterator eit;
  for (eit = T.finite_edges_begin(); eit != T.finite_edges_end(); ++eit)
  {
    Delaunay::Edge e = *eit;
    int i1 = e.first->vertex( (e.second+1)%3 )->info();
    int i2 = e.first->vertex( (e.second+2)%3 )->info();
    edge_descriptor ed; bool inserted;
    vert a = boost::vertex( i1, small_g);
    vert b = boost::vertex( i2, small_g);
    boost::tie(ed, inserted) = boost::add_edge(a,b,small_g);

    /* Add edges for this triangulation to this */
    ADDBACK_EDGES.push_back( std::make_pair(  vertex_mapping[i1],
                                              vertex_mapping[i2] ));
    double w=0;
    net_distance( g, w, vertex_mapping[i1], vertex_mapping[i2] );                         /* Compute the network distance */

    //double w =  sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));             /* Compute euclidean 3d distance */
    sweightmap[ed] = beta * w;                                                            /* Add as a weight */
  }
  DM_("Leaving: ") << __func__ << std::endl ;

#if defined WIN32
  delete[] vertex_mapping;
#endif

  return 1;
}


/*----------------------------------------------------------------------
 * Reads in the graph from the file passed in. The graph input file must 
 * have the format:
 * 1: g <V=num of vertices> <E=num of edges>
 * V lines: <vertex ID> <x> <y> <z> // One line per vertex
 * E lines: <src vID> <tgt vID> <weight:3D Euclidean distance>
 * ----- Example Graph File -------------------------------------
 * g 5 6
 * 0 10.0 5.0 2.0
 * 1 15.0 10.0 7.0
 * 2 20.0 5.0 8.0
 * 3 15.0 0.0 18.0
 * 4 35.0 2.0 10.0
 * 0 1 8.66
 * 0 2 11.66
 * 1 2 7.14
 * 2 3 12.24
 * 2 4 15.42
 * 3 4 21.63
 ----------------------------------------------------------------------*/
bool readgraph(std::istream & in, mygraph_t &g, std::size_t n_vertices,
               std::size_t n_edges,
               boost::property_map
               < mygraph_t, std::size_t VertexProperties::*>::type& id)
{
  DM_("Entering: ") << __func__ << std::endl ;

  for (size_t i = 0; i < n_vertices; ++i) {                                               /* Read in vertices */
    std::size_t vertex_id;
    double x, y, z;
    if (!(in >> vertex_id >> x >> y >> z))
      throw std::runtime_error("Cannot read vertices!");
      DM_( "id ") << vertex_id << "x: "<< x
                << " y: "<< y << " z: "<< z << std::endl;
      vertex_props.push_back( Point3{ x,y,z } );                                          /* Add to external vertex property map */
      id[i] = i;
  } /* End of for vertices */

  WeightMap weightmap = get(boost::edge_weight, g);
  for (size_t i = 0; i < n_edges; ++i) {                                                  /* Read in edges */
    std::size_t a, b;
    double w;
    if (!(in >> a >> b >> w))
      throw std::runtime_error("Cannot read edges!");
      edge_descriptor e; bool inserted;
      boost::tie(e, inserted) = boost::add_edge( boost::vertex(a, g),
                                                 boost::vertex(b, g), g);                 /* Add edge to graph */
      weightmap[e] = w;                                                                   /* Add edge weight */
  } /* End of for edges */

  DM_("Leaving: ") << __func__ << std::endl ;
  return 0;
}


/*----------------------------------------------------------------------
 * Depth Limited Search (DLS) to get a subgraph of the entire graph
 * starting out from a problem vertex
----------------------------------------------------------------------*/
void DLS( mygraph_t &G, vert start, vert goal, int depth_limit,
          std::vector<vert> &v, std::vector<std::pair<vert,vert>> &e )
{
  DM_("Entering: ") << __func__ << std::endl ;

  std::stack<vert> S;                   /* Stack to push vertices on */
  std::map<vert,int> depth_map;         /* Keeps track of the depth of a vertex */
  std::map<vert,int> visited_map;       /* Keeps track of whether this vertex was visited or not */

  S.push(start);
  depth_map[start]    = 0;      /* Root node at depth 0  */
  visited_map[start]  = 1;      /* Root node not visited */

  while( S.size() != 0){
    vert parent       = S.top();
    S.pop();
    int  parent_depth = depth_map[parent];
    v.push_back(parent);

    if(parent == goal)
      break;

    if(parent_depth == depth_limit)
      continue;

    else {
      // boost::graph_traits<mygraph_t>::vertex_iterator i, end;
      adj_iterator ai, aend;
      /* Loop through the adjacent vertices of parent */
      for(boost::tie(ai,aend) = boost::adjacent_vertices(parent, G);
                  ai != aend; ++ai)
      {
        /* Check if this vertex was already visited before */
        if( visited_map[*ai] != 1 ){
          int child_depth = parent_depth + 1;
          S.push(*ai);
          depth_map[ *ai ] = child_depth;
          visited_map[ *ai ] = 1;
          /* Add the edges too */
          e.push_back( std::make_pair(parent,*ai) );
        } /* End of if */
      } /* End of Adjacent vertices */
     }  /* End of else */


  } /* End of while */

  DM_("Leaving: ") << __func__ << std::endl ;
}

/*----------------------------------------------------------------------
 * A* to compute shortest network path
 ----------------------------------------------------------------------*/
 bool net_distance(mygraph_t &g, double & dist, std::size_t s, std::size_t t )
 {
   DM_("Entering: ") << __func__ << std::endl ;

   std::vector<vert> p(num_vertices(g));
   std::vector<cost> d(num_vertices(g));
   vert start = boost::vertex(s, g);                                                      /* Set start */
   vert goal  = boost::vertex(t, g);                                                      /* Set end */

   try {
     boost::astar_search                                                                  /* call A* named parameter interface */
       (g, start,
         distance_heuristic<mygraph_t, cost, Point3*>
          (&vertex_props[0], goal),
         boost::predecessor_map(&p[0]).distance_map(&d[0]).
         visitor(astar_goal_visitor<vert>(goal)));

   } catch(found_goal fg) {                                                               /* found a path to the goal ! */
     std::list<vert> shortest_path;
     for(vert v = goal;; v = p[v]) {
       shortest_path.push_front(v);
       if(p[v] == v)
          break;
      }

     DM_( "Shortest path from ") << start << " to "
           << goal << ": ";
      std::list<vert>::iterator spi = shortest_path.begin();
      DM_( start );
      for(++spi; spi != shortest_path.end(); ++spi)
        DM_( " -> ") << *spi;
      DM_( "Total Distance: ") << d[goal] << std::endl;

     dist = d[goal];
     DM_("Leaving: ") << __func__ << std::endl ;
     return 0;
    }
    std::cout << "Path was not found ! \n";
    DM_("Leaving: ") << __func__ << std::endl ;

    return -1;
 }

/*----------------------------------------------------------------------
 * Dump Graph
 ----------------------------------------------------------------------*/
void dump_graph(mygraph_t &G, boost::property_map< mygraph_t,
                std::size_t VertexProperties::*>::type &id,
                 WeightMap &weightmap)
{
   boost::graph_traits<mygraph_t>::vertex_iterator i, end;
   boost::graph_traits<mygraph_t>::out_edge_iterator ei, edge_end;

   std::cout << "********** MY GRAPH *******************\n";
   for(boost::tie(i,end) = boost::vertices(G); i != end; ++i) {
       std::cout << id[*i];
       std::cout << " #degree : " << boost::degree( *i, G) << std::endl;

       for (boost::tie(ei,edge_end) = boost::out_edges(*i, G); ei != edge_end; ++ei)
          std::cout << " --" << weightmap[*ei] << "--> " << id[target(*ei, G)] << "  ";
     std::cout << std::endl;
   } /* End of for */
 }


/*----------------------------------------------------------------------
 * Dump Graph to visualize
 ----------------------------------------------------------------------*/
void dump_graph_visual( std::string fname, std::vector<vert> &v,
                        std::vector<std::pair<vert,vert>> &e )
{
  std::ofstream g_file;
  g_file.open (fname.c_str());

  g_file << v.size() << " " << e.size() << "\n";

  /* Dump out Vertices */
  for( auto i : v){
     g_file << std::fixed
            << std::setprecision (2)
            << i << " "
            << vertex_props[i].x << " "
            << vertex_props[i].y << "\n";
  }

  /* Dump out Edges */
  for( auto i : e){
    g_file << std::fixed
           << std::setprecision (2)
           << i.first << " " << i.second << "\n";
  }
  g_file.close();
}

/*----------------------------------------------------------------------
 * Dump Graph to visualize
 ----------------------------------------------------------------------*/
void dump_graph_visual( std::string fname, mygraph_t &G,
                        boost::property_map< mygraph_t,
                        std::size_t VertexProperties::*>::type &id )
{
  std::ofstream g_file;
  boost::graph_traits<mygraph_t>::vertex_iterator i, end;
  boost::graph_traits<mygraph_t>::edge_iterator  ei, edge_end;

  g_file.open (fname.c_str());
  int numV = boost::num_vertices(G) - REMOVED_VERTICES;

  g_file << numV  << " "
         << boost::num_edges(G)     << "\n";

  /* Dump out all the vertices */
  for(boost::tie(i,end) = boost::vertices(G); i != end; ++i){
   int degree = boost::degree( *i, G);
   if(degree > 0 )
      g_file << id[*i]     << " "
      << vertex_props[ *i ].x << " "
      << vertex_props[ *i ].y << "\n";
  }

  /* Dump out all the edges */
  for (boost::tie(ei,edge_end) = boost::edges(G); ei != edge_end; ++ei)
    g_file << std::fixed
       << std::setprecision (2)
       << boost::get( id, boost::source(*ei,G))  << " "
       << boost::get( id, boost::target(*ei,G)) << "\n";

  g_file.close();
}

/*----------------------------------------------------------------------
 * Archive the Graph
 ----------------------------------------------------------------------*/
void archive_graph(mygraph_t &G, boost::property_map< mygraph_t,
                std::size_t VertexProperties::*>::type &id,
                WeightMap &weightmap, std::ofstream &of)
{
  double nv = boost::num_vertices(G);
  double ne = boost::num_edges(G);
  
  boost::graph_traits<mygraph_t>::vertex_iterator i, end;
  boost::graph_traits<mygraph_t>::edge_iterator ei, edge_end;

  of << "g "<< nv << " "<< ne << std::endl;
  /* Write out all the vertices */
  for(boost::tie(i,end) = boost::vertices(G); i != end; ++i) {
    of << id[*i]     << " "
    << vertex_props[ *i ].x << " "
    << vertex_props[ *i ].y << " "
    << vertex_props[ *i ].z << std::endl;
  }
  /* Write out all the edges */
  for (boost::tie(ei,edge_end) = boost::edges(G); ei != edge_end; ++ei)
    of << boost::get( id, boost::source(*ei,G)) << " "
    << boost::get( id, boost::target(*ei,G)) << " "
    << weightmap[*ei] << std::endl;     
}


/*----------------------------------------------------------------------
 * Update delta in L(v) : If same v-bar entry already exists for a v,
 *                            If existing delta < new delta
 *                              do nothing
 *                            else
 *                              update existing delta to new delta value.
 *                         else
 *                             Insert new (v-bar,delta)
 ----------------------------------------------------------------------*/

void update_delta_in_L( vert v, vert v_bar, double delta )
{
  mainMap::iterator v_iter;
  inMap::iterator   vbar_iter;

  /* Insert entry into map */
  v_iter = L.find(v);                                                                     /* look for the key i.e. v */
  if (v_iter != L.end()){                                                                 /* found the key, so insert */

      vbar_iter = v_iter->second.find(v_bar);
      if( vbar_iter != v_iter->second.end() ){                                            /* Found this v-bar, so check it's delta */
          double d_prime = vbar_iter->second;

          if( d_prime < delta ){/* Do Nothing */}
          else
            (v_iter->second)[ v_bar ] = delta;                                            /* Update delta value */
      }
      else                                                                                /* Did not find such a v-bar entry! */
        v_iter->second.insert( std::make_pair(v_bar, delta) );
   }
   else {                                                                                 /* Not found v,v-bar entry in main map so new entry */
      L.insert( std::make_pair( v, inMap() ));
      L[v].insert( std::make_pair(v_bar, delta  ));
   }  /* End of else */
} /* End of function */


/*----------------------------------------------------------------------
 * Dump out L(v)
 ----------------------------------------------------------------------*/
void dumpL( mainMap &LL)
{
  mainMap::iterator it;
  inMap::iterator   inner_it;

  std::cout << "********** DUMPING L ************\n";
   for ( it= LL.begin() ; it != LL.end(); it++ ) {
     std::cout << "\n ******* " << (*it).first << std::endl;
     for( inner_it=(*it).second.begin(); inner_it != (*it).second.end(); inner_it++)
       std::cout  << (*inner_it).first << " => " << std::setprecision(4)
                  << (*inner_it).second << std::endl;
   }
}

/*----------------------------------------------------------------------
 * Dump HST
 ----------------------------------------------------------------------*/
void dumpHST(hostMap &HH)
{
  hostMap::iterator pos;
  std::vector<vert>::iterator spos;

  std::cout << "********** DUMPING HST ************\n";
  for(pos= HH.begin(); pos!= HH.end(); ++pos){
    std::cout << pos->first << "|";
    /* Loop through the set values */
    for( spos= (pos->second).begin(); spos != (pos->second).end(); ++spos )
      std::cout << "(" <<  (*spos) << "), ";
    std::cout << "\n";
  }   /* End of for */
}

/*----------------------------------------------------------------------
 * Dump REMOVED vertices
 ----------------------------------------------------------------------*/
void dumpREMOVED(remMap &REM)
{
  remMap::iterator pos;
  std::cout << "********** DUMPING REMOVED VERTICES ************\n";
  for( pos= REM.begin(); pos != REM.end(); ++pos )
    std::cout << "(" << pos->first << "), ";
  std::cout << "\n";
}

/*----------------------------------------------------------------------
 * Get real number of vertices in graph (Subtract cleared vertices)
 ----------------------------------------------------------------------*/
int vertices_left(mygraph_t &g)
{
  int total = boost::num_vertices(g);
  int removed = REMOVED.size();
  return total-removed;
}


/*----------------------------------------------------------------------
 * Computes the Euclidean (i1,i2) - vertex IDs passed in
 ----------------------------------------------------------------------*/
double get_euclidean( vert i1, vert i2)
{
  double x2 = vertex_props[ i1 ].x;
  double y2 = vertex_props[ i1 ].y;
  double z2 = vertex_props[ i1 ].z;
  double x1 = vertex_props[ i2 ].x;
  double y1 = vertex_props[ i2 ].y;
  double z1 = vertex_props[ i2 ].z;
  
  return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));                      /* Compute euclidean 3d distance */
}

/*----------------------------------------------------------------------
 * Computes the Delta : Net_dist(s,t) - Euc_dist(s,t)
 ----------------------------------------------------------------------*/
double get_delta( vert i1, vert i2, mygraph_t &g)
{
  double net_d = 0.0;
  net_distance(g, net_d, i1, i2);
  double e = get_euclidean(i1,i2);
  double delta = fabs(net_d - e );
  
//  std::cout << "--------------------------- [Delta] NetDist ("
//            << i1 << "," << i2 << ") ==>  "
//            << net_d
//            <<" Euc = "<< e << "\n";
//  std::cout << "--------------------------- DELTA ("<< i1
//            << "," << i2 << ") = "<< std::setprecision(4)
//            << delta << std::endl;
  
  return ( delta < 0.001 ) ? 0 : delta;
}

/*----------------------------------------------------------------------
 * Display Statistics
 ----------------------------------------------------------------------*/
void get_stats()
{
  std::cout << "********************************************************** \n";
  std::cout << "Number of REMOVED vertices = " << REMOVED.size() << std::endl;
  hostMap::iterator pos;
  std::vector<vert>::iterator spos;

  int k=0;
  double vals=0;
  for(pos= HST.begin(); pos!= HST.end(); ++pos, ++k){
    /* Loop through the set values */
    for( spos= (pos->second).begin(); spos != (pos->second).end(); ++spos, ++vals )
    { }
  }   /* End of for */
  std::cout << "Avg Number of HOST values = " << vals/k << std::endl;
  std::cout << "********************************************************** \n";
}


bool exist(std::tuple<int, int> value, std::vector<std::tuple<int,int>> v)
{
    // return std::find(v.begin(), v.end(),value)!=v.end();
    for(auto it = v.begin(); it != v.end(); ++it) {
        if(std::get<0>(value) == std::get<0>(*it) && std::get<1>(value) == std::get<1>(*it)){
            return true;
        }
        if(std::get<0>(value) == std::get<1>(*it) && std::get<1>(value) == std::get<0>(*it)){
            return true;
        }
    }
    return false;
}

double distance(Point3 v1, Point3 v2){
    double d[] = {abs(v1.x-v2.x), abs(v1.y-v2.y), abs(v1.z-v2.z)};
    // if (d[0] < d[1]) std::swap(d[0],d[1]);
    // if (d[0] < d[2]) std::swap(d[0],d[2]);
    // double distance = d[0] * sqrt(1.0 + d[1]/d[0] + d[2]/d[0]);
    return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
}

// we assume there is not comment in OFF file
bool generate_graph(std::string in_path, std::string out_path){
    // std::ifstream infile("/etc/terrain_toolkit/models/off/small_terrain.off");
    // std::ofstream outfile("/etc/terrain_toolkit/models/graph/small_terrain.graph");

  std::ifstream infile(in_path);
  std::ofstream outfile(out_path);

  std::vector<std::tuple<int,int>> edge_list;

  if (infile.is_open() && outfile.is_open()) {
      std::string line;
      int num_v, num_f, num_e;
      double x,y,z;
      int k,v1,v2,v3;
      infile >> line;
      if(line.compare("OFF")!=0){
          std::cout<<"not off file!"<<std::endl;
      }
      infile >> num_v >> num_f >> num_e;
      outfile << "g "<< num_v << " "<< num_e << std::endl;
      // read vertex
      for(int i=0; i<num_v; i++){
          infile >> x >> y >> z;
          vertex_props.push_back( Point3{ x,y,z } );
          outfile << i << " " << x << " " << y << " " << z << std::endl;
      }
      // read face
      for(int i=0;i<num_f;i++){
          infile >> k >> v1 >> v2 >> v3;

          if(!exist(std::make_tuple(v1,v2), edge_list)){
              edge_list.push_back(std::make_tuple(v1,v2));
              outfile << v1 << " " << v2 << " " << distance(vertex_props[v1],vertex_props[v2]) << std::endl;
          }

          if(!exist(std::make_tuple(v1,v3), edge_list)){
              edge_list.push_back(std::make_tuple(v1,v3));
              outfile << v1 << " " << v3 << " " << distance(vertex_props[v1],vertex_props[v3]) << std::endl;
          }

          if(!exist(std::make_tuple(v2,v3), edge_list)){
              edge_list.push_back(std::make_tuple(v2,v3));
              outfile << v2 << " " << v3 << " " << distance(vertex_props[v2],vertex_props[v3]) << std::endl;
          }
      }
      infile.close();
      outfile.close();
      return true;
  }
  else{
      std::cout<<"file to open file" <<std::endl;
      return false;
  }
}

bool generateOff(std::string graph_path, float beta_in, std::string off_path)
{
  // init global var
  beta = beta_in;
  REMOVED_VERTICES = 0;
  
  std::ifstream is(graph_path);                                                          /* Read in the graph from file */
  char x;
  if(!(is >> x))
    throw std::runtime_error("Did not read first char in graph file");
  if( x != 'g')
    throw std::runtime_error("Must be g...found something else");

  std::size_t n_vertices, n_edges;
  if (!(is >> n_vertices >> n_edges))                                                     /* Read in #Vertices #Edges in graph */
    throw std::runtime_error("No #V and #E!");

  D_( "#V ") << n_vertices << " #E "<< n_edges << std::endl;

  mygraph_t g(n_vertices);                                                                /* Initialize graph knowing #Vertices! */

  boost::property_map< mygraph_t, std::size_t VertexProperties::*>::type                  /* Get vertex property map */
      id = get(&VertexProperties::index, g);
  WeightMap weightmap = get(boost::edge_weight, g);                                       /* Get edge weights */

  readgraph(is, g, n_vertices, n_edges, id );                                             /* Read in rest of graph */

  vertex_iterator vi, viend;
  int vnum=0;
  for (boost::tie(vi,viend) = vertices(g); vi != viend; ++vi)                             /* Make IDs in vertex properties */
      id[*vi] = vnum++;
  mygraph_t g_orig = g;                                                                   /* Keep a copy of the original graph */

  print_hull(P, ch2d(P, read_points( vertex_props )));

  int K = hull_vertices.size();
   /* Loop through convex hull edges and process them */

  for(int k=0; k< K; k++){
    int a = hull_vertices[k];
    int b = hull_vertices[(k+1) % K];
    //std::cout << a  << " -- " << b << " \n";

    double lx1 = vertex_props[a].x;  /* Line start x-coord */
    double ly1 = vertex_props[a].y;  /* Line start y-coord */

    double lx2 = vertex_props[b].x;  /* Line end x-coord */
    double ly2 = vertex_props[b].y;  /* Line end y-coord */

    int X = a, oldX = a;
    boundary_points_map.insert ( std::pair<int,int>(a,1) );

    while(X != b){

      double mindist = DBL_MAX;
      int    minv    = -1;

      /* Go through each neighbor */
      adj_iterator ai, aend;
      vert v = X;
      /* Iterate through adjacent vertices of a */
      for(boost::tie(ai,aend) = boost::adjacent_vertices(v, g);
              ai != aend; ++ai)
      {
        vert V = *ai;
        if(V == oldX) continue; /* Skip looking at past vertex */

        double px = vertex_props[V].x , py = vertex_props[V].y;
        double dist = FindDistanceToSegment(lx1, ly1, lx2, ly2, px, py);
        /* Find the minimum distance to edge */
        if(dist < mindist){
          mindist = dist;
          minv    = V;
        }
      } /* End of adjacency list search */
      oldX = X;
      X = minv;
      if(X!=b)
        boundary_points_map.insert ( std::pair<int,int>(X,1) );
    } /* End of while X != b */
  } /* End of k convex hull edges */

  INIT_TIMER

  START_TIMER
  D_( "Removing Vertices ... " );
  remove_vertices( g, g_orig, id, weightmap );                                            /* Remove vertices */
  STOP_TIMER("Simplification: Vertex Removal Done!")

  graph_to_off(g, id, off_path);                                                              /* Get OFF from the graph */

  /* MKA: Added before Bound computation */
  get_stats();
 return true;
}




}
#endif // CRITICAL_H

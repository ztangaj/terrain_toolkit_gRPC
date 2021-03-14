/**********************************************************************************
** Filename       : critical.cpp
** Authors        : Manu Kaul
** Last Modified  : 03 Oct 2012
**
** Description	  : Implement the removal of vertices to arrive at critical vertices
**                  in final graph G' and output the results to a flat file
**
**********************************************************************************/

// cross platflam defining
#ifndef defined ( WIN32 )
#define __func__ __FUNCTION__
#endif
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "critical.h"
using namespace simplify;
/*----------------------------------------------------------------------
 Main Function
 ----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  std::string inf, outf;

  DM_("Entering: ") << __func__ << ", argc=" << argc ;

  // Declare the supported options.
  boost::program_options::options_description d("Allowed options for critical.cpp");
  d.add_options()
      ("help","produce this help message")
      ("in",  boost::program_options::value<std::string>(), "load graph from file 'arg'")
      ("out", boost::program_options::value<std::string>(), "output off to 'arg'")
      ("beta",boost::program_options::value<double>(),      "beta = 'arg'")  
      ;

  boost::program_options::variables_map m;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, d), m);
  boost::program_options::notify(m);

  if (m.count("help")){
    //Display the options_description
    std::cout << d << "\n";
    std::cout << "./critical --in G_bounds.dat --out output.txt --beta 2.0\n";
  }

  if(m.count("in") ){
    inf = m["in"].as<std::string>();
    D_("input graph : ") << inf;
  }
  if(m.count("out") ){
    outf = m["out"].as<std::string>();
    D_("output off: ")<< outf;
  }
  if(m.count("beta") ){
    beta = m["beta"].as<double>() ;
    D_( "beta : ")<< beta ;
  }
  else{
    std::cout << "No args set\n";
  }
  std::ifstream is(inf.c_str());                                                          /* Read in the graph from file */


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


#if DEBUG_FLAG && DEBUG_MORE
  D_("---- ORIGINAL GRAPH -----");
  dump_graph( g_orig, id, weightmap);
#endif
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

   #if DEBUG_FLAG && DEBUG_MORE
    DM_(" ---------------- BOUNDARY VERTICES ----------------");
    for(const auto& i : boundary_points_map){
       DM_( i.first ) << " ";
       DM_("\n");
    DM_(" ---------------- BOUNDARY VERTICES [DONE] ---------");
   }
   #endif

   INIT_TIMER

   START_TIMER
   D_( "Removing Vertices ... " );
   remove_vertices( g, g_orig, id, weightmap );                                            /* Remove vertices */
   STOP_TIMER("Simplification: Vertex Removal Done!")


    #if DEBUG_FLAG && DEBUG_MORE
      dump_graph( g, id, weightmap);
      dump_graph_visual( "all_graph_now", g, id );
      dump_graph( g_orig, id, weightmap);
      dumpL(L);
      dumpHST(HST);
      dumpREMOVED(REMOVED);
    #endif

   graph_to_off(g, id, outf);                                                              /* Get OFF from the graph */
  
   /* MKA: Added before Bound computation */
   get_stats();

 return 0;
}
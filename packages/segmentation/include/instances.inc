#include "graph.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations: <captype, tcaptype, flowtype>
// IMPORTANT: 
//    flowtype should be 'larger' than tcaptype 
//    tcaptype should be 'larger' than captype

#ifdef _WIN32
template class Graph<int,int,__int64>;
#else
template class Graph<int,int,long long>;
#endif


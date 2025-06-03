//
//  AdjacencyList.h
//  RichCluster
//
//  Created by Sarah on 6/2/25.
//

#ifndef AdjacencyList_h
#define AdjacencyList_h

#include <unordered_map>
#include <unordered_set>

class AdjacencyList {
public:
  AdjacencyList(int nterms) : n_terms(nterms){};
  
  void addNeighbor(int node, int neighbor);
  bool hasNeighbor(int node, int neighbor);
  
  // return reference to adjList (used to initialize ClusterList)
  const std::unordered_map<int, std::unordered_set<int>>& getAdjList() const {
    return adjList;
  }
  
private:
  std::unordered_map<int, std::unordered_set<int>> adjList;
  int n_terms; // used to initialize size
};

#endif /* AdjacencyList_h */

//
//  LinkageMethod.cpp
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#include <stdio.h>
#include "LinkageMethod.h"

using Cluster = LinkageMethod::Cluster;

double LinkageMethod::computeLinkage(const Cluster& cluster1, const Cluster& cluster2) {
  if (method=="single")
    return LinkageMethod::single(cluster1, cluster2);
  else if (method=="complete")
    return LinkageMethod::complete(cluster1, cluster2);
  else if (method=="average")
    return LinkageMethod::average(cluster1, cluster2);
  else if (method=="ward")
    return LinkageMethod::ward(cluster1, cluster2);
  else if (method=="david")
    return LinkageMethod::david(cluster1, cluster2);
  else
    return 0;
}

double LinkageMethod::single(const Cluster& cluster1, const Cluster& cluster2) {
  double minDist = 100; // a default stupid value
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      if (i==j) 
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distFct(t1, t2);
      if (dist < minDist)
        minDist = dist;
    }
  }
  return minDist;
}

double LinkageMethod::complete(const Cluster& cluster1, const Cluster& cluster2) {
  double maxDist = 0; // a default stupid value
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      if (i==j)
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distFct(t1, t2);
      if (dist > maxDist)
        maxDist = dist;
    } 
  }
  return maxDist;
}

double LinkageMethod::average(const Cluster& cluster1, const Cluster& cluster2) {
  double totalDist = 0;
  int n_terms = 0;
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      if (i==j)
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distFct(t1, t2);
      totalDist += dist;
      n_terms ++;
    }
  }
  return totalDist/n_terms;
} 


double LinkageMethod::ward(const Cluster& cluster1, const Cluster& cluster2) {
  double nA = cluster1.size();
  double nB = cluster2.size();

  if (nA == 0 || nB == 0) {
    return 0.0;
  }

  double delta_AA = 0.0;
  for (int i : cluster1) {
    for (int j : cluster1) {
      if (i == j) continue;
      delta_AA += distFct(i, j);
    }
  }
  delta_AA /= (nA * nA);

  double delta_BB = 0.0;
  for (int i : cluster2) {
    for (int j : cluster2) {
      if (i == j) continue;
      delta_BB += distFct(i, j);
    }
  }
  delta_BB /= (nB * nB);

  double delta_AB = 0.0;
  for (int i : cluster1) {
    for (int j : cluster2) {
      delta_AB += distFct(i, j);
    }
  }
  delta_AB /= (nA * nB);

  return (nA * nB) / (nA + nB) * (2 * delta_AB - delta_AA - delta_BB);
}

double LinkageMethod::david(const Cluster& cluster1, const Cluster& cluster2) {
    int linkCount = 0;
    for (int i : cluster1) {
        for (int j : cluster2) {
            if (distFct(i, j) < cutoff) {
                linkCount++;
            }
        }
    }
    return static_cast<double>(linkCount);
}

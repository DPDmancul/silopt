// Copyright 2022 Davide Peressoni
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package silopt;

import silopt.Clustering.Clustering;
import silopt.Clustering.Cluster;
import silopt.Clustering.Point;

import java.util.List;
import javafx.util.Pair;

import org.apache.commons.math3.random.MersenneTwister;

public class KMeansPP<P extends Point<P>> extends Optimizable<P> {
  private final int k;
  private final int n;
  private final Cluster centers = new Cluster();
  private final Cluster not_centers = new Cluster();

  private MersenneTwister rand;

  KMeansPP(List<P> points, int k, int seed) {
    super(new Clustering<P>());
    this.k = k;
    this.rand = new MersenneTwister(seed);
    clustering.points = points;
    n = points.size();
    for (int i = 0; i < n; ++i)
      not_centers.add(i);
  }

  public Pair<Double, Integer> optimize() {
    int c = rand.nextInt(n);
    centers.add(c);
    not_centers.remove(c);
    clustering.clusters.add(new Cluster(c));

    for (int i = 1; i < k ; ++i) {

      // Compute probabilities
      double[] p = new double[n];
      double total = not_centers
        .stream()
        .mapToDouble(q -> Math.pow(clustering.d(q, centers), 2))
        .sum();
      for (int j : not_centers)
        p[j] = Math.pow(clustering.d(j, centers), 2) / total;

      // Sample new center
      double sample = rand.nextDouble();
      double sum_p = 0;

      for (int j : not_centers) {
        sum_p += p[j];
        if (sample <= sum_p) {
          centers.add(j);
          not_centers.remove(j);
          clustering.clusters.add(new Cluster(j));
          break;
        }
      }

    }

    // Find the center which minimizes the distance
    int[] index = new int[clustering.points.size()];

    for (int e : not_centers) {
      int argmin = -1;
      double min = Double.MAX_VALUE;

      for (int i = 0; i < k; ++i) {
        // The cluster contains only its center
        final Cluster c_i = clustering.clusters.get(i);
        final double d = clustering.d(e, c_i);
        if (d < min){
          min = d;
          argmin = i;
        }
      }
      index[e] = argmin;
    }

    // Assign not chosen points to nearest cluster
    for (int e : not_centers)
      clustering.clusters.get(index[e]).add(e);

    return null;
  }

  public Cluster get_centers() {
    return centers;
  }
}

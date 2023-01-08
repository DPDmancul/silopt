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

package silopt.Clustering;

import silopt.Sil;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

public class Clustering<P extends Point<P>> {
  /** Vector of clusters */
  public List<Cluster> clusters = new ArrayList<>();
  /** Vector of all points */
  public List<P> points = new ArrayList<>();

  public Clustering() {}

  /**
   * Create a new clustering
   *
   * @param v an iterable of iterables of points. Each inner iterable represents a cluster.
   */
  public <J extends Iterable<P>, I extends Iterable<J>> Clustering(final I v) {
    for (final Iterable<P> c : v) {
      Cluster cluster = new Cluster();
      for (final P p : c) {
          cluster.add(points.size());
          points.add(p);
      }
      clusters.add(cluster);
    }
  }

  public <I extends Iterable<P>> Clustering(final I... v) {
    this(Arrays.asList(v));
  }

  public Clustering(P[]... v) {
    this(Arrays.asList(v).stream().map(Arrays::asList).toList());
  }

  /**
   * Computes distance between two points
   *
   * @param p the first point index
   * @param q the second point index
   * @return distance
   */
  public double d(final int p, final int q) {
    return points.get(p).d(points.get(q));
  }

  /**
   * Computes distance between two points
   *
   * @param p the first point
   * @param q the second point index
   * @return distance
   */
  public double d(final P p, final int q) {
    return p.d(points.get(q));
  }

  /**
   * Computes distance between a point and a cluster
   *
   * @param p the point
   * @param c the cluster
   * @return sum of distances
   */
  public double d(final P p, final Cluster c) {
    return c.stream()
      .mapToDouble(q -> d(p, q))
      .sum();
  }

  /**
   * Computes distance between a point and a cluster
   *
   * @param p the point index
   * @param c the cluster
   * @return sum of distances
   */
  public double d(final int p, final Cluster c) {
    return d(points.get(p), c);
  }

  /**
   * Compute the exact weights of a cluster
   *
   * @param cluster reference cluster
   * @return vector of weights
   **/
  public final double[] get_weights(Cluster cluster) {
    return points.stream().mapToDouble(e -> d(e, cluster)).toArray();
  };

  /**
   * Compute the weights of a cluster
   *
   * @param cluster reference cluster index
   * @return vector of weights
   **/
  public final double[] get_weights(int cluster) {
    return get_weights(clusters.get(cluster));
  };

  /**
   * Computes silhouette of current clustering
   *
   * @return silhouette of the clustering
   */
  public final double silhouette() {
    double[][] weights = clusters.stream()
      .map(this::get_weights)
      .toArray(double[][]::new);
    return Sil.sil(this, weights);
  }
}

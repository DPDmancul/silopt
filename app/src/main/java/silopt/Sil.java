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
import java.util.stream.IntStream;
import java.util.function.Function;

import org.apache.commons.math3.random.MersenneTwister;

/** Utilities to compute (approximated) silhouette */
public class Sil {

  /**
   * Compute a PPS sample of a cluster
   *
   * @param cluster the cluster to sample from
   * @param points the vector of points of the clustering
   * @param k the number of clusters
   * @param t the expected sample size
   * @param delta the required probability
   * @param seed seed for MersenneTwister random generator
   * @return a PPS sample of cluster of expected size t. The returned object is a list of pairs (index of point, probability).
   */
  public final static <P extends Point<P>> List<Pair<Integer, Double>> sample(final Cluster cluster, final List<P> points, final int k, final int t, final double delta, int seed) {
    if (cluster.size() <= t)
      return cluster
        .stream()
        .map(p -> new Pair<>(p, 1.))
        .toList();

    MersenneTwister rand = new MersenneTwister(seed);

    double unif = 1. / cluster.size();
    Bernoulli prob = new Bernoulli(rand, 2. * unif * Math.log(2. * k / delta));

    List<Integer> s_0 = cluster
      .stream()
      .filter(e -> prob.sample())
      .toList();

    List<Double> weights = s_0
      .stream()
      .map(p -> cluster
        .stream()
        .mapToDouble(e -> points.get(p).d(points.get(e)))
        .sum()
      )
      .toList();

    Function<Integer, Double> probability = p -> Math.min(1., t * Math.max(unif,
      //γ
      IntStream.range(0, s_0.size())
        .mapToDouble(i -> points.get(p).d(points.get(s_0.get(i))) / weights.get(i))
        .max().orElse(0.)
    ));

    return cluster.stream()
      .map(p -> new Pair<>(p, probability.apply(p)))
      .filter(pair -> Bernoulli.sample(rand, pair.getValue()))
      .toList();
  }

  /**
   * Compute the silhouette of a clustering
   *
   * @param clustering the clustering to be evaluated
   * @param weights a matrix of weights involving clusters and points
   * @return the silhouette
   */
  public final static <P extends Point<P>> double sil(final Clustering<P> clustering, final double[][] weights) {
    return IntStream.range(0, clustering.clusters.size())
      .mapToDouble(i -> {
        final Cluster c_i = clustering.clusters.get(i);
        if (c_i.size() <= 1)
          return 0.;
        double a_den = c_i.size() - 1;
        return c_i.stream()
          .mapToDouble(x -> {
            final double a = weights[i][x] / a_den;
            final double b = IntStream.range(0, clustering.clusters.size())
              .filter(j -> i != j)
              .mapToDouble(j -> weights[j][x] / clustering.clusters.get(j).size())
              .min().orElse(0.);
            return (b - a) / Math.max(a, b);
          })
          .sum();
      })
      .sum() / clustering.points.size();
  }
}

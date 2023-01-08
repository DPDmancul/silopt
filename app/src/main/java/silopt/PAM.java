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

import silopt.Clustering.Cluster;
import silopt.Clustering.Clustering;
import silopt.Clustering.Point;

import java.util.List;
import java.util.ArrayList;
import java.util.stream.IntStream;
import java.util.stream.Collectors;

import javafx.util.Pair;

/** k-medoids algorithm (Partitioning around medoids) */
public class PAM<P extends Point<P>> extends Optimizable<P> {
  private final int k;
  private final int n;
  private final List<Integer> medoids;

  PAM(List<P> points, Cluster centers) {
    super(new Clustering<P>());
    k = centers.size();
    n = points.size();
    clustering.points = points;
    clustering.clusters = IntStream.range(0, k)
      .mapToObj(i -> new Cluster())
      .toList();
    medoids = new ArrayList<>(centers);
  }

  public Pair<Double, Integer> optimize() {

    int iterations = 0;

    double old_obj;
    double new_obj = Double.MAX_VALUE;
    List<Integer> clusters = null;

    do {
      ++ iterations;
      old_obj = new_obj;

      int best_i= -1;
      int best_e = -1;

      for (int i = 0; i < k; ++i)
        for (int e = 0; e < n; ++e) {
          if(medoids.contains(e)) continue;

          // System.err.printf("i=%d, e=%d%n", i, e);

          // swap medoids[i] and e
          int old_medoid = medoids.get(i);
          medoids.set(i, e);

          // compute new objective function value
          Pair<List<Integer>, Double> obj = clustering.points.stream()
            .map(p -> IntStream.range(0,k)
              .mapToObj(j -> new Pair<>(j, clustering.d(p, medoids.get(j))))
              .min((a,b) -> (int) Math.signum(a.getValue() - b.getValue()))
              .get()
            )
            .collect(Collectors.teeing(
              Collectors.mapping(Pair::getKey, Collectors.toList()),
              Collectors.summingDouble(pair -> pair.getValue()),
              Pair::new
            ));

          if(obj.getValue() < new_obj) {
            new_obj = obj.getValue();
            clusters = obj.getKey();
            best_i = i;
            best_e = e;
          }

          // restore medoids[i]
          medoids.set(i, old_medoid);
        }

      if (best_i >= 0 && best_e >= 0)
        medoids.set(best_i, best_e);

      System.err.println("Iteration: " + iterations + " (" + new_obj + ")\t\t\t ");

    } while (new_obj < old_obj);

    for (int i = 0; i < n; ++i)
      clustering.clusters.get(clusters.get(i)).add(i);

    return new Pair<>(new_obj, iterations);
  }
}

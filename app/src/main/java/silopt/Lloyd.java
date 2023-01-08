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
import java.util.stream.IntStream;

import javafx.util.Pair;

public class Lloyd<P extends Point<P>> extends Optimizable<P> {
  private final int k;
  private final int n;

  Lloyd(Clustering<P> clustering) {
    super(clustering);
    k = clustering.clusters.size();
    n = clustering.points.size();
  }

  public Pair<Double, Integer> optimize() {

    int iterations = 0;

    double old_obj;
    double new_obj = Double.MAX_VALUE;

    do {
      ++ iterations;
      old_obj = new_obj;
      new_obj = 0;

      // Compute centers
      List<P> centers = clustering.clusters.stream()
        .map(c -> c.stream()
          .map(clustering.points::get)
          .reduce(Point::add).get()
          .ratio(c.size())
        )
        .toList();

      clustering.clusters = IntStream.range(0, k)
        .mapToObj(i -> new Cluster())
        .toList();

      // Assign points to nearest cluster
      for (int e = 0; e < n; ++e) {
        int argmin = -1;
        double min = Double.MAX_VALUE;

        for (int i = 0; i < k; ++i) {
          final double d = clustering.d(centers.get(i), e);
          if (d < min){
            min = d;
            argmin = i;
          }
        }

        new_obj += Math.pow(min, 2);
        clustering.clusters.get(argmin).add(e);
      }

      System.err.println("Iteration: " + iterations + " (" + new_obj + ")\t\t\t ");

    } while (new_obj < old_obj);

    return new Pair<>(new_obj, iterations);
  }
}

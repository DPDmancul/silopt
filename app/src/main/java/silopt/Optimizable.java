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
import silopt.Clustering.Point;

import javafx.util.Pair;

/**
 * Generic abstract class for all optimization algorithms
 */
public abstract class Optimizable<P extends Point<P>> {

  protected final Clustering<P> clustering;

  Optimizable(Clustering<P> c) {
    clustering = c;
  }

  /**
   * Optimize in place the clustering for the silhouette
   *
   * @return internal silhouette of the optimized clustering and number of iterations
   */
  public abstract Pair<Double, Integer> optimize();

  /**
   * Computes exact silhouette of current clustering
   *
   * @return exact silhouette of the clustering
   */
  public double silhouette(){
    return clustering.silhouette();
  }
}


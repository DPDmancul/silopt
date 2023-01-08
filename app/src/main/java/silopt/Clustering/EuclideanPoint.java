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

/** Vector representing a point in an Euclidean space */
public class EuclideanPoint extends CartesianPoint<EuclideanPoint> {

  public EuclideanPoint(double... p) {
    super(p);
  }

  public EuclideanPoint(CartesianVector v) {
    super(v);
  }

  /**
   * Creates new EuclideanPoint from a CartesianVector
   *
   * @param v the CartesianVector
   */
  public EuclideanPoint fromCartesianVector(CartesianVector v){
    return new EuclideanPoint(v);
  }

  /** Build an EuclideanPoint  */
  protected EuclideanPoint builder(double... p) {
    return new EuclideanPoint(p);
  }

  /**
   * Compute distance w.r.t. another point
   *
   * If the two points have different size the program will exit
   *
   * @param other the other point
   * @return distance
   */
  public final double d(final EuclideanPoint other) {
    assert_same_dim(other);

    double sqdist = 0;
    for (int i = 0; i < dim(); ++i) {
      double diff = vector[i] - other.vector[i];
      sqdist += diff * diff;
    }

    return Math.sqrt(sqdist);
  }
}

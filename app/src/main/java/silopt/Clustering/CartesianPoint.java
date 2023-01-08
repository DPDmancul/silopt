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

/** Vector representing a point in an real space */
public abstract class CartesianPoint<This extends CartesianPoint<This>> extends CartesianVector implements Point<This> {
  public CartesianPoint(double... p) {
    super(p);
  }

  public CartesianPoint(CartesianVector v) {
    super(v.vector);
  }

  /**
   * Creates new CartesianPoint from a CartesianVector
   *
   * usually returns `new CartesianPoint(v)`
   *
   * @param v the CartesianVector
   */
  public abstract This fromCartesianVector(CartesianVector v);

  /**
   * Asserts the dimension of this point is the same of another point
   *
   * @param other the point to compare dimension with
   */
  public final void assert_same_dim(final This other) {
    assert_same_dim(other, "Points");
  }

  /**
   * Compute distance w.r.t. another point
   *
   * If the two points have different size the program will exit
   *
   * @param other the other point
   * @return distance
   */
  public abstract double d(final This other);

  /**
   * Computes sum with another point
   *
   * If the two points have different size the program will exit
   *
   * @param other the other point
   * @return sum
   */
  public This add(final This other){
    return fromCartesianVector(super.add(other));
  }

  /**
   * Computes ratio by a scalar of this point
   *
   * @param divisor the scalar divisor
   * @return ratio
   */
  public This ratio(final double divisor) {
    return fromCartesianVector(super.ratio(divisor));
  }
}

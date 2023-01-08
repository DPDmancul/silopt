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

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/** Vector in a real space */
public class CartesianVector {
  protected double[] vector;

  public CartesianVector(double... p) {
    vector = p;
  }

  /**
   * Get a component of the vector
   *
   * @param i component index
   * @return the requested component
   */
  public final double get(int i) {
    return vector[i];
  }
  /**
   * Set a component of the vector
   *
   * @param i component index
   * @param v the value to give to the component
   */
  public final void set(int i, double v) {
    vector[i] = v;
  }

  /**
   * Get the dimension of the vector
   */
  public final int dim() {
    return vector.length;
  }

  /**
   * Asserts the dimension of this vector is the same of another
   *
   * @param other the vector to compare dimension with
   * @param name name of the objects (e.g. Vectors, Points)
   */
  public final void assert_same_dim(final CartesianVector other, final String name) {
    if (vector.length != other.vector.length) {
      System.err.printf("%s with different size: %d and %d%n", name, vector.length, other.vector.length);
      System.exit(1);
    }
  }
  public final void assert_same_dim(final CartesianVector other) {
    assert_same_dim(other, "Vectors");
  }

  /**
   * Computes sum with another vector
   *
   * If the two vectors have different size the program will exit
   *
   * @param other the other vector
   * @return sum
   */
  public CartesianVector add(final CartesianVector other){
    assert_same_dim(other);

    return new CartesianVector(
      IntStream.range(0, vector.length)
        .mapToDouble(i -> vector[i] + other.vector[i])
        .toArray()
    );
  }

  /**
   * Computes ratio by a scalar of this vector
   *
   * @param divisor the scalar divisor
   * @return ratio
   */
  public CartesianVector ratio(final double divisor){
    return new CartesianVector(
      DoubleStream.of(vector)
        .map(x -> x / divisor)
        .toArray()
    );
  }

  @Override
  public String toString() {
    return Arrays.toString(vector);
  }

  public String toDatasetRow() {
    return DoubleStream.of(vector)
      .mapToObj(Double::toString)
      .collect(Collectors.joining(","));
  }
}

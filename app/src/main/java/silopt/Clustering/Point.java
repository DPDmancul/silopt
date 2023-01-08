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

public interface Point<This> {
  /**
   * Computes distance with another point
   *
   * @param other the other point
   * @return distance
   */
  public double d(final This other);

  /**
   * Computes sum with another points
   *
   * @param other the other point
   * @return sum
   */
  public This add(final This other);

  /**
   * Computes ratio by a scalar of this point
   *
   * @param divisor the scalar divisor
   * @return ratio
   */
  public This ratio(final double divisor);
}

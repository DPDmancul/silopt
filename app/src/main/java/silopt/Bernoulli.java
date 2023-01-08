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

import org.apache.commons.math3.random.BitsStreamGenerator;
import org.apache.commons.math3.random.MersenneTwister;

/** Bernoulli distribution */
public class Bernoulli {

  /** success probability */
  public double p;

  private BitsStreamGenerator rand;

  /**
   * Create a Bernoulli distribution
   *
   * If the probability is less than zero it will be treated as zero
   * If the probability is greater than one it will be treated as one
   *
   * @param rand random generator
   * @param p success probability
   **/
  public Bernoulli(BitsStreamGenerator rand, double p){
    this.p = p;
    this.rand = rand;
  }

  /**
   * Create a Bernoulli distribution
   *
   * If the probability is less than zero it will be treated as zero
   * If the probability is greater than one it will be treated as one
   *
   * @param p success probability
   * @param seed seed for MersenneTwister random generator
   **/
  public Bernoulli(double p, int seed){
    this.p = p;
    this.rand = new MersenneTwister(seed);
  }

  /**
   * Sample a Ber(p)
   *
   * If the probability is less than zero it will be treated as zero
   * If the probability is greater than one it will be treated as one
   *
   * @param p success probability
   * @return a sampled Bernoulli value
   **/
  public final static boolean sample(BitsStreamGenerator rand, double p) {
    return rand.nextDouble() < p;
  }

  /**
   * Sample a Ber(p)
   *
   * @return a sampled Bernoulli value
   **/
  public final boolean sample() {
    return sample(this.rand, this.p);
  }
}

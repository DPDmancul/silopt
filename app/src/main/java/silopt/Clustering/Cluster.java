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

import java.util.HashSet;

/** Each cluster contains the indexes of the contained points. */
public class Cluster extends HashSet<Integer> {
  public Cluster() {
    super();
  }

  public Cluster(int... points) {
    for (int p : points)
      this.add(p);
  }
}

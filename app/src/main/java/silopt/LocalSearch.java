package silopt;

import silopt.Clustering.Clustering;
import silopt.Clustering.Cluster;
import silopt.Clustering.Point;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import javafx.util.Pair;

public class LocalSearch {

  /** Record the moving of a point form cluster Ci to Cj */
  private record Move(int i, int j, int point) {}

  /**
   * Abstract class with the core of all algorithms
   */
  private abstract static class BaseAlg<P extends Point<P>> extends Optimizable<P> {
    protected final int k;
    protected final int n;
    protected final double heta;
    protected final double eps_bar;

    /**
     * Points which will be moved to another cluster, organized by their clusters
     * - it contains all the clusters, in standard cases
     * - it contains the samples of the clusters, in the coreset based algorithms
     */
    protected List<Cluster> movable;

    /**
     * @param clustering clustering to optimize
     * @param heta local search minimum improvement
     */
    public BaseAlg(Clustering<P> clustering, double heta) {
      // The false positive threshold eps_bar is set to zero
      // since we a re using exact silhouette
      this(clustering, heta, 0.);
    }

    protected BaseAlg(Clustering<P> clustering, double heta, double eps_bar) {
      this(clustering, heta, nk -> eps_bar);
    }

    protected BaseAlg(Clustering<P> clustering, double heta, Function<Integer, Double> eps_bar) {
      super(clustering);
      this.heta = heta;

      /** Number of clusters */
      k = clustering.clusters.size();
      /** Number of points */
      n = clustering.points.size();

      this.eps_bar = eps_bar.apply(n * k);
    }

    /**
     * Compute the weights of a cluster
     *
     * @param cluster reference cluster
     * @return vector of weights
     **/
    protected double[] get_weights(Cluster cluster) {
      return clustering.get_weights(cluster);
    };

    /**
     * Compute initial weights and set movables points
     *
     * @returns initial weights
     */
    protected double[][] init() {
      movable = clustering.clusters;

      return clustering.clusters.stream()
        .map(this::get_weights)
        .toArray(double[][]::new);
    }

    /**
     * Compute new weights between a cluster and a point
     *
     * @param cluster the reference cluster
     * @param old the old weights
     * @param moved the moved point
     * @param added states if the point was added or removed from the cluster
     * @returns updated weights
     */
    protected double[] update_weights(Cluster cluster, double[] old, int moved, boolean added) {
      final P x = clustering.points.get(moved);

      return IntStream.range(0, n)
        .mapToDouble(i -> added ?
            old[i] + clustering.d(x, i)
          :
            old[i] - clustering.d(x, i)
        )
        .toArray();
    }

    /**
     * What to do after applying a moving
     * (e.g. update samples)
     */
    protected void apply(Move move) {}

    /**
     * Optimize in place the clustering for the silhouette
     *
     * @return internal silhouette of the optimized clustering and number of iterations
     */
    public Pair<Double, Integer> optimize() {
      double[][] w_curr = init();

      // Silhouette of the current optimal clustering
      double[][] w_opt = w_curr.clone();
      double s_opt = Sil.sil(clustering, w_opt);
      // System.err.println(s_opt);

      double s_curr;
      int r = 0;
      do {
        // Silhouette of the reference clustering
        s_curr = s_opt;

        ++r;
        System.err.println("Iteration: " + r + " (" + s_curr + ")\t\t\t ");

        // Current optimal clustering
        // Only the point move is stored
        Move c_opt = null;

        int moved_points = 0;

        for (int i = 0; i < k; ++i) {
          final Cluster c_i = clustering.clusters.get(i);
          final double[] w_i = w_curr[i];
          // For each point x in the clustering
          // toArray avoids ConcurrentModificationException
          for (final int x : movable.get(i).toArray(Integer[]::new)) {
            ++moved_points;

            c_i.remove(x); // Ci'
            w_curr[i] = update_weights(c_i, w_curr[i], x, false); // W'

            // For each cluster, except the one containing x
            for (int j = 0; j < k; ++j) {
              if (j == i) continue;
              // System.err.printf("Moving: %d -> %d (%.2f%%)\t\t \r", i, j, (100. * moved_points) / n);

              final Cluster c_j = clustering.clusters.get(j);
              final double[] w_j = w_curr[j];

              c_j.add(x); //Cj'
              w_curr[j] = update_weights(c_j, w_curr[j], x, true); // W'

              double s_prime = Sil.sil(clustering, w_curr);

              if (s_prime > s_opt && s_prime > s_curr + eps_bar) {
                // Move x from Ci to Cj
                c_opt = new Move (i, j, x);
                w_opt = w_curr.clone();
                s_opt = s_prime;
              }

              // Restore Cj
              c_j.remove(x);
              w_curr[j] = w_j;
            }

            // Restore Ci
            c_i.add(x);
            w_curr[i] = w_i;
          }
        }

        // Update the reference clustering for the next iteration
        if (c_opt != null) {
          // System.err.printf("Moving %s from C%d to C%d; sil = %f; impr = %f%n", clustering.points.get(c_opt.point()).toString(), c_opt.i(), c_opt.j(), s_opt, s_opt - s_curr);
          clustering.clusters.get(c_opt.i()).remove(c_opt.point());
          clustering.clusters.get(c_opt.j()).add(c_opt.point());
          apply(c_opt);
          // System.err.println(s_opt + "\t" + silhouette());
        }
        w_curr = w_opt.clone();


      } while (s_opt - s_curr > heta + eps_bar);

      if (s_opt == s_curr)
        System.out.printf("Local maximum (internal metric = %f)%n", s_opt);
      else
        System.out.printf("Improvement under threshold (improvement = %f)%n", s_opt - s_curr);

      return new Pair<>(s_opt, r);
    }
  }

  /**
   * Extended core for algorithms which use samples
   */
  private abstract static class BaseAlgSamplable<P extends Point<P>> extends BaseAlg<P> {
    protected final int t;
    protected final double delta;
    protected final int seed;

    /** Suitable constant */
    protected static final double c = 0; // Do not use false positive threshold

    /**
     * @param clustering clustering to optimize
     * @param t sample expected size
     * @param delta probability bound
     * @param heta local search minimum improvement
     * @param approx whether to use or not approximated silhouette
     * @param seed seed for MersenneTwister random generator
     */
    public BaseAlgSamplable(Clustering<P> clustering, final int t, final double delta, final double heta, boolean approx, int seed) {
      super(clustering, heta, nk -> {
        if (approx) {
          /** Weights relative error bound */
          // t = ceil(c/2ε² ln(4nk/δ))
          final double eps = Math.sqrt(c / (2 * t) * Math.log(4 * nk / delta));

          /** False positive threshold */
          return 8 * eps / (1 - eps);
        } else
          return 0.;
      });
      this.t = t;
      this.delta = delta;
      this.seed = seed;
    }

    /**
     * Samples a given cluster
     *
     * @param cluster the cluster to be sampled
     * @return the sample with probabilities
     */
    protected List<Pair<Integer, Double>> sample(Cluster cluster) {
      return Sil.sample(cluster, clustering.points, k, t, delta, seed);
    }

    /**
     * Converts a sample in a cluster
     *
     * @param sample the sample
     * @return the sample in cluster format
     */
    protected final Cluster sample_as_cluster(List<Pair<Integer, Double>> sample) {
      return sample.stream()
        .map(Pair::getKey)
        .collect(Collectors.toCollection(Cluster::new));
    }

    /**
     * Samples a given cluster
     *
     * @param cluster the cluster to be sampled
     * @return the sample in cluster format
     */
    protected final Cluster sample_as_cluster(Cluster cluster) {
      return sample_as_cluster(sample(cluster));
    }

    /**
     * Compute the approximated weights from a sample
     *
     * @param sample sample of the cluster
     * @return vector of approximated weights
     **/
    protected double[] get_weights(List<Pair<Integer, Double>> sample) {
      return IntStream.range(0, n).mapToDouble(i ->
        sample.stream()
          .mapToDouble(pair -> clustering.d(i, pair.getKey()) / pair.getValue())
          .sum()
      ).toArray();
    };
  }

  public static class ExactMemoization<P extends Point<P>> extends BaseAlg<P> {
    /**
     * @param clustering clustering to optimize
     * @param heta local search minimum improvement
     */
    public ExactMemoization(Clustering<P> clustering, double heta) {
      super(clustering, heta);
    }
  }

  public static class ApproxMemoization<P extends Point<P>> extends BaseAlgSamplable<P> {
    /**
     * @param clustering clustering to optimize
     * @param t sample expected size
     * @param delta probability bound
     * @param heta local search minimum improvement
     * @param seed seed for MersenneTwister random generator
     */
    public ApproxMemoization(Clustering<P> clustering, final int t, final double delta, final double heta, final int seed) {
      super(clustering, t, delta, heta, true, seed);
    }

    // Use approximated weights
    @Override
    protected double[] get_weights(Cluster cluster) {
      return get_weights(sample(cluster));
    }

    // Recompute approximate weights at update
    @Override
    protected double[] update_weights(Cluster cluster, double[] old, int moved, boolean added) {
      return get_weights(cluster);
    }
  }

  public static class ExactCoresetMemoization<P extends Point<P>> extends BaseAlgSamplable<P> {
    /**
     * @param clustering clustering to optimize
     * @param t sample expected size
     * @param delta probability bound
     * @param heta local search minimum improvement
     * @param seed seed for MersenneTwister random generator
     */
    public ExactCoresetMemoization(Clustering<P> clustering, final int t, final double delta, final double heta, final int seed) {
      this(clustering, t, delta, heta, false, seed);
    }

    protected ExactCoresetMemoization(Clustering<P> clustering, final int t, final double delta, final double heta, boolean approx, int seed) {
      super(clustering, t, delta, heta, approx, seed);
    }

    @Override
    protected double[][] init() {
      // Use coreset
      movable = clustering.clusters.stream()
        .map(super::sample_as_cluster)
        .collect(Collectors.toList());

      // Compute exact weights
      return clustering.clusters.stream()
        .map(super::get_weights)
        .toArray(double[][]::new);
    }

    @Override
    protected void apply(Move move) {
      int[] indexes = {move.i(), move.j()};
      for (int i : indexes)
        movable.set(i, sample_as_cluster(clustering.clusters.get(i)));
    }

  }

  public static class ApproxCoresetMemoization<P extends Point<P>> extends ExactCoresetMemoization<P> {
    /**
     * @param clustering clustering to optimize
     * @param t sample expected size
     * @param delta probability bound
     * @param heta local search minimum improvement
     * @param seed seed for MersenneTwister random generator
     */
    public ApproxCoresetMemoization(Clustering<P> clustering, final int t, final double delta, final double heta, final int seed) {
      super (clustering, t, delta, heta, true, seed);
    }

    @Override
    protected double[][] init() {
      // sample -> approximated weights
      List<List<Pair<Integer, Double>>> samples = clustering.clusters.stream()
        .map(super::sample)
        .toList();

      // Use coreset
      movable = samples.stream()
        .map(super::sample_as_cluster)
        .collect(Collectors.toList());

      // Reuse sample to compute approximated weights
      return samples.stream()
        .map(super::get_weights)
        .toArray(double[][]::new);
    }
  }

  public static class SampleOpt<P extends Point<P>> extends BaseAlgSamplable<P> {
    /**
     * @param clustering clustering to optimize
     * @param t sample expected size
     * @param delta probability bound
     * @param heta local search minimum improvement
     * @param seed seed for MersenneTwister random generator
     */
    public SampleOpt(Clustering<P> clustering, final int t, final double delta, final double heta, int seed) {
      super(clustering, t, delta, heta, false, seed);
    }

    /**
     * Optimize in place the clustering for the silhouette
     *
     * WARNING: return silhouette of sample, not of clustering
     *
     * @return internal silhouette of the optimized clustering and number of iterations
     */
    @Override
    public Pair<Double, Integer> optimize() {

      // Compute sample
      movable = clustering.clusters.stream()
        .map(super::sample_as_cluster)
        .collect(Collectors.toList());

      List<Integer> point_translation = new ArrayList<>();
      Clustering<P> sample = new Clustering<P>();
      boolean[] sampled = new boolean[clustering.points.size()];
      for (final Cluster c : movable) {
        Cluster cluster = new Cluster();
        for (final int i : c) {
          sampled[i] = true;
          cluster.add(sample.points.size());
          point_translation.add(i);
          sample.points.add(clustering.points.get(i));
        }
        sample.clusters.add(cluster);
      }

      // Apply algorithm for optimizing Exact silhouette on sample
      Pair<Double, Integer> result = new ExactMemoization<P>(sample, heta).optimize();

      clustering.clusters = sample.clusters.stream()
        .map(c -> c.stream()
          .map(point_translation::get)
          .collect(Collectors.toCollection(Cluster::new))
        )
        .toList();

      // Find the cluster which minimizes the mean distance
      int[] index = new int[clustering.points.size()];

      for (int e = 0; e < clustering.points.size(); ++e) {
        if (sampled[e]) continue;

        int argmin = -1;
        double min = Double.MAX_VALUE;

        for (int i = 0; i < k; ++i) {
          final Cluster c_i = clustering.clusters.get(i);
          final double d = clustering.d(e, c_i) / c_i.size();
          if (d < min){
            min = d;
            argmin = i;
          }
        }
        index[e] = argmin;
      }

      // Assign not sampled points to nearest cluster
      for (int e = 0; e < clustering.points.size(); ++e)
        if (!sampled[e])
          clustering.clusters.get(index[e]).add(e);

      return result;
    }

  }
}

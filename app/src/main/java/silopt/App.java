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

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;
import javafx.util.Pair;

import com.opencsv.CSVReader;

import org.apache.commons.math3.random.MersenneTwister;

import silopt.Clustering.Clustering;
import silopt.Clustering.Cluster;
import silopt.Clustering.Point;
import silopt.Clustering.CartesianVector;
import silopt.Clustering.CartesianPoint;
import silopt.Clustering.EuclideanPoint;
import silopt.Clustering.ManhattanPoint;

public class App {

  /**
   * Run an algorithm and print the results
   *
   * @param name The name of the algorithm to be printed
   * @param alg The class implementing the algorithm
   */
  private static <P extends Point<P>> void run(String name, Optimizable<P> alg) {
      System.err.println("\n" + name);

      // Execute the algorithm
      long startTime = System.currentTimeMillis();
      // result: (final value of internal optimized metric, number of iterations)
      Pair<Double, Integer> result = alg.optimize();
      long endTime = System.currentTimeMillis();

      // Compute the silhouette of the resulting clustering
      double sil = alg.silhouette();

      System.out.print(name + "\t\t");
      if (result == null)
        System.out.print("-\t\t\t" + sil + "\t-");
      else
        System.out.print(result.getKey() + "\t" + sil + "\t" + result.getValue());
      System.out.println("\t\t" + (endTime - startTime) + " ms");
  }

  /**
   * Obtain points from dataset file
   *
   * @param dataset Name of the dataset, stored in `datasets/<dataset>.data`
   * @param limit Maximum number of points to load
   * @param point_type The class into which represent points. The distance used to compute distances will depend on this class.
   *                   Must have a constructor accepting a `CartesianVector`.
   * @return List of points
   */
  private static <P extends CartesianVector> List<P> get_points(String dataset, long limit, Class<P> point_type) {
    List<P> result = null;
    // open dataset file
    try (CSVReader file = new CSVReader(new FileReader("../datasets/" + dataset + ".data"))) {
      result = StreamSupport.stream(file.spliterator(), false)
        .map(line -> {
          // each line represents a point
          P point = null;
          try{
            // construct a point for the given class
            point = point_type.getConstructor(CartesianVector.class).newInstance(
              // parse the line as a Cartesian vector
              new CartesianVector(
                Arrays.stream(line)
                  .mapToDouble(Double::parseDouble)
                  .toArray()
                )
              );
          } catch (Exception e) {
            System.err.println("Error casting vector to point");
            System.err.println(e);
            System.exit(-1);
          }
          return point;
        })
        .limit(limit)
        .toList();
    } catch (Exception e) {
      System.err.println("Error reading points");
      System.err.println(e);
      System.exit(-1);
    }
    return result;
  }

  /** Get all points from the dataset */
  private static <P extends CartesianVector> List<P> get_points(String dataset, Class<P> point_type) {
    return get_points(dataset, Long.MAX_VALUE, point_type);
  }

  /**
   * Get points for a specific cluster of a datastet
   *
   * @param dataset Name of the dataset
   * @param k The number of clusters
   * @param i The cluster number
   * @param point_type The class into which represent points
   * @return List of points
   */
  private static <P extends CartesianVector> List<P> get_points(String dataset, int k, int i, Class<P> point_type) {
    return get_points(dataset + "/" + k + "/" + i, point_type);
  }

  /**
   * Load a clustering of a dataset
   *
   * @param dataset Name of the dataset
   * @param k The number of clusters
   * @param point_type The class into which represent points
   * @return The clustering
   */
  private static <P extends CartesianPoint<P>> Clustering<P> get_clustering(String dataset, int k, Class<P> point_type) {
    return new Clustering<>(
      IntStream.range(0, k)
        .mapToObj(i -> get_points(dataset, k, i, point_type))
        .toList()
    );
  }

  /**
   * Run an experiment with all algorithms
   *
   * @param point_type The class into which represent points
   * @param dataset Name of the dataset
   * @param k The number of clusters
   * @param delta Probability of the error bound
   * @param heta Threshold value
   * @param max_n Maximum number of points
   * @param t... Value(s) of the sample expected size to test
   */
  private static <P extends CartesianPoint<P>> void run_for(Class<P> point_type, final String dataset, final int k, final double delta, final double heta, final long max_n, final int... t) {
    // Do not execute more than 3 times PAM, since it is slow
    int pam_watchdog = 3;

    // Repeat 5 times and then return median, mean and variance
    // We use fixed random seeds for reproducibility of the experiments
    for (int seed : new int[]{361471, 623465, 7147, 487722, 250750}) {
      if (pam_watchdog == 0)
        break;

      // Obtain all points of the dataset
      List<P> points = get_points(dataset, max_n, point_type);

      System.out.printf("Dataset: %s, n = %d (d = %d)%nk = %d, δ = %f, η = %f%n", dataset, points.size(), points.get(0).dim(), k, delta, heta);
      System.err.printf("Dataset: %s, n = %d%nk = %d, δ = %f, η = %f%n", dataset, points.size(), k, delta, heta);
      System.out.printf("Seed: %d%n", seed);
      System.err.printf("Seed: %d%n", seed);
      System.out.println("Optimization algorithm\tInternal metric\t\tExact silhouette\tIterations\tTime");

      // Execute k-means++
      KMeansPP<P> kmeanspp = new KMeansPP<>(points, k, seed);
      run("k-means++", kmeanspp);

      // Store k-means++ result, which is the starting point of other algorithms
      for (int i = 0; i < k; ++i) {
        try {
          File folder = new File("../datasets/" + dataset + "/" + k );
          folder.mkdirs();
          try (FileWriter file = new FileWriter(folder + "/" + i + ".data")) {
            for(int p : kmeanspp.clustering.clusters.get(i))
              file.write(points.get(p).toDatasetRow() + "\n");
          }
        } catch (Exception e) {
          System.err.println("Error saving clustering");
          System.err.println(e);
          System.exit(-1);
        }
      }

      // Run Lloyd's algorithm (or PAM if we are not using Euclidean sitances)
      if (point_type.getName().equals(EuclideanPoint.class.getName()))
        run("Lloyd k-means", new Lloyd<>(get_clustering(dataset, k, point_type)));
      else if (points.size() <= 10_000) {
        --pam_watchdog;
        run("PAM k-medoids", new PAM<>(points, kmeanspp.get_centers()));
      }

      if (k * points.size() <= 2 * 5_000)
        run("ExactMemoization", new LocalSearch.ExactMemoization<>(get_clustering(dataset, k, point_type), heta));

      // For each value of expected sample size
      for (int this_t : t) {
        System.out.println("t = " + this_t);
        System.err.println("t = " + this_t);
        if (k * points.size() <= 2 * 5_000) {
          // run("ApproxMemoization", new LocalSearch.ApproxMemoization<>(get_clustering(), T, DELTA, HETA, seed));
          run("ExactCoresetMemo", new LocalSearch.ExactCoresetMemoization<>(get_clustering(dataset, k, point_type), this_t, delta, heta, seed));
          run("ApproxCoresetMemo", new LocalSearch.ApproxCoresetMemoization<>(get_clustering(dataset, k, point_type), this_t, delta, heta, seed));
        }
        run("SampleOptimization", new LocalSearch.SampleOpt<>(get_clustering(dataset, k, point_type), this_t, delta, heta, seed));
        run("UniformSample", new LocalSearch.SampleOpt<>(get_clustering(dataset, k, point_type), this_t, delta, heta, seed) {
          // Override the sampling strategy

          protected List<Pair<Integer, Double>> sample(Cluster cluster) {
            if (cluster.size() <= t)
              return cluster
                .stream()
                .map(p -> new Pair<>(p, 1.))
                .toList();

            MersenneTwister rand = new MersenneTwister(this.seed);

            double unif = t / (double) cluster.size();
            Bernoulli prob = new Bernoulli(rand, unif);

            return cluster.stream()
              .filter(e -> prob.sample())
              .map(p -> new Pair<>(p, unif))
              .toList();
          }
        });
      }

      System.out.println("\n");
    }
  }

  public static void main(String[] args) {
    final int T = 5;
    final double DELTA = 0.1;
    final double HETA = 0.001;
    final long MAX_N = 1_000;


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////// ↓↓↓ EDIT HERE ↓↓↓ /////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////


    // https://archive.ics.uci.edu/ml/datasets/Covertype
    run_for(EuclideanPoint.class, "covtype", 5, DELTA, HETA, MAX_N, T);

    // https://archive.ics.uci.edu/ml/datasets/Query+Analytics+Workloads+Dataset
    run_for(EuclideanPoint.class, "radius", 5, DELTA, HETA, MAX_N, T);

    // https://archive.ics.uci.edu/ml/datasets/selfBACK
    run_for(EuclideanPoint.class, "selfback_026", 9, DELTA, HETA, MAX_N, T);

    // https://archive.ics.uci.edu/ml/datasets/cloud
    run_for(EuclideanPoint.class, "cloud", 10, DELTA, HETA, MAX_N, T);

    // Artificial
    run_for(EuclideanPoint.class, "concentric", 2, DELTA, HETA, 5_000, T);
    run_for(EuclideanPoint.class, "concentric", 2, DELTA, HETA, 7_000, T);
    run_for(EuclideanPoint.class, "concentric", 2, DELTA, HETA, 10_000, T, 10 * T, 100 * T);
    run_for(EuclideanPoint.class, "concentric", 2, DELTA, HETA, 20_000, T, 10 * T, 100 * T);

    final int[] t = new int[20];
    for (int i = 0; i < 20; ++i)
      t[i] = 5 + 5 * i;

    for (String dataset : new String[] {
      "radius",

      "concentric_big",
      "concentric_big_noise",

      "spherical",
      "spherical_noise",

      "gaussian_10k",
      "gaussian",
      "gaussian_noise",

      "higgs",
    })
      for (int n : new int[] {10_000, 100_000})
        for (int k : new int[] {5, 10})
          if (dataset.equals("higgs"))
            run_for(ManhattanPoint.class, dataset, k, DELTA, HETA, n/10, t);
          else
            run_for(EuclideanPoint.class, dataset, k, DELTA, HETA, n, t);
  }
}

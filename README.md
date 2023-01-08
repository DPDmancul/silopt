# A Local-Search Approach to Silhouette-Based Clustering

This repository includes a Java implementation of the algorithms proposed in the
thesis "A Local-Search Approach to Silhouette-Based Clustering" of Peressoni and supervisors Pietracaprina, Pucci, Vandin and Sarpe (<https://thesis.unipd.it/handle/20.500.12608/40252>).


## Usage

Requirements:
- `openjdk17`

1. Put datasets in `datasets` folder (see below)
2. Set up the desired experiments in `app/src/main/java/silopt/App.java` (in the `main` method)
3. Run `./gradlew run`

The results will be printed on `stdout`, while some logging will be printed on
`stderr`. To collect them into two files you can run in bash:

```bash
{ ./gradlew run | tee results.stdout; } 2> >(tee results.stderr)
```

and in fish:

```fish
begin; ./gradlew run | tee results.stdout; end 2>| tee results.stderr
```

## Datasets

The datasets have to be added into `datasets` folder, in the main dir. Each
dataset is a file with extension `.data`, which is indeed a csv. Each row
represents a point, which components are comma separated.

## Classes

- `silopt.Sil`: utilities to compute silhouette and approximated silhouette
- `silopt.Bernoulli`: implements a Bernoulli distributed random generator, based on MersenneTwister random generator
- Clustering
  * `silopt.Clustering.Clustering`: represents a clustering of a set of points
  * `silopt.Clustering.Cluster`: a cluster in a clustering
  * Points
    + `silopt.Clustering.Point`: interface for generic points
    + `silopt.Clustering.CartesianVector`: a vector in a real space
    + `silopt.Clustering.CartesianPoint`: abstract point in a real space
    + `silopt.Clustering.EuclideanPoint`: point in real space with Euclidean distance
    + `silopt.Clustering.ManhattanPoint`: point in real space with Manhattan distance
- Algorithms
  * `silopt.Optimizable`: abstract for generic optimization algorithms on
    clustering
  * `silopt.KMeansPP`: optimizable implementing k-means++ 
  * `silopt.Lloyd`: optimizable implementing Lloyd's algorithm
  * `silopt.PAM`: optimizable implementing PAM
  * `silopt.LocalSearch`: contains optimizables implementing the local-search algorithms proposed in the thesis
    + `silopt.LocalSearch.ExactMemoization`
    + `silopt.LocalSearch.ApproxMemoization`
    + `silopt.LocalSearch.ExactCoresetMemoization`
    + `silopt.LocalSearch.ApproxCoresetMemoization`
    + `silopt.LocalSearch.SampleOpt`

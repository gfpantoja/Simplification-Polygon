# Simplification-Polygon

This repository consists of complementary material from the article titled "Convex Decomposition and Simplification of a Polygon," developed by German Pantoja-Benavides, Francisco Parrerño Torres, and David Álvarez Martínez.

## Instances

The instances were generated and used by Fernández, J., Cánovas, L., & Pelegrı́n, B. (2000). Algorithms for the decomposition of a polygon into convex polygons. European Journal of Operational Research, 121(2), 330–342. https://doi.org/https://doi.org/10.1016/S0377-2217(99)00033-8

If you use these instances, please give proper credit by citing them.

## What does this project do?

The algorithms developed in the framework of the article solve two problems:

Convex Decomposition of a Polygon: The objective is to decompose a polygon in non-necessarily disjointed convex sub-polygons, such that the number of convex sub-polygons and the number of edges of those sub-polygons are minimized.
Simplification of a Polygon: The objective is to generate a polygon containing the original polygon, and the number of convex sub-polygons and edges of those sub-polygons are minimized. The simplification is controlled with an area-based similarity parameter that goes from 0% (bounding box of the original polygon) to 100% (the original polygon).

These algorithms aim to yield the minimum number of convex sub-polygons with the fewest edges. They do not return optimum solutions.

## Motivation

This project surged as a motivation to decrease the size of mathematical models and/or accelerate methodologies that solve Two-Dimensional Cutting and Packing Problems with Irregular items. These problems require that the items do not overlap, which is easily solvable for convex polygons. Therefore, irregular items are decomposed in a group of convex sub-polygons to verify the overlapping constraint. These verifications depend on the number of sub-polygons and edges of those sub-polygons.

## How to use?

The simplification algorithm has two parameters the instance (-ins) and the number of threads for multithreading (-nThreads). Some examples are shown as follows:

PiezasBurdas4 -ins P125_6 -nThreads 10

PiezasBurdas4 -ins P50_30 -nThreads 10

The instance must be a .txt file placed within a Folder called "Instances". The results are .txt files placed within a folder called "Solutions".

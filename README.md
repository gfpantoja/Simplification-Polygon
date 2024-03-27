# Simplification-Polygon

This is the repository of the article titled "Convex Decomposition and Simplification of a Polygon" developed by German Pantoja-Benavides, Francisco Parrerño Torres, and David Álvarez Martínez.

## What does this project do?

This is an algorithm that solves two problems:
1. Convex Decomposition of a polygon in which returns non-necessarily disjointed convex sub-polygons.
2. Simplification of a polygon according to an area-based similarity parameter. a 0% similarity indicates the Bounding Box of the polygon, a 100% similarity indicates the original Polygon.

These algorithms aim to yield the minimum number of convex sub-polygons with the fewest edges. These algorithms do NOT return optimum solutions.

## Motivation

This project surged as a motivation to decrease the size of mathematical models or accelerate methodologies that solve Two-Dimensional Cutting and Packing Problems with Irregular items. These problems require thar the items do not overlap, which is easyly solvable for convex polygons. Therefore the irregular items can be decomposed in a group of convex sub-polygons to verify the overlapping constraint.

## How to use?



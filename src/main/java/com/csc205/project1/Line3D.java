package com.csc205.project1;

import java.util.Objects;
import java.util.Optional;
import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Line3D represents a line segment in three-dimensional space defined by two endpoints.
 * 
 * This class demonstrates advanced object-oriented design patterns and geometric algorithms:
 * 
 * 1. COMPOSITION PATTERN: Line3D is composed of two Point3D objects, demonstrating how
 *    complex geometric entities can be built from simpler primitives. This follows the
 *    "has-a" relationship rather than inheritance.
 * 
 * 2. IMMUTABLE VALUE OBJECT: Like Point3D, Line3D is immutable after construction,
 *    ensuring thread safety and preventing accidental modification during geometric
 *    computations.
 * 
 * 3. FACTORY METHOD PATTERN: Multiple static factory methods provide different ways
 *    to construct lines while maintaining type safety and clear semantics.
 * 
 * 4. STRATEGY PATTERN: Different distance calculation methods (point-to-line,
 *    line-to-line) implement different geometric algorithms while maintaining
 *    consistent interfaces.
 * 
 * 5. OPTIONAL PATTERN: Methods that might not have valid results (like line intersections)
 *    use Optional<T> to avoid null pointer exceptions and make the API more expressive.
 * 
 * FOUNDATIONAL PRINCIPLES FOR DATA STRUCTURES AND ALGORITHMS:
 * 
 * - GEOMETRIC ALGEBRA: Implements fundamental 3D geometric algorithms including
 *   parametric line representations, vector projections, and distance calculations.
 * 
 * - NUMERICAL STABILITY: Careful handling of edge cases (parallel lines, zero-length
 *   lines) and floating-point precision issues critical for robust geometric algorithms.
 * 
 * - ALGORITHMIC EFFICIENCY: Optimized algorithms for common operations like distance
 *   calculations that avoid unnecessary computations and use geometric properties
 *   to improve performance.
 * 
 * - SPATIAL REASONING: Provides building blocks for spatial data structures like
 *   bounding volume hierarchies, spatial partitioning trees, and collision detection
 *   systems.
 * 
 * @author Generated Code
 * @version 1.0
 */
public class Line3D {
    
    private static final Logger logger = Logger.getLogger(Line3D.class.getName());
    
    // Epsilon for floating-point comparisons - critical for geometric algorithms
    private static final double EPSILON = 1e-9;
    
    private final Point3D start;
    private final Point3D end;
    
    /**
     * Creates a new Line3D with the specified start and end points.
     * 
     * This constructor follows defensive programming principles by validating inputs
     * and creating defensive copies. The validation ensures that degenerate cases
     * (null points) are caught early rather than causing failures in later geometric
     * computations.
     * 
     * @param start the starting point of the line segment
     * @param end the ending point of the line segment
     * @throws IllegalArgumentException if either point is null
     */
    public Line3D(Point3D start, Point3D end) {
        if (start == null || end == null) {
            logger.severe("Attempted to create Line3D with null points: start=" + 
                         start + ", end=" + end);
            throw new IllegalArgumentException("Line endpoints cannot be null");
        }
        
        // Defensive copying to ensure immutability
        this.start = Point3D.from(start);
        this.end = Point3D.from(end);
        
        if (isDegenerate()) {
            logger.warning("Created degenerate line (start and end points are the same)");
        }
        
        logger.info("Created Line3D from " + start + " to " + end + 
                   " with length " + String.format("%.6f", length()));
    }
    
    /**
     * Creates a line from the origin to the specified endpoint.
     * 
     * This static factory method provides a convenient way to create lines starting
     * from the origin, which is common in vector mathematics and physics simulations.
     * Factory methods improve readability and provide semantic meaning to object creation.
     * 
     * @param endpoint the end point of the line
     * @return a new Line3D from origin to the endpoint
     * @throws IllegalArgumentException if the endpoint is null
     */
    public static Line3D fromOrigin(Point3D endpoint) {
        if (endpoint == null) {
            logger.severe("Attempted to create line from origin with null endpoint");
            throw new IllegalArgumentException("Endpoint cannot be null");
        }
        
        logger.info("Creating Line3D from origin to " + endpoint);
        return new Line3D(Point3D.origin(), endpoint);
    }
    
    /**
     * Creates a line from a point in a specified direction with a given length.
     * 
     * This factory method demonstrates parametric line construction using vector
     * mathematics. It's particularly useful for algorithms that need to create
     * lines based on direction vectors, such as ray tracing or pathfinding.
     * 
     * @param start the starting point
     * @param direction the direction vector (will be normalized)
     * @param length the desired length of the line
     * @return a new Line3D in the specified direction
     * @throws IllegalArgumentException if start or direction is null, or length is negative
     */
    public static Line3D fromDirection(Point3D start, Point3D direction, double length) {
        if (start == null || direction == null) {
            logger.severe("Attempted to create line with null parameters");
            throw new IllegalArgumentException("Start point and direction cannot be null");
        }
        
        if (length < 0 || Double.isNaN(length) || Double.isInfinite(length)) {
            logger.severe("Invalid length provided: " + length);
            throw new IllegalArgumentException("Length must be a non-negative finite number");
        }
        
        if (Math.abs(direction.magnitude()) < EPSILON) {
            logger.severe("Cannot create line from zero direction vector");
            throw new IllegalArgumentException("Direction vector cannot be zero");
        }
        
        Point3D normalizedDirection = direction.normalize();
        Point3D end = start.translate(
            normalizedDirection.getX() * length,
            normalizedDirection.getY() * length,
            normalizedDirection.getZ() * length
        );
        
        logger.info("Created Line3D from direction with length " + length);
        return new Line3D(start, end);
    }
    
    /**
     * Returns the starting point of this line.
     * 
     * Since Line3D is immutable and Point3D is also immutable, we can safely
     * return the internal point reference without risk of external modification.
     * 
     * @return the starting point
     */
    public Point3D getStart() {
        return start;
    }
    
    /**
     * Returns the ending point of this line.
     * 
     * @return the ending point
     */
    public Point3D getEnd() {
        return end;
    }
    
    /**
     * Calculates and returns the length of this line segment.
     * 
     * The length is calculated using the Euclidean distance formula between
     * the start and end points. This is a fundamental operation for many
     * geometric algorithms including perimeter calculations, arc length
     * parameterization, and spatial partitioning.
     * 
     * @return the length of the line segment
     */
    public double length() {
        double len = start.distanceTo(end);
        logger.info("Calculated line length: " + len);
        return len;
    }
    
    /**
     * Returns the squared length of this line segment without taking square root.
     * 
     * This optimization is crucial for performance-critical algorithms where
     * only relative lengths matter, such as sorting line segments by length
     * or finding the shortest/longest line in a collection.
     * 
     * @return the squared length of the line segment
     */
    public double lengthSquared() {
        double lenSq = start.distanceSquaredTo(end);
        logger.info("Calculated squared line length: " + lenSq);
        return lenSq;
    }
    
    /**
     * Checks if this line is degenerate (start and end points are the same).
     * 
     * Degenerate lines can cause division by zero errors in many geometric
     * algorithms, so detecting them early is crucial for robust computational
     * geometry implementations.
     * 
     * @return true if the line is degenerate
     */
    public boolean isDegenerate() {
        boolean degenerate = start.approximatelyEquals(end);
        if (degenerate) {
            logger.info("Line is degenerate (zero length)");
        }
        return degenerate;
    }
    
    /**
     * Returns the direction vector of this line (from start to end).
     * 
     * The direction vector is fundamental for many geometric operations including
     * dot products, cross products, and parametric line equations. This method
     * returns the unnormalized direction vector.
     * 
     * @return a Point3D representing the direction vector
     */
    public Point3D getDirection() {
        Point3D direction = end.translate(-start.getX(), -start.getY(), -start.getZ());
        logger.info("Calculated direction vector: " + direction);
        return direction;
    }
    
    /**
     * Returns the normalized direction vector (unit vector) of this line.
     * 
     * Normalized direction vectors are essential for algorithms that require
     * unit vectors, such as lighting calculations, collision detection, and
     * vector field operations.
     * 
     * @return a Point3D representing the normalized direction vector
     * @throws IllegalStateException if the line is degenerate
     */
    public Point3D getUnitDirection() {
        if (isDegenerate()) {
            logger.severe("Cannot get unit direction of degenerate line");
            throw new IllegalStateException("Cannot normalize direction of degenerate line");
        }
        
        Point3D unitDir = getDirection().normalize();
        logger.info("Calculated unit direction vector: " + unitDir);
        return unitDir;
    }
    
    /**
     * Returns the midpoint of this line segment.
     * 
     * The midpoint calculation is useful for spatial partitioning algorithms,
     * bounding box calculations, and geometric center computations. It uses
     * the standard midpoint formula: midpoint = (start + end) / 2.
     * 
     * @return the midpoint of the line segment
     */
    public Point3D getMidpoint() {
        Point3D midpoint = new Point3D(
            (start.getX() + end.getX()) / 2.0,
            (start.getY() + end.getY()) / 2.0,
            (start.getZ() + end.getZ()) / 2.0
        );
        
        logger.info("Calculated midpoint: " + midpoint);
        return midpoint;
    }
    
    /**
     * Calculates a point on this line at the specified parameter t.
     * 
     * This method implements the parametric line equation: P(t) = start + t * direction
     * where t=0 gives the start point, t=1 gives the end point, and values outside
     * [0,1] extend the line beyond its endpoints. This is fundamental for many
     * geometric algorithms including ray tracing and line interpolation.
     * 
     * @param t the parameter value (0.0 = start, 1.0 = end)
     * @return the point on the line at parameter t
     * @throws IllegalArgumentException if t is NaN or infinite
     */
    public Point3D getPointAt(double t) {
        if (Double.isNaN(t) || Double.isInfinite(t)) {
            logger.severe("Invalid parameter t: " + t);
            throw new IllegalArgumentException("Parameter t must be a finite number");
        }
        
        Point3D direction = getDirection();
        Point3D point = start.translate(
            direction.getX() * t,
            direction.getY() * t,
            direction.getZ() * t
        );
        
        logger.info("Calculated point at parameter t=" + t + ": " + point);
        return point;
    }
    
    /**
     * Calculates the shortest distance from a point to this line.
     * 
     * This method implements the point-to-line distance algorithm using vector
     * projection. The algorithm:
     * 1. Project the point onto the infinite line
     * 2. Clamp the projection to the line segment
     * 3. Calculate distance from point to the clamped projection
     * 
     * This is essential for collision detection, spatial queries, and geometric
     * proximity algorithms.
     * 
     * @param point the point to calculate distance from
     * @return the shortest distance from the point to this line segment
     * @throws IllegalArgumentException if the point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        if (isDegenerate()) {
            double distance = start.distanceTo(point);
            logger.info("Distance from point to degenerate line: " + distance);
            return distance;
        }
        
        Point3D direction = getDirection();
        Point3D toPoint = point.translate(-start.getX(), -start.getY(), -start.getZ());
        
        // Project point onto the line using dot product
        double t = toPoint.dotProduct(direction) / direction.dotProduct(direction);
        
        // Clamp t to [0,1] to stay within line segment
        t = Math.max(0.0, Math.min(1.0, t));
        
        Point3D closestPoint = getPointAt(t);
        double distance = point.distanceTo(closestPoint);
        
        logger.info("Distance from point " + point + " to line: " + distance + 
                   " (closest point: " + closestPoint + ", t=" + t + ")");
        
        return distance;
    }
    
    /**
     * Finds the closest point on this line segment to a given point.
     * 
     * This method returns the actual closest point rather than just the distance,
     * which is useful for algorithms that need to know the exact location of
     * the closest approach, such as collision response and constraint solving.
     * 
     * @param point the reference point
     * @return the closest point on this line segment to the given point
     * @throws IllegalArgumentException if the point is null
     */
    public Point3D closestPointTo(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to find closest point to null");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        if (isDegenerate()) {
            logger.info("Closest point on degenerate line is the start point");
            return Point3D.from(start);
        }
        
        Point3D direction = getDirection();
        Point3D toPoint = point.translate(-start.getX(), -start.getY(), -start.getZ());
        
        double t = toPoint.dotProduct(direction) / direction.dotProduct(direction);
        t = Math.max(0.0, Math.min(1.0, t));
        
        Point3D closest = getPointAt(t);
        logger.info("Closest point on line to " + point + " is " + closest);
        
        return closest;
    }
    
    /**
     * Calculates the shortest distance between this line and another line.
     * 
     * This method implements the 3D line-to-line distance algorithm, which is
     * one of the most complex geometric calculations. The algorithm handles:
     * - Parallel lines (distance between parallel planes)
     * - Skew lines (shortest distance via common perpendicular)
     * - Intersecting lines (distance = 0)
     * 
     * This is fundamental for collision detection, interference checking, and
     * spatial relationship analysis in 3D applications.
     * 
     * @param other the other line to calculate distance to
     * @return the shortest distance between the two line segments
     * @throws IllegalArgumentException if the other line is null
     */
    public double distanceToLine(Line3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate distance to null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        // Handle degenerate cases
        if (this.isDegenerate() && other.isDegenerate()) {
            double distance = this.start.distanceTo(other.start);
            logger.info("Distance between two degenerate lines: " + distance);
            return distance;
        }
        
        if (this.isDegenerate()) {
            double distance = other.distanceToPoint(this.start);
            logger.info("Distance from degenerate line to line: " + distance);
            return distance;
        }
        
        if (other.isDegenerate()) {
            double distance = this.distanceToPoint(other.start);
            logger.info("Distance from line to degenerate line: " + distance);
            return distance;
        }
        
        Point3D d1 = this.getDirection();
        Point3D d2 = other.getDirection();
        Point3D w0 = this.start.translate(-other.start.getX(), -other.start.getY(), -other.start.getZ());
        
        double a = d1.dotProduct(d1);
        double b = d1.dotProduct(d2);
        double c = d2.dotProduct(d2);
        double d = d1.dotProduct(w0);
        double e = d2.dotProduct(w0);
        
        double denominator = a * c - b * b;
        
        // Check if lines are parallel
        if (Math.abs(denominator) < EPSILON) {
            logger.info("Lines are parallel, calculating distance between parallel lines");
            
            // Lines are parallel - find distance between them
            double t2 = Math.max(0.0, Math.min(1.0, e / c));
            Point3D point2 = other.getPointAt(t2);
            return this.distanceToPoint(point2);
        }
        
        // Lines are skew - find closest approach
        double t1 = (b * e - c * d) / denominator;
        double t2 = (a * e - b * d) / denominator;
        
        // Clamp parameters to line segments
        t1 = Math.max(0.0, Math.min(1.0, t1));
        t2 = Math.max(0.0, Math.min(1.0, t2));
        
        Point3D point1 = this.getPointAt(t1);
        Point3D point2 = other.getPointAt(t2);
        
        double distance = point1.distanceTo(point2);
        
        logger.info("Distance between lines: " + distance + 
                   " (closest points: " + point1 + " and " + point2 + ")");
        
        return distance;
    }
    
    /**
     * Finds the intersection point of this line with another line, if it exists.
     * 
     * This method determines if two line segments intersect and returns the
     * intersection point if found. The algorithm checks:
     * 1. If lines are coplanar
     * 2. If the intersection point lies within both line segments
     * 
     * Returns Optional.empty() if lines don't intersect, demonstrating modern
     * Java practices for handling potentially absent results.
     * 
     * @param other the other line to intersect with
     * @return Optional containing the intersection point, or empty if no intersection
     * @throws IllegalArgumentException if the other line is null
     */
    public Optional<Point3D> intersectionWith(Line3D other) {
        if (other == null) {
            logger.severe("Attempted to find intersection with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        if (this.isDegenerate() || other.isDegenerate()) {
            logger.info("Cannot find intersection with degenerate line");
            return Optional.empty();
        }
        
        // Check if lines intersect within tolerance
        double distance = this.distanceToLine(other);
        if (distance > EPSILON) {
            logger.info("Lines do not intersect (distance: " + distance + ")");
            return Optional.empty();
        }
        
        Point3D d1 = this.getDirection();
        Point3D d2 = other.getDirection();
        Point3D w0 = this.start.translate(-other.start.getX(), -other.start.getY(), -other.start.getZ());
        
        double a = d1.dotProduct(d1);
        double b = d1.dotProduct(d2);
        double c = d2.dotProduct(d2);
        double d = d1.dotProduct(w0);
        double e = d2.dotProduct(w0);
        
        double denominator = a * c - b * b;
        
        if (Math.abs(denominator) < EPSILON) {
            logger.info("Lines are parallel - checking for overlap");
            return Optional.empty();
        }
        
        double t1 = (b * e - c * d) / denominator;
        
        // Check if intersection is within this line segment
        if (t1 < 0.0 || t1 > 1.0) {
            logger.info("Intersection point outside first line segment");
            return Optional.empty();
        }
        
        double t2 = (a * e - b * d) / denominator;
        
        // Check if intersection is within other line segment
        if (t2 < 0.0 || t2 > 1.0) {
            logger.info("Intersection point outside second line segment");
            return Optional.empty();
        }
        
        Point3D intersection = this.getPointAt(t1);
        logger.info("Found intersection at " + intersection + " (t1=" + t1 + ", t2=" + t2 + ")");
        
        return Optional.of(intersection);
    }
    
    /**
     * Checks if this line is parallel to another line.
     * 
     * Two lines are parallel if their direction vectors are parallel (cross product
     * is zero). This check is fundamental for many geometric algorithms and
     * optimization strategies.
     * 
     * @param other the other line to check parallelism with
     * @return true if the lines are parallel
     * @throws IllegalArgumentException if the other line is null
     */
    public boolean isParallelTo(Line3D other) {
        if (other == null) {
            logger.severe("Attempted to check parallelism with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        if (this.isDegenerate() || other.isDegenerate()) {
            logger.warning("Cannot determine parallelism with degenerate line");
            return false;
        }
        
        Point3D d1 = this.getUnitDirection();
        Point3D d2 = other.getUnitDirection();
        
        // Lines are parallel if cross product is zero (within epsilon)
        Point3D cross = d1.crossProduct(d2);
        boolean parallel = cross.magnitude() < EPSILON;
        
        logger.info("Lines are " + (parallel ? "parallel" : "not parallel"));
        return parallel;
    }
    
    /**
     * Checks if this line is perpendicular to another line.
     * 
     * Two lines are perpendicular if their direction vectors are perpendicular
     * (dot product is zero). This is useful for geometric constraint solving
     * and orthogonality checks.
     * 
     * @param other the other line to check perpendicularity with
     * @return true if the lines are perpendicular
     * @throws IllegalArgumentException if the other line is null
     */
    public boolean isPerpendicularTo(Line3D other) {
        if (other == null) {
            logger.severe("Attempted to check perpendicularity with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        if (this.isDegenerate() || other.isDegenerate()) {
            logger.warning("Cannot determine perpendicularity with degenerate line");
            return false;
        }
        
        Point3D d1 = this.getDirection();
        Point3D d2 = other.getDirection();
        
        double dot = d1.dotProduct(d2);
        boolean perpendicular = Math.abs(dot) < EPSILON;
        
        logger.info("Lines are " + (perpendicular ? "perpendicular" : "not perpendicular") +
                   " (dot product: " + dot + ")");
        return perpendicular;
    }
    
    /**
     * Returns a new line that is the reverse of this line (end becomes start).
     * 
     * Line reversal is useful for path algorithms, orientation consistency,
     * and geometric transformations where direction matters.
     * 
     * @return a new Line3D with reversed direction
     */
    public Line3D reverse() {
        logger.info("Creating reversed line");
        return new Line3D(end, start);
    }
    
    /**
     * Translates this line by the specified offset vector.
     * 
     * Translation moves both endpoints by the same offset, maintaining the
     * line's direction and length. This is fundamental for coordinate system
     * transformations and spatial manipulations.
     * 
     * @param offset the translation vector
     * @return a new translated Line3D
     * @throws IllegalArgumentException if the offset is null
     */
    public Line3D translate(Point3D offset) {
        if (offset == null) {
            logger.severe("Attempted to translate line with null offset");
            throw new IllegalArgumentException("Offset vector cannot be null");
        }
        
        Point3D newStart = start.translate(offset.getX(), offset.getY(), offset.getZ());
        Point3D newEnd = end.translate(offset.getX(), offset.getY(), offset.getZ());
        
        logger.info("Translated line by offset " + offset);
        return new Line3D(newStart, newEnd);
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Line3D line3D = (Line3D) obj;
        return Objects.equals(start, line3D.start) && Objects.equals(end, line3D.end);
    }
    
    @Override
    public int hashCode() {
        return Objects.hash(start, end);
    }
    
    @Override
    public String toString() {
        return String.format("Line3D[%s -> %s, length=%.6f]", start, end, length());
    }
}
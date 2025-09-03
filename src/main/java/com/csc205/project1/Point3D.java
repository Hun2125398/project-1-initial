package com.csc205.project1;

import java.util.Objects;
import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Point3D represents a point in three-dimensional space with x, y, and z coordinates.
 * 
 * This class demonstrates several foundational object-oriented design patterns:
 * 
 * 1. VALUE OBJECT PATTERN: Point3D is immutable after construction, representing a value
 *    rather than an entity. This ensures thread safety and prevents accidental modification.
 * 
 * 2. ENCAPSULATION: Private fields with public accessor methods protect internal state
 *    while providing controlled access to coordinate values.
 * 
 * 3. BUILDER PATTERN (Static Factory Methods): Multiple constructors and static factory
 *    methods provide flexible object creation while maintaining clear interfaces.
 * 
 * 4. TEMPLATE METHOD PATTERN: Mathematical operations follow consistent patterns of
 *    validation, computation, and logging.
 * 
 * FOUNDATIONAL PRINCIPLES FOR DATA STRUCTURES AND ALGORITHMS:
 * 
 * - IMMUTABILITY: Prevents side effects in algorithms, enabling safe sharing between
 *   data structures and concurrent operations.
 * 
 * - MATHEMATICAL PRECISION: Proper handling of floating-point operations with epsilon
 *   comparisons for geometric algorithms.
 * 
 * - DEFENSIVE PROGRAMMING: Input validation and null checks prevent runtime errors
 *   in complex algorithms and data structure operations.
 * 
 * - SEPARATION OF CONCERNS: Each method has a single responsibility, making the class
 *   suitable for composition in larger algorithmic frameworks.
 * 
 * @author Generated Code
 * @version 1.0
 */
public class Point3D {
    
    private static final Logger logger = Logger.getLogger(Point3D.class.getName());
    
    // Epsilon for floating-point comparisons
    private static final double EPSILON = 1e-9;
    
    private final double x;
    private final double y;
    private final double z;
    
    /**
     * Creates a new Point3D with the specified coordinates.
     * 
     * This constructor follows the Spring Framework convention of providing clear
     * parameter documentation and validation. The immutable design ensures that
     * once created, the point's coordinates cannot be modified, which is essential
     * for geometric algorithms where point stability is crucial.
     * 
     * @param x the x-coordinate of the point
     * @param y the y-coordinate of the point
     * @param z the z-coordinate of the point
     * @throws IllegalArgumentException if any coordinate is NaN or infinite
     */
    public Point3D(double x, double y, double z) {
        if (Double.isNaN(x) || Double.isInfinite(x) ||
            Double.isNaN(y) || Double.isInfinite(y) ||
            Double.isNaN(z) || Double.isInfinite(z)) {
            logger.severe("Attempted to create Point3D with invalid coordinates: " +
                         "x=" + x + ", y=" + y + ", z=" + z);
            throw new IllegalArgumentException("Coordinates must be finite numbers");
        }
        
        this.x = x;
        this.y = y;
        this.z = z;
        
        logger.info("Created new Point3D at coordinates (" + x + ", " + y + ", " + z + ")");
    }
    
    /**
     * Creates a Point3D at the origin (0, 0, 0).
     * 
     * This static factory method provides a convenient way to create a point at the origin,
     * which is frequently used as a reference point in geometric algorithms. Following
     * Spring's naming conventions, this method clearly indicates its purpose through its name.
     * 
     * @return a new Point3D at the origin
     */
    public static Point3D origin() {
        logger.info("Creating Point3D at origin");
        return new Point3D(0, 0, 0);
    }
    
    /**
     * Creates a Point3D from another point (copy constructor pattern).
     * 
     * This method demonstrates defensive copying, a crucial pattern for maintaining
     * immutability in data structures. It's particularly useful when working with
     * algorithms that need to preserve original point data while creating variations.
     * 
     * @param other the point to copy from
     * @return a new Point3D with the same coordinates as the input point
     * @throws IllegalArgumentException if the input point is null
     */
    public static Point3D from(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to create Point3D from null reference");
            throw new IllegalArgumentException("Cannot create point from null reference");
        }
        
        logger.info("Creating Point3D copy from existing point at (" + 
                   other.x + ", " + other.y + ", " + other.z + ")");
        return new Point3D(other.x, other.y, other.z);
    }
    
    /**
     * Returns the x-coordinate of this point.
     * 
     * Simple accessor method following JavaBean conventions. The immutable design
     * means we can safely return the coordinate value without risk of external modification.
     * 
     * @return the x-coordinate
     */
    public double getX() {
        return x;
    }
    
    /**
     * Returns the y-coordinate of this point.
     * 
     * @return the y-coordinate
     */
    public double getY() {
        return y;
    }
    
    /**
     * Returns the z-coordinate of this point.
     * 
     * @return the z-coordinate
     */
    public double getZ() {
        return z;
    }
    
    /**
     * Calculates the Euclidean distance from this point to another point.
     * 
     * This method implements the fundamental distance formula in 3D space:
     * distance = √[(x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²]
     * 
     * The distance calculation is essential for many geometric algorithms including
     * nearest neighbor searches, clustering algorithms, and spatial partitioning
     * data structures like KD-trees and octrees.
     * 
     * @param other the target point to calculate distance to
     * @return the Euclidean distance between this point and the target point
     * @throws IllegalArgumentException if the target point is null
     */
    public double distanceTo(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        
        logger.info("Calculated distance from (" + this.x + ", " + this.y + ", " + this.z + 
                   ") to (" + other.x + ", " + other.y + ", " + other.z + ") = " + distance);
        
        return distance;
    }
    
    /**
     * Calculates the squared distance to another point without taking the square root.
     * 
     * This optimization is crucial for performance-critical algorithms where only
     * relative distances matter (like sorting by distance). Avoiding the square root
     * operation can significantly improve performance in data structure operations
     * such as priority queues and comparison-based sorting algorithms.
     * 
     * @param other the target point
     * @return the squared Euclidean distance
     * @throws IllegalArgumentException if the target point is null
     */
    public double distanceSquaredTo(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate squared distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        
        double distanceSquared = dx * dx + dy * dy + dz * dz;
        
        logger.info("Calculated squared distance: " + distanceSquared);
        
        return distanceSquared;
    }
    
    /**
     * Calculates the Manhattan distance (L1 distance) to another point.
     * 
     * Manhattan distance is the sum of absolute differences of coordinates:
     * distance = |x₁-x₂| + |y₁-y₂| + |z₁-z₂|
     * 
     * This distance metric is useful in algorithms where movement is constrained
     * to axis-aligned directions, and in certain machine learning algorithms
     * where L1 distance provides better performance than Euclidean distance.
     * 
     * @param other the target point
     * @return the Manhattan distance
     * @throws IllegalArgumentException if the target point is null
     */
    public double manhattanDistanceTo(Point3D other) {
        if (other == null) {
            logger.severe("Attempted to calculate Manhattan distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        double distance = Math.abs(this.x - other.x) + 
                         Math.abs(this.y - other.y) + 
                         Math.abs(this.z - other.z);
        
        logger.info("Calculated Manhattan distance: " + distance);
        
        return distance;
    }
    
    /**
     * Rotates this point around the X-axis by the specified angle and returns a new point.
     * 
     * Rotation around the X-axis follows the right-hand rule:
     * y' = y*cos(θ) - z*sin(θ)
     * z' = y*sin(θ) + z*cos(θ)
     * x remains unchanged
     * 
     * This method is fundamental for 3D transformations in computer graphics,
     * robotics, and geometric algorithms. The immutable design ensures the
     * original point remains unchanged, which is essential for transformation chains.
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateX(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle provided: " + angleRadians);
            return Point3D.from(this);
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newY = y * cos - z * sin;
        double newZ = y * sin + z * cos;
        
        logger.info("Rotating point around X-axis by " + Math.toDegrees(angleRadians) + " degrees");
        
        return new Point3D(x, newY, newZ);
    }
    
    /**
     * Rotates this point around the Y-axis by the specified angle and returns a new point.
     * 
     * Rotation around the Y-axis:
     * x' = x*cos(θ) + z*sin(θ)
     * z' = -x*sin(θ) + z*cos(θ)
     * y remains unchanged
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateY(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle provided: " + angleRadians);
            return Point3D.from(this);
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = x * cos + z * sin;
        double newZ = -x * sin + z * cos;
        
        logger.info("Rotating point around Y-axis by " + Math.toDegrees(angleRadians) + " degrees");
        
        return new Point3D(newX, y, newZ);
    }
    
    /**
     * Rotates this point around the Z-axis by the specified angle and returns a new point.
     * 
     * Rotation around the Z-axis:
     * x' = x*cos(θ) - y*sin(θ)
     * y' = x*sin(θ) + y*cos(θ)
     * z remains unchanged
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateZ(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle provided: " + angleRadians);
            return Point3D.from(this);
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = x * cos - y * sin;
        double newY = x * sin + y * cos;
        
        logger.info("Rotating point around Z-axis by " + Math.toDegrees(angleRadians) + " degrees");
        
        return new Point3D(newX, newY, z);
    }
    
    /**
     * Translates this point by the specified offsets and returns a new point.
     * 
     * Translation is a fundamental geometric transformation that moves a point
     * by adding offset values to each coordinate. This operation is essential
     * for spatial data structures and geometric algorithms that need to adjust
     * coordinate systems or perform relative positioning.
     * 
     * @param dx the offset in the x direction
     * @param dy the offset in the y direction
     * @param dz the offset in the z direction
     * @return a new Point3D representing the translated point
     * @throws IllegalArgumentException if any offset is NaN or infinite
     */
    public Point3D translate(double dx, double dy, double dz) {
        if (Double.isNaN(dx) || Double.isInfinite(dx) ||
            Double.isNaN(dy) || Double.isInfinite(dy) ||
            Double.isNaN(dz) || Double.isInfinite(dz)) {
            logger.severe("Attempted translation with invalid offsets: " +
                         "dx=" + dx + ", dy=" + dy + ", dz=" + dz);
            throw new IllegalArgumentException("Translation offsets must be finite numbers");
        }
        
        logger.info("Translating point by offset (" + dx + ", " + dy + ", " + dz + ")");
        
        return new Point3D(x + dx, y + dy, z + dz);
    }
    
    /**
     * Calculates the magnitude (length) of the vector from origin to this point.
     * 
     * The magnitude represents the distance from the origin (0,0,0) to this point,
     * calculated using the Euclidean distance formula. This is fundamental for
     * normalizing vectors and calculating vector-based algorithms.
     * 
     * @return the magnitude of the position vector
     */
    public double magnitude() {
        double mag = Math.sqrt(x * x + y * y + z * z);
        logger.info("Calculated magnitude: " + mag);
        return mag;
    }
    
    /**
     * Returns a normalized version of this point (unit vector in the same direction).
     * 
     * Normalization creates a vector with magnitude 1 pointing in the same direction.
     * This is crucial for direction calculations, dot product operations, and
     * algorithms requiring unit vectors such as lighting calculations and
     * collision detection.
     * 
     * @return a new Point3D representing the normalized vector
     * @throws IllegalStateException if this point is at the origin (cannot normalize zero vector)
     */
    public Point3D normalize() {
        double mag = magnitude();
        
        if (Math.abs(mag) < EPSILON) {
            logger.severe("Attempted to normalize zero vector");
            throw new IllegalStateException("Cannot normalize zero vector");
        }
        
        logger.info("Normalizing vector with magnitude " + mag);
        
        return new Point3D(x / mag, y / mag, z / mag);
    }
    
    /**
     * Calculates the dot product with another point (treating both as vectors).
     * 
     * The dot product is fundamental for many geometric algorithms:
     * - Angle calculations (cos θ = a·b / |a||b|)
     * - Projection operations
     * - Determining if vectors are perpendicular (dot product = 0)
     * 
     * @param other the other vector
     * @return the dot product
     * @throws IllegalArgumentException if the other point is null
     */
    public double dotProduct(Point3D other) {
        if (other == null) {
            logger.severe("Attempted dot product with null vector");
            throw new IllegalArgumentException("Cannot calculate dot product with null vector");
        }
        
        double result = this.x * other.x + this.y * other.y + this.z * other.z;
        logger.info("Calculated dot product: " + result);
        
        return result;
    }
    
    /**
     * Calculates the cross product with another point (treating both as vectors).
     * 
     * The cross product produces a vector perpendicular to both input vectors,
     * essential for calculating surface normals, determining orientation,
     * and many 3D geometric algorithms.
     * 
     * @param other the other vector
     * @return a new Point3D representing the cross product
     * @throws IllegalArgumentException if the other point is null
     */
    public Point3D crossProduct(Point3D other) {
        if (other == null) {
            logger.severe("Attempted cross product with null vector");
            throw new IllegalArgumentException("Cannot calculate cross product with null vector");
        }
        
        double newX = this.y * other.z - this.z * other.y;
        double newY = this.z * other.x - this.x * other.z;
        double newZ = this.x * other.y - this.y * other.x;
        
        logger.info("Calculated cross product");
        
        return new Point3D(newX, newY, newZ);
    }
    
    /**
     * Checks if this point is equal to another point within a small epsilon tolerance.
     * 
     * This method addresses floating-point precision issues that are crucial
     * when implementing geometric algorithms. Direct equality comparison of
     * floating-point numbers can fail due to rounding errors, making epsilon
     * comparison essential for reliable geometric computations.
     * 
     * @param other the point to compare with
     * @return true if the points are approximately equal
     */
    public boolean approximatelyEquals(Point3D other) {
        if (other == null) {
            return false;
        }
        
        boolean isEqual = Math.abs(this.x - other.x) < EPSILON &&
                         Math.abs(this.y - other.y) < EPSILON &&
                         Math.abs(this.z - other.z) < EPSILON;
        
        if (isEqual) {
            logger.info("Points are approximately equal");
        } else {
            logger.info("Points are not approximately equal");
        }
        
        return isEqual;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Point3D point3D = (Point3D) obj;
        return Double.compare(point3D.x, x) == 0 &&
               Double.compare(point3D.y, y) == 0 &&
               Double.compare(point3D.z, z) == 0;
    }
    
    @Override
    public int hashCode() {
        return Objects.hash(x, y, z);
    }
    
    @Override
    public String toString() {
        return String.format("Point3D(%.6f, %.6f, %.6f)", x, y, z);
    }
}
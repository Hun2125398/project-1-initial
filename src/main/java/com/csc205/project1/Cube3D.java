package com.csc205.project1;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Cube3D represents a cube in three-dimensional space with advanced geometric operations
 * suitable for 3D graphics applications, game engines, and computational geometry.
 * 
 * This class demonstrates sophisticated object-oriented design patterns and geometric algorithms:
 * 
 * 1. COMPOSITE PATTERN: Cube3D is composed of Point3D vertices and Line3D edges,
 *    demonstrating how complex 3D objects can be built from simpler geometric primitives.
 *    This hierarchical composition enables modular geometric operations.
 * 
 * 2. IMMUTABLE VALUE OBJECT: The cube maintains immutability after construction,
 *    ensuring thread safety and preventing accidental modifications during complex
 *    3D transformations and rendering operations.
 * 
 * 3. BUILDER PATTERN: Multiple static factory methods provide different construction
 *    strategies (axis-aligned, centered, from bounds) while maintaining clear semantics
 *    and type safety.
 * 
 * 4. TEMPLATE METHOD PATTERN: Transformation operations (rotation, translation, scaling)
 *    follow consistent patterns of validation → computation → reconstruction,
 *    enabling reliable chaining of geometric operations.
 * 
 * 5. FACADE PATTERN: Complex 3D operations are exposed through simple interfaces,
 *    hiding the underlying mathematical complexity while providing powerful functionality
 *    for graphics programming.
 * 
 * 6. FLYWEIGHT PATTERN: Edge and face calculations are computed on-demand rather than
 *    stored, reducing memory overhead while maintaining performance for graphics applications.
 * 
 * FOUNDATIONAL PRINCIPLES FOR DATA STRUCTURES AND ALGORITHMS:
 * 
 * - GEOMETRIC TRANSFORMATIONS: Implements fundamental 3D transformation matrices
 *   and operations essential for computer graphics, robotics, and spatial algorithms.
 * 
 * - SPATIAL COHERENCE: Provides bounding box operations and spatial queries that
 *   are foundational for spatial data structures like octrees, BSP trees, and
 *   collision detection systems.
 * 
 * - TOPOLOGICAL CONSISTENCY: Maintains proper vertex ordering and face orientation
 *   crucial for rendering algorithms, mesh processing, and geometric modeling.
 * 
 * - NUMERICAL ROBUSTNESS: Careful handling of floating-point precision and
 *   transformation composition to prevent accumulation of numerical errors in
 *   graphics pipelines.
 * 
 * @author Generated Code
 * @version 1.0
 */
public class Cube3D {
    
    private static final Logger logger = Logger.getLogger(Cube3D.class.getName());
    
    // Epsilon for floating-point comparisons
    private static final double EPSILON = 1e-9;
    
    // Cube vertices in a specific order for consistent topology
    // Order: bottom face (CCW from -Z), then top face (CCW from +Z)
    private final Point3D[] vertices;
    private final Point3D center;
    private final double sideLength;
    
    /**
     * Creates a cube from 8 vertices in the standard vertex ordering.
     * 
     * This constructor expects vertices in a specific order that maintains
     * topological consistency for graphics operations:
     * - Vertices 0-3: bottom face (z = min), counter-clockwise when viewed from -Z
     * - Vertices 4-7: top face (z = max), counter-clockwise when viewed from +Z
     * 
     * This ordering is crucial for proper face normal calculations, texture mapping,
     * and rendering operations in 3D graphics pipelines.
     * 
     * @param vertices array of 8 Point3D objects representing cube vertices
     * @throws IllegalArgumentException if vertices array is invalid
     */
    public Cube3D(Point3D[] vertices) {
        if (vertices == null || vertices.length != 8) {
            logger.severe("Invalid vertices array: expected 8 vertices, got " + 
                         (vertices == null ? "null" : vertices.length));
            throw new IllegalArgumentException("Cube requires exactly 8 vertices");
        }
        
        for (int i = 0; i < 8; i++) {
            if (vertices[i] == null) {
                logger.severe("Null vertex at index " + i);
                throw new IllegalArgumentException("Vertex " + i + " cannot be null");
            }
        }
        
        // Defensive copying to ensure immutability
        this.vertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            this.vertices[i] = Point3D.from(vertices[i]);
        }
        
        this.center = calculateCenter();
        this.sideLength = calculateSideLength();
        
        if (!isValidCube()) {
            logger.warning("Created cube with non-uniform edge lengths or invalid geometry");
        }
        
        logger.info("Created Cube3D with center " + center + " and side length " + 
                   String.format("%.6f", sideLength));
    }
    
    /**
     * Creates an axis-aligned cube with the specified center point and side length.
     * 
     * This factory method is the most common way to create cubes for graphics applications.
     * It generates an axis-aligned cube (sides parallel to coordinate axes) which is
     * optimal for spatial partitioning algorithms, collision detection, and rendering
     * optimizations.
     * 
     * @param center the center point of the cube
     * @param sideLength the length of each side
     * @return a new axis-aligned Cube3D
     * @throws IllegalArgumentException if center is null or side length is invalid
     */
    public static Cube3D createAxisAligned(Point3D center, double sideLength) {
        if (center == null) {
            logger.severe("Attempted to create cube with null center");
            throw new IllegalArgumentException("Center point cannot be null");
        }
        
        if (sideLength <= 0 || Double.isNaN(sideLength) || Double.isInfinite(sideLength)) {
            logger.severe("Invalid side length: " + sideLength);
            throw new IllegalArgumentException("Side length must be a positive finite number");
        }
        
        double halfSide = sideLength / 2.0;
        Point3D[] vertices = new Point3D[8];
        
        // Bottom face vertices (z = center.z - halfSide)
        vertices[0] = new Point3D(center.getX() - halfSide, center.getY() - halfSide, center.getZ() - halfSide);
        vertices[1] = new Point3D(center.getX() + halfSide, center.getY() - halfSide, center.getZ() - halfSide);
        vertices[2] = new Point3D(center.getX() + halfSide, center.getY() + halfSide, center.getZ() - halfSide);
        vertices[3] = new Point3D(center.getX() - halfSide, center.getY() + halfSide, center.getZ() - halfSide);
        
        // Top face vertices (z = center.z + halfSide)
        vertices[4] = new Point3D(center.getX() - halfSide, center.getY() - halfSide, center.getZ() + halfSide);
        vertices[5] = new Point3D(center.getX() + halfSide, center.getY() - halfSide, center.getZ() + halfSide);
        vertices[6] = new Point3D(center.getX() + halfSide, center.getY() + halfSide, center.getZ() + halfSide);
        vertices[7] = new Point3D(center.getX() - halfSide, center.getY() + halfSide, center.getZ() + halfSide);
        
        logger.info("Created axis-aligned cube with center " + center + " and side length " + sideLength);
        return new Cube3D(vertices);
    }
    
    /**
     * Creates a unit cube centered at the origin with side length 1.0.
     * 
     * Unit cubes are frequently used as reference objects in 3D graphics,
     * providing a standard primitive for testing, development, and as a basis
     * for scaled transformations. This is particularly useful in shader development
     * and procedural geometry generation.
     * 
     * @return a new unit Cube3D at the origin
     */
    public static Cube3D createUnitCube() {
        logger.info("Creating unit cube at origin");
        return createAxisAligned(Point3D.origin(), 1.0);
    }
    
    /**
     * Creates a cube that encompasses the given bounding box defined by min and max points.
     * 
     * This factory method is essential for spatial algorithms that need to create
     * bounding volumes around sets of objects. It's commonly used in octree construction,
     * frustum culling, and collision detection broad-phase algorithms.
     * 
     * @param min the minimum corner of the bounding box
     * @param max the maximum corner of the bounding box
     * @return a new Cube3D that encompasses the bounding box
     * @throws IllegalArgumentException if min or max are null, or if max is not greater than min
     */
    public static Cube3D fromBounds(Point3D min, Point3D max) {
        if (min == null || max == null) {
            logger.severe("Attempted to create cube from bounds with null points");
            throw new IllegalArgumentException("Min and max points cannot be null");
        }
        
        if (max.getX() <= min.getX() || max.getY() <= min.getY() || max.getZ() <= min.getZ()) {
            logger.severe("Invalid bounds: max must be greater than min in all dimensions");
            throw new IllegalArgumentException("Max point must be greater than min point in all dimensions");
        }
        
        // Calculate the maximum dimension to create a cube
        double xSize = max.getX() - min.getX();
        double ySize = max.getY() - min.getY();
        double zSize = max.getZ() - min.getZ();
        double maxSize = Math.max(Math.max(xSize, ySize), zSize);
        
        Point3D center = new Point3D(
            (min.getX() + max.getX()) / 2.0,
            (min.getY() + max.getY()) / 2.0,
            (min.getZ() + max.getZ()) / 2.0
        );
        
        logger.info("Created cube from bounds with center " + center + " and side length " + maxSize);
        return createAxisAligned(center, maxSize);
    }
    
    /**
     * Returns a defensive copy of the vertices array.
     * 
     * This method provides access to the cube's vertices while maintaining
     * immutability. The returned array is a copy, preventing external modification
     * of the cube's internal state.
     * 
     * @return a copy of the vertices array
     */
    public Point3D[] getVertices() {
        Point3D[] copy = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            copy[i] = Point3D.from(vertices[i]);
        }
        return copy;
    }
    
    /**
     * Returns the center point of this cube.
     * 
     * The center is calculated as the average of all vertices, providing
     * the geometric centroid which is essential for rotation operations
     * and spatial queries.
     * 
     * @return the center point of the cube
     */
    public Point3D getCenter() {
        return Point3D.from(center);
    }
    
    /**
     * Returns the side length of this cube.
     * 
     * For non-axis-aligned or deformed cubes, this returns the average
     * edge length. For perfect cubes, all edges should have the same length.
     * 
     * @return the side length of the cube
     */
    public double getSideLength() {
        return sideLength;
    }
    
    /**
     * Calculates and returns the volume of this cube.
     * 
     * Volume calculation is fundamental for physics simulations, spatial
     * partitioning algorithms, and resource management in graphics applications.
     * For a perfect cube, volume = sideLength³.
     * 
     * @return the volume of the cube
     */
    public double getVolume() {
        double volume = Math.pow(sideLength, 3);
        logger.info("Calculated cube volume: " + volume);
        return volume;
    }
    
    /**
     * Calculates and returns the surface area of this cube.
     * 
     * Surface area is important for texture mapping calculations, heat transfer
     * simulations, and rendering optimizations. For a perfect cube,
     * surface area = 6 × sideLength².
     * 
     * @return the surface area of the cube
     */
    public double getSurfaceArea() {
        double area = 6 * Math.pow(sideLength, 2);
        logger.info("Calculated cube surface area: " + area);
        return area;
    }
    
    /**
     * Calculates and returns the perimeter length of this cube.
     * 
     * The perimeter represents the total length of all edges. A cube has 12 edges,
     * so for a perfect cube, perimeter = 12 × sideLength. This is useful for
     * wireframe rendering and geometric analysis.
     * 
     * @return the total perimeter length of all edges
     */
    public double getPerimeter() {
        double totalLength = 0;
        List<Line3D> edges = getEdges();
        
        for (Line3D edge : edges) {
            totalLength += edge.length();
        }
        
        logger.info("Calculated cube perimeter: " + totalLength);
        return totalLength;
    }
    
    /**
     * Returns all 12 edges of the cube as Line3D objects.
     * 
     * The edges are organized as:
     * - 4 edges for the bottom face
     * - 4 edges for the top face  
     * - 4 vertical edges connecting bottom and top faces
     * 
     * This representation is essential for wireframe rendering, collision detection,
     * and geometric analysis algorithms.
     * 
     * @return a list of 12 Line3D objects representing the cube's edges
     */
    public List<Line3D> getEdges() {
        List<Line3D> edges = new ArrayList<>();
        
        // Bottom face edges (0-1-2-3-0)
        edges.add(new Line3D(vertices[0], vertices[1]));
        edges.add(new Line3D(vertices[1], vertices[2]));
        edges.add(new Line3D(vertices[2], vertices[3]));
        edges.add(new Line3D(vertices[3], vertices[0]));
        
        // Top face edges (4-5-6-7-4)
        edges.add(new Line3D(vertices[4], vertices[5]));
        edges.add(new Line3D(vertices[5], vertices[6]));
        edges.add(new Line3D(vertices[6], vertices[7]));
        edges.add(new Line3D(vertices[7], vertices[4]));
        
        // Vertical edges connecting bottom and top faces
        edges.add(new Line3D(vertices[0], vertices[4]));
        edges.add(new Line3D(vertices[1], vertices[5]));
        edges.add(new Line3D(vertices[2], vertices[6]));
        edges.add(new Line3D(vertices[3], vertices[7]));
        
        logger.info("Generated 12 cube edges");
        return edges;
    }
    
    /**
     * Returns the 6 faces of the cube, each represented as an array of 4 vertices.
     * 
     * Faces are ordered consistently for proper normal calculation and rendering:
     * - Bottom face (indices 0,1,2,3)
     * - Top face (indices 4,5,6,7)
     * - Front, Right, Back, Left faces in order
     * 
     * This is fundamental for mesh generation, collision detection, and rendering
     * operations that require face-based processing.
     * 
     * @return a list of 6 arrays, each containing 4 vertices representing a face
     */
    public List<Point3D[]> getFaces() {
        List<Point3D[]> faces = new ArrayList<>();
        
        // Bottom face (z = min)
        faces.add(new Point3D[]{vertices[0], vertices[1], vertices[2], vertices[3]});
        
        // Top face (z = max)
        faces.add(new Point3D[]{vertices[4], vertices[5], vertices[6], vertices[7]});
        
        // Front face (y = min)
        faces.add(new Point3D[]{vertices[0], vertices[1], vertices[5], vertices[4]});
        
        // Right face (x = max)
        faces.add(new Point3D[]{vertices[1], vertices[2], vertices[6], vertices[5]});
        
        // Back face (y = max)
        faces.add(new Point3D[]{vertices[2], vertices[3], vertices[7], vertices[6]});
        
        // Left face (x = min)
        faces.add(new Point3D[]{vertices[3], vertices[0], vertices[4], vertices[7]});
        
        logger.info("Generated 6 cube faces");
        return faces;
    }
    
    /**
     * Calculates the face normal for a given face index.
     * 
     * Face normals are essential for lighting calculations, collision detection,
     * and determining face orientation in 3D graphics. The normal is calculated
     * using the cross product of two edge vectors of the face.
     * 
     * @param faceIndex the index of the face (0-5)
     * @return the normalized face normal vector
     * @throws IllegalArgumentException if face index is invalid
     */
    public Point3D getFaceNormal(int faceIndex) {
        if (faceIndex < 0 || faceIndex >= 6) {
            logger.severe("Invalid face index: " + faceIndex);
            throw new IllegalArgumentException("Face index must be between 0 and 5");
        }
        
        List<Point3D[]> faces = getFaces();
        Point3D[] face = faces.get(faceIndex);
        
        // Calculate normal using cross product of two edges
        Point3D edge1 = face[1].translate(-face[0].getX(), -face[0].getY(), -face[0].getZ());
        Point3D edge2 = face[3].translate(-face[0].getX(), -face[0].getY(), -face[0].getZ());
        
        Point3D normal = edge1.crossProduct(edge2).normalize();
        
        logger.info("Calculated normal for face " + faceIndex + ": " + normal);
        return normal;
    }
    
    /**
     * Translates the cube by the specified offset vector.
     * 
     * Translation is a fundamental 3D transformation that moves the entire cube
     * without changing its orientation or size. This operation is essential for
     * positioning objects in 3D space and implementing camera movements.
     * 
     * @param offset the translation vector
     * @return a new translated Cube3D
     * @throws IllegalArgumentException if offset is null
     */
    public Cube3D translate(Point3D offset) {
        if (offset == null) {
            logger.severe("Attempted to translate cube with null offset");
            throw new IllegalArgumentException("Offset vector cannot be null");
        }
        
        Point3D[] newVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            newVertices[i] = vertices[i].translate(offset.getX(), offset.getY(), offset.getZ());
        }
        
        logger.info("Translated cube by offset " + offset);
        return new Cube3D(newVertices);
    }
    
    /**
     * Rotates the cube around the X-axis by the specified angle.
     * 
     * Rotation around the X-axis is one of the fundamental 3D transformations.
     * The rotation is performed around the cube's center point to maintain
     * spatial coherence. This is essential for object orientation in 3D graphics.
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new rotated Cube3D
     */
    public Cube3D rotateX(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ", returning original cube");
            return new Cube3D(this.vertices);
        }
        
        Point3D[] newVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            // Translate to origin, rotate, then translate back
            Point3D translated = vertices[i].translate(-center.getX(), -center.getY(), -center.getZ());
            Point3D rotated = translated.rotateX(angleRadians);
            newVertices[i] = rotated.translate(center.getX(), center.getY(), center.getZ());
        }
        
        logger.info("Rotated cube around X-axis by " + Math.toDegrees(angleRadians) + " degrees");
        return new Cube3D(newVertices);
    }
    
    /**
     * Rotates the cube around the Y-axis by the specified angle.
     * 
     * Y-axis rotation is commonly used for horizontal rotation in 3D applications,
     * such as rotating objects around their vertical axis or implementing turntable
     * camera movements.
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new rotated Cube3D
     */
    public Cube3D rotateY(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ", returning original cube");
            return new Cube3D(this.vertices);
        }
        
        Point3D[] newVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            Point3D translated = vertices[i].translate(-center.getX(), -center.getY(), -center.getZ());
            Point3D rotated = translated.rotateY(angleRadians);
            newVertices[i] = rotated.translate(center.getX(), center.getY(), center.getZ());
        }
        
        logger.info("Rotated cube around Y-axis by " + Math.toDegrees(angleRadians) + " degrees");
        return new Cube3D(newVertices);
    }
    
    /**
     * Rotates the cube around the Z-axis by the specified angle.
     * 
     * Z-axis rotation affects the X and Y coordinates while leaving Z unchanged.
     * This is commonly used for 2D-style rotation within a 3D context and for
     * implementing roll movements in flight simulations.
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new rotated Cube3D
     */
    public Cube3D rotateZ(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.warning("Invalid rotation angle: " + angleRadians + ", returning original cube");
            return new Cube3D(this.vertices);
        }
        
        Point3D[] newVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            Point3D translated = vertices[i].translate(-center.getX(), -center.getY(), -center.getZ());
            Point3D rotated = translated.rotateZ(angleRadians);
            newVertices[i] = rotated.translate(center.getX(), center.getY(), center.getZ());
        }
        
        logger.info("Rotated cube around Z-axis by " + Math.toDegrees(angleRadians) + " degrees");
        return new Cube3D(newVertices);
    }
    
    /**
     * Scales the cube by the specified factor around its center.
     * 
     * Uniform scaling changes the size of the cube while maintaining its proportions
     * and center position. This is fundamental for level-of-detail algorithms,
     * animation systems, and procedural geometry generation.
     * 
     * @param scaleFactor the scaling factor (1.0 = no change)
     * @return a new scaled Cube3D
     * @throws IllegalArgumentException if scale factor is invalid
     */
    public Cube3D scale(double scaleFactor) {
        if (scaleFactor <= 0 || Double.isNaN(scaleFactor) || Double.isInfinite(scaleFactor)) {
            logger.severe("Invalid scale factor: " + scaleFactor);
            throw new IllegalArgumentException("Scale factor must be a positive finite number");
        }
        
        Point3D[] newVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            // Scale around center point
            Point3D translated = vertices[i].translate(-center.getX(), -center.getY(), -center.getZ());
            Point3D scaled = new Point3D(
                translated.getX() * scaleFactor,
                translated.getY() * scaleFactor,
                translated.getZ() * scaleFactor
            );
            newVertices[i] = scaled.translate(center.getX(), center.getY(), center.getZ());
        }
        
        logger.info("Scaled cube by factor " + scaleFactor);
        return new Cube3D(newVertices);
    }
    
    /**
     * Returns the axis-aligned bounding box of this cube.
     * 
     * The AABB is fundamental for spatial algorithms, collision detection,
     * and rendering optimizations. Even for rotated cubes, the AABB provides
     * a quick spatial approximation for broad-phase algorithms.
     * 
     * @return an array containing [minPoint, maxPoint] of the bounding box
     */
    public Point3D[] getAxisAlignedBoundingBox() {
        double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE, minZ = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE, maxY = Double.MIN_VALUE, maxZ = Double.MIN_VALUE;
        
        for (Point3D vertex : vertices) {
            minX = Math.min(minX, vertex.getX());
            minY = Math.min(minY, vertex.getY());
            minZ = Math.min(minZ, vertex.getZ());
            maxX = Math.max(maxX, vertex.getX());
            maxY = Math.max(maxY, vertex.getY());
            maxZ = Math.max(maxZ, vertex.getZ());
        }
        
        Point3D min = new Point3D(minX, minY, minZ);
        Point3D max = new Point3D(maxX, maxY, maxZ);
        
        logger.info("Calculated AABB: min=" + min + ", max=" + max);
        return new Point3D[]{min, max};
    }
    
    /**
     * Checks if this cube intersects with another cube using AABB intersection.
     * 
     * AABB intersection is a fast broad-phase collision detection algorithm
     * that's essential for spatial partitioning and collision systems.
     * For more precise collision detection, additional algorithms would be needed.
     * 
     * @param other the other cube to check intersection with
     * @return true if the cubes' bounding boxes intersect
     * @throws IllegalArgumentException if other cube is null
     */
    public boolean intersects(Cube3D other) {
        if (other == null) {
            logger.severe("Attempted intersection test with null cube");
            throw new IllegalArgumentException("Other cube cannot be null");
        }
        
        Point3D[] thisBounds = this.getAxisAlignedBoundingBox();
        Point3D[] otherBounds = other.getAxisAlignedBoundingBox();
        
        boolean intersects = 
            thisBounds[0].getX() <= otherBounds[1].getX() && thisBounds[1].getX() >= otherBounds[0].getX() &&
            thisBounds[0].getY() <= otherBounds[1].getY() && thisBounds[1].getY() >= otherBounds[0].getY() &&
            thisBounds[0].getZ() <= otherBounds[1].getZ() && thisBounds[1].getZ() >= otherBounds[0].getZ();
        
        logger.info("Cube intersection test result: " + intersects);
        return intersects;
    }
    
    /**
     * Checks if a point is inside this cube (within its bounding box).
     * 
     * Point-in-cube testing is fundamental for spatial queries, collision detection,
     * and selection algorithms in 3D applications. This implementation uses
     * AABB testing for efficiency.
     * 
     * @param point the point to test
     * @return true if the point is inside the cube's bounding box
     * @throws IllegalArgumentException if point is null
     */
    public boolean contains(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to test containment of null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        Point3D[] bounds = getAxisAlignedBoundingBox();
        
        boolean contains = 
            point.getX() >= bounds[0].getX() && point.getX() <= bounds[1].getX() &&
            point.getY() >= bounds[0].getY() && point.getY() <= bounds[1].getY() &&
            point.getZ() >= bounds[0].getZ() && point.getZ() <= bounds[1].getZ();
        
        logger.info("Point containment test for " + point + ": " + contains);
        return contains;
    }
    
    /**
     * Calculates the distance from a point to the surface of this cube.
     * 
     * Distance calculations are essential for proximity queries, force field
     * computations, and collision response algorithms. This implementation
     * calculates the distance to the cube's bounding box.
     * 
     * @param point the point to calculate distance from
     * @return the minimum distance from the point to the cube's surface
     * @throws IllegalArgumentException if point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        if (contains(point)) {
            logger.info("Point is inside cube, distance = 0");
            return 0.0;
        }
        
        Point3D[] bounds = getAxisAlignedBoundingBox();
        
        // Calculate distance to AABB
        double dx = Math.max(0, Math.max(bounds[0].getX() - point.getX(), 
                                        point.getX() - bounds[1].getX()));
        double dy = Math.max(0, Math.max(bounds[0].getY() - point.getY(), 
                                        point.getY() - bounds[1].getY()));
        double dz = Math.max(0, Math.max(bounds[0].getZ() - point.getZ(), 
                                        point.getZ() - bounds[1].getZ()));
        
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        logger.info("Distance from point " + point + " to cube: " + distance);
        return distance;
    }
    
    /**
     * Checks if this cube is axis-aligned (sides parallel to coordinate axes).
     * 
     * Axis-aligned cubes have optimized algorithms for many operations including
     * collision detection, spatial partitioning, and rendering. This check is
     * useful for selecting appropriate algorithms based on cube orientation.
     * 
     * @return true if the cube is axis-aligned
     */
    public boolean isAxisAligned() {
        // Check if all edges are parallel to coordinate axes
        List<Line3D> edges = getEdges();
        
        for (Line3D edge : edges) {
            Point3D direction = edge.getDirection();
            
            // Count how many components are effectively zero
            int zeroComponents = 0;
            if (Math.abs(direction.getX()) < EPSILON) zeroComponents++;
            if (Math.abs(direction.getY()) < EPSILON) zeroComponents++;
            if (Math.abs(direction.getZ()) < EPSILON) zeroComponents++;
            
            // For axis-aligned cube, each edge should be parallel to exactly one axis
            if (zeroComponents != 2) {
                logger.info("Cube is not axis-aligned");
                return false;
            }
        }
        
        logger.info("Cube is axis-aligned");
        return true;
    }
    
    /**
     * Returns the closest vertex to a given point.
     * 
     * Finding the closest vertex is useful for snapping operations, vertex selection
     * in modeling applications, and proximity-based algorithms. This operation is
     * fundamental for user interaction in 3D editors.
     * 
     * @param point the reference point
     * @return the closest vertex to the given point
     * @throws IllegalArgumentException if point is null
     */
    public Point3D getClosestVertex(Point3D point) {
        if (point == null) {
            logger.severe("Attempted to find closest vertex to null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        Point3D closest = vertices[0];
        double minDistance = point.distanceTo(vertices[0]);
        
        for (int i = 1; i < 8; i++) {
            double distance = point.distanceTo(vertices[i]);
            if (distance < minDistance) {
                minDistance = distance;
                closest = vertices[i];
            }
        }
        
        logger.info("Closest vertex to " + point + " is " + closest + 
                   " at distance " + minDistance);
        return Point3D.from(closest);
    }
    
    /**
     * Calculates the diagonal length of the cube.
     * 
     * The space diagonal connects two opposite corners of the cube and represents
     * the maximum distance between any two points within the cube. For a perfect
     * cube with side length s, the diagonal length is s√3.
     * 
     * @return the length of the space diagonal
     */
    public double getDiagonalLength() {
        // Find the maximum distance between any two vertices
        double maxDistance = 0;
        
        for (int i = 0; i < 8; i++) {
            for (int j = i + 1; j < 8; j++) {
                double distance = vertices[i].distanceTo(vertices[j]);
                if (distance > maxDistance) {
                    maxDistance = distance;
                }
            }
        }
        
        logger.info("Calculated cube diagonal length: " + maxDistance);
        return maxDistance;
    }
    
    /**
     * Creates a wireframe representation of the cube for debug visualization.
     * 
     * Wireframe representation is essential for debug rendering, CAD applications,
     * and situations where you need to visualize the cube's structure without
     * solid faces. This method returns a formatted string suitable for logging
     * or text-based visualization.
     * 
     * @return a string representation of the cube's wireframe
     */
    public String getWireframeRepresentation() {
        StringBuilder wireframe = new StringBuilder();
        wireframe.append("Cube3D Wireframe:\n");
        wireframe.append("Center: ").append(center).append("\n");
        wireframe.append("Side Length: ").append(String.format("%.6f", sideLength)).append("\n");
        wireframe.append("Vertices:\n");
        
        for (int i = 0; i < 8; i++) {
            wireframe.append("  V").append(i).append(": ").append(vertices[i]).append("\n");
        }
        
        wireframe.append("Edges:\n");
        List<Line3D> edges = getEdges();
        for (int i = 0; i < edges.size(); i++) {
            Line3D edge = edges.get(i);
            wireframe.append("  E").append(i).append(": ")
                     .append(edge.getStart()).append(" -> ").append(edge.getEnd())
                     .append(" (length: ").append(String.format("%.6f", edge.length())).append(")\n");
        }
        
        logger.info("Generated wireframe representation");
        return wireframe.toString();
    }
    
    /**
     * Private helper method to calculate the center point from vertices.
     * 
     * The center is calculated as the arithmetic mean of all vertices,
     * providing the geometric centroid of the cube.
     * 
     * @return the calculated center point
     */
    private Point3D calculateCenter() {
        double sumX = 0, sumY = 0, sumZ = 0;
        
        for (Point3D vertex : vertices) {
            sumX += vertex.getX();
            sumY += vertex.getY();
            sumZ += vertex.getZ();
        }
        
        return new Point3D(sumX / 8.0, sumY / 8.0, sumZ / 8.0);
    }
    
    /**
     * Private helper method to calculate the average side length.
     * 
     * For perfect cubes, all edges should have the same length.
     * This method calculates the average edge length to handle
     * slightly deformed cubes due to floating-point precision.
     * 
     * @return the calculated average side length
     */
    private double calculateSideLength() {
        List<Line3D> edges = getEdges();
        double totalLength = 0;
        
        for (Line3D edge : edges) {
            totalLength += edge.length();
        }
        
        return totalLength / edges.size();
    }
    
    /**
     * Private helper method to validate cube geometry.
     * 
     * This method performs basic validation to ensure the cube has
     * reasonable geometry, checking for consistent edge lengths
     * and proper vertex relationships.
     * 
     * @return true if the cube has valid geometry
     */
    private boolean isValidCube() {
        List<Line3D> edges = getEdges();
        
        // Check if all edges have similar lengths (within tolerance)
        double firstLength = edges.get(0).length();
        
        for (Line3D edge : edges) {
            if (Math.abs(edge.length() - firstLength) > EPSILON * firstLength) {
                return false;
            }
        }
        
        return true;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Cube3D cube3D = (Cube3D) obj;
        return Arrays.equals(vertices, cube3D.vertices);
    }
    
    @Override
    public int hashCode() {
        return Arrays.hashCode(vertices);
    }
    
    @Override
    public String toString() {
        return String.format("Cube3D[center=%s, sideLength=%.6f, volume=%.6f]", 
                           center, sideLength, getVolume());
    }
}
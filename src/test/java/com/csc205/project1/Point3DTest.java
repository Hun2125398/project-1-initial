package com.csc205.project1;

import org.junit.jupiter.api.*;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.params.provider.CsvSource;
import static org.junit.jupiter.api.Assertions.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.logging.Handler;
import java.util.logging.LogRecord;
import java.util.ArrayList;
import java.util.List;

/**
 * Comprehensive test suite for Point3D class demonstrating advanced testing patterns
 * and geometric algorithm validation.
 * 
 * This test class demonstrates several testing design patterns and principles:
 * 
 * 1. ARRANGE-ACT-ASSERT PATTERN: Each test method follows the AAA pattern for clarity
 *    and maintainability, making tests easy to understand and modify.
 * 
 * 2. BOUNDARY VALUE TESTING: Tests edge cases including zero values, negative values,
 *    very large values, and floating-point precision limits to ensure robustness.
 * 
 * 3. EQUIVALENCE PARTITIONING: Groups test cases into equivalence classes to ensure
 *    comprehensive coverage without redundant testing.
 * 
 * 4. PARAMETERIZED TESTING: Uses JUnit 5 parameterized tests to verify behavior
 *    across multiple input values efficiently.
 * 
 * 5. EXCEPTION TESTING: Validates that appropriate exceptions are thrown for invalid
 *    inputs, ensuring defensive programming principles are maintained.
 * 
 * 6. MATHEMATICAL PROPERTY TESTING: Verifies geometric and mathematical properties
 *    such as distance symmetry, vector operations, and transformation invariants.
 * 
 * 7. LOGGING VERIFICATION: Tests that appropriate log messages are generated,
 *    ensuring observability and debugging support in production systems.
 * 
 * TESTING PRINCIPLES FOR GEOMETRIC ALGORITHMS:
 * 
 * - FLOATING-POINT PRECISION: Uses epsilon-based comparisons for reliable testing
 *   of floating-point calculations in geometric operations.
 * 
 * - GEOMETRIC INVARIANTS: Verifies mathematical properties that must hold true
 *   for geometric operations (e.g., distance symmetry, vector magnitude properties).
 * 
 * - TRANSFORMATION CONSISTENCY: Ensures that geometric transformations preserve
 *   expected properties and relationships between points.
 * 
 * - PERFORMANCE CONSIDERATIONS: Includes tests for optimized operations like
 *   squared distance calculations to ensure they provide correct results.
 * 
 * @author Generated Code
 * @version 1.0
 */
@DisplayName("Point3D Test Suite")
class Point3DTest {
    
    private static final double EPSILON = 1e-9;
    private static final double DELTA = 1e-6; // For assertEquals with floating-point
    
    // Test fixture points for consistent testing
    private Point3D origin;
    private Point3D unitX;
    private Point3D unitY;
    private Point3D unitZ;
    private Point3D positivePoint;
    private Point3D negativePoint;
    
    // Custom log handler for testing log output
    private TestLogHandler logHandler;
    private Logger logger;
    
    /**
     * Custom log handler for capturing and verifying log messages during tests.
     * This demonstrates how to test logging behavior in unit tests.
     */
    private static class TestLogHandler extends Handler {
        private final List<LogRecord> logRecords = new ArrayList<>();
        
        @Override
        public void publish(LogRecord record) {
            logRecords.add(record);
        }
        
        @Override
        public void flush() {}
        
        @Override
        public void close() throws SecurityException {}
        
        public List<LogRecord> getLogRecords() {
            return new ArrayList<>(logRecords);
        }
        
        public void clear() {
            logRecords.clear();
        }
        
        public boolean containsMessage(String message) {
            return logRecords.stream()
                    .anyMatch(record -> record.getMessage().contains(message));
        }
        
        public boolean containsLevel(Level level) {
            return logRecords.stream()
                    .anyMatch(record -> record.getLevel().equals(level));
        }
    }
    
    @BeforeEach
    void setUp() {
        // Initialize test fixture points
        origin = Point3D.origin();
        unitX = new Point3D(1.0, 0.0, 0.0);
        unitY = new Point3D(0.0, 1.0, 0.0);
        unitZ = new Point3D(0.0, 0.0, 1.0);
        positivePoint = new Point3D(3.0, 4.0, 5.0);
        negativePoint = new Point3D(-2.0, -3.0, -1.0);
        
        // Set up logging test infrastructure
        logger = Logger.getLogger(Point3D.class.getName());
        logHandler = new TestLogHandler();
        logger.addHandler(logHandler);
        logger.setLevel(Level.ALL);
    }
    
    @AfterEach
    void tearDown() {
        // Clean up logging infrastructure
        logger.removeHandler(logHandler);
    }
    
    // ============================================================================
    // CONSTRUCTOR TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Constructor should create point with valid coordinates")
    void testValidConstruction() {
        // Arrange
        double x = 1.5, y = -2.7, z = 3.14;
        
        // Act
        Point3D point = new Point3D(x, y, z);
        
        // Assert
        assertEquals(x, point.getX(), DELTA);
        assertEquals(y, point.getY(), DELTA);
        assertEquals(z, point.getZ(), DELTA);
        assertTrue(logHandler.containsMessage("Created new Point3D"));
        assertTrue(logHandler.containsLevel(Level.INFO));
    }
    
    @ParameterizedTest
    @DisplayName("Constructor should reject NaN coordinates")
    @ValueSource(doubles = {Double.NaN})
    void testConstructorRejectsNaN(double invalidValue) {
        // Test each coordinate position with NaN
        assertAll(
            () -> assertThrows(IllegalArgumentException.class, 
                () -> new Point3D(invalidValue, 0, 0)),
            () -> assertThrows(IllegalArgumentException.class, 
                () -> new Point3D(0, invalidValue, 0)),
            () -> assertThrows(IllegalArgumentException.class, 
                () -> new Point3D(0, 0, invalidValue))
        );
        assertTrue(logHandler.containsLevel(Level.SEVERE));
    }
    
    @ParameterizedTest
    @DisplayName("Constructor should reject infinite coordinates")
    @ValueSource(doubles = {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY})
    void testConstructorRejectsInfinity(double invalidValue) {
        // Test each coordinate position with infinity
        assertAll(
            () -> assertThrows(IllegalArgumentException.class, 
                () -> new Point3D(invalidValue, 0, 0)),
            () -> assertThrows(IllegalArgumentException.class, 
                () -> new Point3D(0, invalidValue, 0)),
            () -> assertThrows(IllegalArgumentException.class, 
                () -> new Point3D(0, 0, invalidValue))
        );
        assertTrue(logHandler.containsLevel(Level.SEVERE));
    }
    
    // ============================================================================
    // FACTORY METHOD TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Origin factory method should create point at (0,0,0)")
    void testOriginFactory() {
        // Act
        Point3D originPoint = Point3D.origin();
        
        // Assert
        assertEquals(0.0, originPoint.getX(), DELTA);
        assertEquals(0.0, originPoint.getY(), DELTA);
        assertEquals(0.0, originPoint.getZ(), DELTA);
        assertTrue(logHandler.containsMessage("Creating Point3D at origin"));
    }
    
    @Test
    @DisplayName("From factory method should create copy of existing point")
    void testFromFactory() {
        // Act
        Point3D copy = Point3D.from(positivePoint);
        
        // Assert
        assertEquals(positivePoint.getX(), copy.getX(), DELTA);
        assertEquals(positivePoint.getY(), copy.getY(), DELTA);
        assertEquals(positivePoint.getZ(), copy.getZ(), DELTA);
        assertNotSame(positivePoint, copy); // Ensure it's a different object
        assertTrue(logHandler.containsMessage("Creating Point3D copy"));
    }
    
    @Test
    @DisplayName("From factory method should reject null input")
    void testFromFactoryRejectsNull() {
        // Act & Assert
        IllegalArgumentException exception = assertThrows(
            IllegalArgumentException.class, 
            () -> Point3D.from(null)
        );
        assertEquals("Cannot create point from null reference", exception.getMessage());
        assertTrue(logHandler.containsLevel(Level.SEVERE));
    }
    
    // ============================================================================
    // DISTANCE CALCULATION TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Distance calculation should be symmetric")
    void testDistanceSymmetry() {
        // Act
        double distance1 = positivePoint.distanceTo(negativePoint);
        double distance2 = negativePoint.distanceTo(positivePoint);
        
        // Assert
        assertEquals(distance1, distance2, DELTA);
        assertTrue(distance1 > 0);
    }
    
    @Test
    @DisplayName("Distance to self should be zero")
    void testDistanceToSelf() {
        // Act
        double distance = positivePoint.distanceTo(positivePoint);
        
        // Assert
        assertEquals(0.0, distance, DELTA);
    }
    
    @ParameterizedTest
    @DisplayName("Distance calculation should work for known values")
    @CsvSource({
        "0, 0, 0, 1, 0, 0, 1.0",      // Unit distance along X
        "0, 0, 0, 0, 1, 0, 1.0",      // Unit distance along Y
        "0, 0, 0, 0, 0, 1, 1.0",      // Unit distance along Z
        "0, 0, 0, 1, 1, 1, 1.732050807568877",  // Distance to (1,1,1)
        "1, 2, 3, 4, 6, 8, 5.830951894845301"   // General case
    })
    void testKnownDistances(double x1, double y1, double z1, 
                           double x2, double y2, double z2, 
                           double expectedDistance) {
        // Arrange
        Point3D point1 = new Point3D(x1, y1, z1);
        Point3D point2 = new Point3D(x2, y2, z2);
        
        // Act
        double actualDistance = point1.distanceTo(point2);
        
        // Assert
        assertEquals(expectedDistance, actualDistance, DELTA);
    }
    
    @Test
    @DisplayName("Distance to null should throw exception")
    void testDistanceToNull() {
        // Act & Assert
        IllegalArgumentException exception = assertThrows(
            IllegalArgumentException.class,
            () -> positivePoint.distanceTo(null)
        );
        assertEquals("Cannot calculate distance to null point", exception.getMessage());
        assertTrue(logHandler.containsLevel(Level.SEVERE));
    }
    
    @Test
    @DisplayName("Squared distance should equal distance squared")
    void testSquaredDistance() {
        // Act
        double distance = positivePoint.distanceTo(negativePoint);
        double squaredDistance = positivePoint.distanceSquaredTo(negativePoint);
        
        // Assert
        assertEquals(distance * distance, squaredDistance, DELTA);
    }
    
    @Test
    @DisplayName("Manhattan distance should sum absolute differences")
    void testManhattanDistance() {
        // Arrange
        Point3D point1 = new Point3D(1, 2, 3);
        Point3D point2 = new Point3D(4, 6, 8);
        double expectedManhattan = Math.abs(4-1) + Math.abs(6-2) + Math.abs(8-3);
        
        // Act
        double actualManhattan = point1.manhattanDistanceTo(point2);
        
        // Assert
        assertEquals(expectedManhattan, actualManhattan, DELTA);
    }
    
    // ============================================================================
    // ROTATION TESTS
    // ============================================================================
    
    @Test
    @DisplayName("X-axis rotation should preserve X coordinate")
    void testXRotationPreservesX() {
        // Arrange
        double angle = Math.PI / 4; // 45 degrees
        
        // Act
        Point3D rotated = positivePoint.rotateX(angle);
        
        // Assert
        assertEquals(positivePoint.getX(), rotated.getX(), DELTA);
        assertNotEquals(positivePoint.getY(), rotated.getY(), DELTA);
        assertNotEquals(positivePoint.getZ(), rotated.getZ(), DELTA);
    }
    
    @Test
    @DisplayName("Y-axis rotation should preserve Y coordinate")
    void testYRotationPreservesY() {
        // Arrange
        double angle = Math.PI / 3; // 60 degrees
        
        // Act
        Point3D rotated = positivePoint.rotateY(angle);
        
        // Assert
        assertNotEquals(positivePoint.getX(), rotated.getX(), DELTA);
        assertEquals(positivePoint.getY(), rotated.getY(), DELTA);
        assertNotEquals(positivePoint.getZ(), rotated.getZ(), DELTA);
    }
    
    @Test
    @DisplayName("Z-axis rotation should preserve Z coordinate")
    void testZRotationPreservesZ() {
        // Arrange
        double angle = Math.PI / 6; // 30 degrees
        
        // Act
        Point3D rotated = positivePoint.rotateZ(angle);
        
        // Assert
        assertNotEquals(positivePoint.getX(), rotated.getX(), DELTA);
        assertNotEquals(positivePoint.getY(), rotated.getY(), DELTA);
        assertEquals(positivePoint.getZ(), rotated.getZ(), DELTA);
    }
    
    @Test
    @DisplayName("360-degree rotation should return to original position")
    void testFullRotationReturnsToOriginal() {
        // Arrange
        double fullRotation = 2 * Math.PI;
        
        // Act
        Point3D rotatedX = positivePoint.rotateX(fullRotation);
        Point3D rotatedY = positivePoint.rotateY(fullRotation);
        Point3D rotatedZ = positivePoint.rotateZ(fullRotation);
        
        // Assert
        assertTrue(positivePoint.approximatelyEquals(rotatedX));
        assertTrue(positivePoint.approximatelyEquals(rotatedY));
        assertTrue(positivePoint.approximatelyEquals(rotatedZ));
    }
    
    @ParameterizedTest
    @DisplayName("Invalid rotation angles should be handled gracefully")
    @ValueSource(doubles = {Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY})
    void testInvalidRotationAngles(double invalidAngle) {
        // Act
        Point3D rotatedX = positivePoint.rotateX(invalidAngle);
        Point3D rotatedY = positivePoint.rotateY(invalidAngle);
        Point3D rotatedZ = positivePoint.rotateZ(invalidAngle);
        
        // Assert - should return copy of original point
        assertTrue(positivePoint.approximatelyEquals(rotatedX));
        assertTrue(positivePoint.approximatelyEquals(rotatedY));
        assertTrue(positivePoint.approximatelyEquals(rotatedZ));
        assertTrue(logHandler.containsLevel(Level.WARNING));
    }
    
    // ============================================================================
    // TRANSLATION TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Translation should add offsets to coordinates")
    void testTranslation() {
        // Arrange
        double dx = 2.5, dy = -1.7, dz = 3.8;
        
        // Act
        Point3D translated = positivePoint.translate(dx, dy, dz);
        
        // Assert
        assertEquals(positivePoint.getX() + dx, translated.getX(), DELTA);
        assertEquals(positivePoint.getY() + dy, translated.getY(), DELTA);
        assertEquals(positivePoint.getZ() + dz, translated.getZ(), DELTA);
    }
    
    @Test
    @DisplayName("Zero translation should return equivalent point")
    void testZeroTranslation() {
        // Act
        Point3D translated = positivePoint.translate(0, 0, 0);
        
        // Assert
        assertEquals(positivePoint.getX(), translated.getX(), DELTA);
        assertEquals(positivePoint.getY(), translated.getY(), DELTA);
        assertEquals(positivePoint.getZ(), translated.getZ(), DELTA);
    }
    
    @ParameterizedTest
    @DisplayName("Invalid translation offsets should throw exception")
    @ValueSource(doubles = {Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY})
    void testInvalidTranslationOffsets(double invalidOffset) {
        // Act & Assert
        assertAll(
            () -> assertThrows(IllegalArgumentException.class,
                () -> positivePoint.translate(invalidOffset, 0, 0)),
            () -> assertThrows(IllegalArgumentException.class,
                () -> positivePoint.translate(0, invalidOffset, 0)),
            () -> assertThrows(IllegalArgumentException.class,
                () -> positivePoint.translate(0, 0, invalidOffset))
        );
        assertTrue(logHandler.containsLevel(Level.SEVERE));
    }
    
    // ============================================================================
    // VECTOR OPERATION TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Magnitude should calculate vector length correctly")
    void testMagnitude() {
        // Arrange
        Point3D point = new Point3D(3, 4, 0); // 3-4-5 triangle
        
        // Act
        double magnitude = point.magnitude();
        
        // Assert
        assertEquals(5.0, magnitude, DELTA);
    }
    
    @Test
    @DisplayName("Magnitude of origin should be zero")
    void testOriginMagnitude() {
        // Act
        double magnitude = origin.magnitude();
        
        // Assert
        assertEquals(0.0, magnitude, DELTA);
    }
    
    @Test
    @DisplayName("Normalization should create unit vector")
    void testNormalization() {
        // Act
        Point3D normalized = positivePoint.normalize();
        
        // Assert
        assertEquals(1.0, normalized.magnitude(), DELTA);
        
        // Verify direction is preserved (same ratios)
        double originalMagnitude = positivePoint.magnitude();
        assertEquals(positivePoint.getX() / originalMagnitude, normalized.getX(), DELTA);
        assertEquals(positivePoint.getY() / originalMagnitude, normalized.getY(), DELTA);
        assertEquals(positivePoint.getZ() / originalMagnitude, normalized.getZ(), DELTA);
    }
    
    @Test
    @DisplayName("Normalizing zero vector should throw exception")
    void testNormalizeZeroVector() {
        // Act & Assert
        IllegalStateException exception = assertThrows(
            IllegalStateException.class,
            () -> origin.normalize()
        );
        assertEquals("Cannot normalize zero vector", exception.getMessage());
        assertTrue(logHandler.containsLevel(Level.SEVERE));
    }
    
    @Test
    @DisplayName("Dot product should be commutative")
    void testDotProductCommutativity() {
        // Act
        double dot1 = positivePoint.dotProduct(negativePoint);
        double dot2 = negativePoint.dotProduct(positivePoint);
        
        // Assert
        assertEquals(dot1, dot2, DELTA);
    }
    
    @Test
    @DisplayName("Dot product with perpendicular vectors should be zero")
    void testDotProductPerpendicular() {
        // Act
        double dotXY = unitX.dotProduct(unitY);
        double dotXZ = unitX.dotProduct(unitZ);
        double dotYZ = unitY.dotProduct(unitZ);
        
        // Assert
        assertEquals(0.0, dotXY, DELTA);
        assertEquals(0.0, dotXZ, DELTA);
        assertEquals(0.0, dotYZ, DELTA);
    }
    
    @Test
    @DisplayName("Dot product with self should equal magnitude squared")
    void testDotProductWithSelf() {
        // Act
        double dotProduct = positivePoint.dotProduct(positivePoint);
        double magnitudeSquared = Math.pow(positivePoint.magnitude(), 2);
        
        // Assert
        assertEquals(magnitudeSquared, dotProduct, DELTA);
    }
    
    @Test
    @DisplayName("Cross product should be anti-commutative")
    void testCrossProductAntiCommutativity() {
        // Act
        Point3D cross1 = positivePoint.crossProduct(negativePoint);
        Point3D cross2 = negativePoint.crossProduct(positivePoint);
        
        // Assert - cross2 should be negative of cross1
        assertEquals(-cross1.getX(), cross2.getX(), DELTA);
        assertEquals(-cross1.getY(), cross2.getY(), DELTA);
        assertEquals(-cross1.getZ(), cross2.getZ(), DELTA);
    }
    
    @Test
    @DisplayName("Cross product of unit vectors should follow right-hand rule")
    void testCrossProductRightHandRule() {
        // Act
        Point3D crossXY = unitX.crossProduct(unitY); // Should be +Z
        Point3D crossYZ = unitY.crossProduct(unitZ); // Should be +X
        Point3D crossZX = unitZ.crossProduct(unitX); // Should be +Y
        
        // Assert
        assertTrue(crossXY.approximatelyEquals(unitZ));
        assertTrue(crossYZ.approximatelyEquals(unitX));
        assertTrue(crossZX.approximatelyEquals(unitY));
    }
    
    @Test
    @DisplayName("Cross product with parallel vectors should be zero")
    void testCrossProductParallel() {
        // Arrange
        Point3D parallelPoint = new Point3D(6.0, 8.0, 10.0); // 2 * positivePoint
        
        // Act
        Point3D cross = positivePoint.crossProduct(parallelPoint);
        
        // Assert
        assertEquals(0.0, cross.magnitude(), DELTA);
    }
    
    @Test
    @DisplayName("Vector operations with null should throw exception")
    void testVectorOperationsWithNull() {
        // Act & Assert
        assertAll(
            () -> assertThrows(IllegalArgumentException.class,
                () -> positivePoint.dotProduct(null)),
            () -> assertThrows(IllegalArgumentException.class,
                () -> positivePoint.crossProduct(null))
        );
        assertTrue(logHandler.containsLevel(Level.SEVERE));
    }
    
    // ============================================================================
    // EQUALITY AND COMPARISON TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Approximately equals should handle floating-point precision")
    void testApproximatelyEquals() {
        // Arrange
        Point3D point1 = new Point3D(1.0, 2.0, 3.0);
        Point3D point2 = new Point3D(1.0000000001, 2.0000000001, 3.0000000001);
        Point3D point3 = new Point3D(1.1, 2.1, 3.1);
        
        // Act & Assert
        assertTrue(point1.approximatelyEquals(point2)); // Within epsilon
        assertFalse(point1.approximatelyEquals(point3)); // Outside epsilon
        assertFalse(point1.approximatelyEquals(null)); // Null check
    }
    
    @Test
    @DisplayName("Equals should work correctly with exact values")
    void testEquals() {
        // Arrange
        Point3D point1 = new Point3D(1.0, 2.0, 3.0);
        Point3D point2 = new Point3D(1.0, 2.0, 3.0);
        Point3D point3 = new Point3D(1.1, 2.0, 3.0);
        
        // Act & Assert
        assertTrue(point1.equals(point2));
        assertTrue(point1.equals(point1)); // Reflexive
        assertFalse(point1.equals(point3));
        assertFalse(point1.equals(null));
        assertFalse(point1.equals("not a point"));
    }
    
    @Test
    @DisplayName("Hash code should be consistent with equals")
    void testHashCode() {
        // Arrange
        Point3D point1 = new Point3D(1.0, 2.0, 3.0);
        Point3D point2 = new Point3D(1.0, 2.0, 3.0);
        
        // Act & Assert
        assertEquals(point1.hashCode(), point2.hashCode());
    }
    
    // ============================================================================
    // STRING REPRESENTATION TESTS
    // ============================================================================
    
    @Test
    @DisplayName("ToString should provide readable representation")
    void testToString() {
        // Act
        String representation = positivePoint.toString();
        
        // Assert
        assertNotNull(representation);
        assertTrue(representation.contains("Point3D"));
        assertTrue(representation.contains("3.000000"));
        assertTrue(representation.contains("4.000000"));
        assertTrue(representation.contains("5.000000"));
    }
    
    // ============================================================================
    // PERFORMANCE AND EDGE CASE TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Operations should handle very large coordinates")
    void testLargeCoordinates() {
        // Arrange
        Point3D largePoint = new Point3D(1e10, 1e10, 1e10);
        
        // Act & Assert - Should not throw exceptions
        assertDoesNotThrow(() -> {
            largePoint.magnitude();
            largePoint.distanceTo(origin);
            largePoint.normalize();
            largePoint.translate(1e5, 1e5, 1e5);
        });
    }
    
    @Test
    @DisplayName("Operations should handle very small coordinates")
    void testSmallCoordinates() {
        // Arrange
        Point3D smallPoint = new Point3D(1e-10, 1e-10, 1e-10);
        
        // Act & Assert - Should not throw exceptions

        assertDoesNotThrow(() -> {
            smallPoint.magnitude();
            smallPoint.distanceTo(origin);
            smallPoint.normalize();
            smallPoint.translate(1e-5, 1e-5, 1e-5);
        });
    }
    
    @Test
    @DisplayName("Chained operations should maintain precision")
    void testChainedOperations() {
        // Arrange
        Point3D original = new Point3D(1.0, 2.0, 3.0);
        
        // Act - Chain multiple operations
        Point3D result = original
            .translate(1.0, 0.0, 0.0)
            .rotateZ(Math.PI / 2)
            .translate(-1.0, 0.0, 0.0)
            .rotateZ(-Math.PI / 2);
        
        // Assert - Should approximately return to original
        assertTrue(original.approximatelyEquals(result));
    }
    
    // ============================================================================
    // LOGGING VERIFICATION TESTS
    // ============================================================================
    
    @Test
    @DisplayName("Methods should generate appropriate log messages")
    void testLoggingBehavior() {
        // Act - Perform various operations
        Point3D point1 = new Point3D(1, 2, 3);
        Point3D point2 = Point3D.origin();
        point1.distanceTo(point2);
        point1.magnitude();
        point1.normalize();
        
        // Assert - Check that appropriate log messages were generated
        assertTrue(logHandler.containsMessage("Created new Point3D"));
        assertTrue(logHandler.containsMessage("Creating Point3D at origin"));
        assertTrue(logHandler.containsMessage("Calculated distance"));
        assertTrue(logHandler.containsMessage("Calculated magnitude"));
        assertTrue(logHandler.containsMessage("Normalizing vector"));
    }
}
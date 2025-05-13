import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Optimized implementation of JSystem using a structure-of-arrays approach
 * with cache-aligned data in direct ByteBuffers.
 */
public class JOptimizedArraySystem implements JSystem {
    // Cache line size (typical for modern processors)
    private static final int CACHE_LINE_SIZE = 64;
    
    // Particle data arrays using direct byte buffers for better memory alignment
    private final ByteBuffer positionBuffer;
    private final ByteBuffer velocityBuffer;
    private final ByteBuffer massRadiusBuffer;
    
    private final int numBodies;
    
    /**
     * Creates a new optimized particle system with the specified number of bodies
     */
    public JOptimizedArraySystem(int nb) {
        this.numBodies = nb;
        
        // Allocate direct ByteBuffers with native byte order for best performance
        // Position: x,y,z for each particle (3 doubles per particle)
        positionBuffer = ByteBuffer.allocateDirect(nb * 3 * Double.BYTES)
                                  .order(ByteOrder.nativeOrder());
        
        // Velocity: vx,vy,vz for each particle (3 doubles per particle)
        velocityBuffer = ByteBuffer.allocateDirect(nb * 3 * Double.BYTES)
                                  .order(ByteOrder.nativeOrder());
        
        // Mass and radius: m,r for each particle (2 doubles per particle)
        // Pad to align with cache lines for better performance
        int paddingPerParticle = (CACHE_LINE_SIZE - (2 * Double.BYTES) % CACHE_LINE_SIZE) % CACHE_LINE_SIZE;
        massRadiusBuffer = ByteBuffer.allocateDirect(nb * (2 * Double.BYTES + paddingPerParticle))
                                   .order(ByteOrder.nativeOrder());
    }

    @Override
    public int numBodies() {
        return numBodies;
    }

    @Override
    public double p(int index, int dim) {
        return positionBuffer.getDouble((index * 3 + dim) * Double.BYTES);
    }

    @Override
    public double v(int index, int dim) {
        return velocityBuffer.getDouble((index * 3 + dim) * Double.BYTES);
    }

    @Override
    public double r(int index) {
        return massRadiusBuffer.getDouble(index * (2 * Double.BYTES + getPadding()) + Double.BYTES);
    }

    @Override
    public double m(int index) {
        return massRadiusBuffer.getDouble(index * (2 * Double.BYTES + getPadding()));
    }

    @Override
    public void init(int index, double x, double y, double z, double vx, double vy, double vz, double rad, double mass) {
        // Set position
        positionBuffer.putDouble((index * 3 + 0) * Double.BYTES, x);
        positionBuffer.putDouble((index * 3 + 1) * Double.BYTES, y);
        positionBuffer.putDouble((index * 3 + 2) * Double.BYTES, z);
        
        // Set velocity
        velocityBuffer.putDouble((index * 3 + 0) * Double.BYTES, vx);
        velocityBuffer.putDouble((index * 3 + 1) * Double.BYTES, vy);
        velocityBuffer.putDouble((index * 3 + 2) * Double.BYTES, vz);
        
        // Set mass and radius
        int offset = index * (2 * Double.BYTES + getPadding());
        massRadiusBuffer.putDouble(offset, mass);
        massRadiusBuffer.putDouble(offset + Double.BYTES, rad);
    }

    @Override
    public void incP(int index, int dim, double delta) {
        int offset = (index * 3 + dim) * Double.BYTES;
        double current = positionBuffer.getDouble(offset);
        positionBuffer.putDouble(offset, current + delta);
    }

    @Override
    public void incV(int index, int dim, double delta) {
        int offset = (index * 3 + dim) * Double.BYTES;
        double current = velocityBuffer.getDouble(offset);
        velocityBuffer.putDouble(offset, current + delta);
    }
    
    /**
     * Calculate padding bytes per particle for proper alignment
     */
    private int getPadding() {
        return (CACHE_LINE_SIZE - (2 * Double.BYTES) % CACHE_LINE_SIZE) % CACHE_LINE_SIZE;
    }
} 
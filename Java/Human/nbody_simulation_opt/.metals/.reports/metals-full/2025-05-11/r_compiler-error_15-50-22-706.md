file://<WORKSPACE>/src/main/java/JClassSystem.java
### java.util.NoSuchElementException: next on empty iterator

occurred in the presentation compiler.

presentation compiler configuration:


action parameters:
uri: file://<WORKSPACE>/src/main/java/JClassSystem.java
text:
```scala
public class JClassSystem implements JSystem {
  static class Particle {
    public Particle(double x, double y, double z, double vx, double vy, double vz, double rad, double mass) {
      p[0] = x;
      p[1] = y;
      p[2] = z;
      v[0] = vx;
      v[1] = vy;
      v[2] = vz;
      r = rad;
      m = mass;
    }
    public double[] p = new double[3];
    public double[] v = new double[3];
    public double r;
    public double m;
  }

  private Particle[] particles;

  public JClassSystem(int nb) {
    particles = new Particle[nb];
  }

  @Override
  public int numBodies() {
    return particles.length;
  }

  @Override
  public double p(int index, int dim) {
    return particles[index].p[dim];
  }

  @Override
  public double v(int index, int dim) {
    return particles[index].v[dim];
  }

  @Override
  public double r(int index) {
    return particles[index].r;
  }

  @Override
  public double m(int index) {
    return particles[index].m;
  }

  @Override
  public void init(int index, double x, double y, double z, double vx, double vy, double vz, double rad, double mass) {
    particles[index] = new Particle(x, y, z, vx, vy, vz, rad, mass);
  }

  @Override
  public void incP(int index, int dim, double delta) {
    particles[index].p[dim] += delta;
  }

  @Override
  public void incV(int index, int dim, double delta) {
    particles[index].v[dim] += delta;
  }
}

```



#### Error stacktrace:

```
scala.collection.Iterator$$anon$19.next(Iterator.scala:973)
	scala.collection.Iterator$$anon$19.next(Iterator.scala:971)
	scala.collection.mutable.MutationTracker$CheckedIterator.next(MutationTracker.scala:76)
	scala.collection.IterableOps.head(Iterable.scala:222)
	scala.collection.IterableOps.head$(Iterable.scala:222)
	scala.collection.AbstractIterable.head(Iterable.scala:935)
	dotty.tools.dotc.interactive.InteractiveDriver.run(InteractiveDriver.scala:164)
	dotty.tools.pc.CachingDriver.run(CachingDriver.scala:45)
	dotty.tools.pc.WithCompilationUnit.<init>(WithCompilationUnit.scala:31)
	dotty.tools.pc.SimpleCollector.<init>(PcCollector.scala:351)
	dotty.tools.pc.PcSemanticTokensProvider$Collector$.<init>(PcSemanticTokensProvider.scala:63)
	dotty.tools.pc.PcSemanticTokensProvider.Collector$lzyINIT1(PcSemanticTokensProvider.scala:63)
	dotty.tools.pc.PcSemanticTokensProvider.Collector(PcSemanticTokensProvider.scala:63)
	dotty.tools.pc.PcSemanticTokensProvider.provide(PcSemanticTokensProvider.scala:88)
	dotty.tools.pc.ScalaPresentationCompiler.semanticTokens$$anonfun$1(ScalaPresentationCompiler.scala:111)
```
#### Short summary: 

java.util.NoSuchElementException: next on empty iterator
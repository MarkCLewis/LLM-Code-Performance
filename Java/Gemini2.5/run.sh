java -Xms8G -Xmx8G \
     -XX:+UseParallelGC \
     -XX:+UseAVX=2 \
     -XX:+AggressiveOpts \
     -XX:+AlwaysPreTouch \
     -Djava.util.concurrent.ForkJoinPool.common.parallelism=8 \
     -cp . NBodySimulationOptimized
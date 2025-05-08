java -Xms8G -Xmx8G \
     -XX:+UseParallelGC \
     -XX:+AlwaysPreTouch \
     -Djava.util.concurrent.ForkJoinPool.common.parallelism=8 \
     -cp . $1
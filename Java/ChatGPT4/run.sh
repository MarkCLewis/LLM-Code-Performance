java \
  -XX:+UseParallelGC \
  -XX:+UseNUMA \
  -XX:+UseCompressedOops \
  -XX:+AlwaysPreTouch \
  -XX:+UnlockExperimentalVMOptions \
  -XX:+UnlockDiagnosticVMOptions \
  -XX:+UseVectorizedMismatchIntrinsic \
  -XX:+OptimizeStringConcat \
  -Xms16G -Xmx16G \
  $1
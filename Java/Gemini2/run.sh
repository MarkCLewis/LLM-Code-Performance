# Suggested run command.
java -server -Xms8g -Xmx16g -XX:+UseG1GC -XX:ParallelGCThreads=16 -XX:+UseNUMA -XX:+AggressiveOpts -XX:+PrintGCDetails NBodySimulation

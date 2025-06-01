import java.io.File

val lineRegex = """([A-Za-z0-9\.]+) ([A-Za-z0-9\._]+)|(real|user)\s+(\d+)m(\d+.\d+)s""".r

def stats(times: List[Double]): (Double, Double) = {
  val avg = times.sum / times.length
  val std = math.sqrt(times.map(t => (t - avg)*(t - avg)).sum / times.length)
  (avg, std)
}

def processTimeFile(file: File): Unit = {
  println(s"$file")
  var author = ""
  var version = ""
  var realSecs: List[Double] = Nil
  var userSecs: List[Double] = Nil
  val source = io.Source.fromFile(file)
  val lines = source.getLines()
  for (case lineRegex(a, v, ru, m, s) <- lines) {
    if (a != null) {
      if (realSecs.nonEmpty) {
        val (ravg, rstd) = stats(realSecs)
        val (uavg, ustd) = stats(userSecs)
        println(f"$author $version \\shortstack[c]{$$$ravg%1.2f \\pm $rstd%1.2f$$ \\\\ $$${uavg/ravg}%1.2f$$}")
        realSecs = Nil
        userSecs = Nil
      }
      author = a
      version = v
    } else {
      if (ru == "real") {
        realSecs ::= 60*m.toDouble + s.toDouble
      } else {
        userSecs ::= 60*m.toDouble + s.toDouble
      }
    }
  }
  val (ravg, rstd) = stats(realSecs)
  val (uavg, ustd) = stats(userSecs)
  println(f"$author $version \\shortstack[c]{$$$ravg%1.2f \\pm $rstd%1.2f$$ \\\\ $$${uavg/ravg}%1.2f$$}")
  source.close()
}

def recurseDirs(dir: File): Unit = {
  for (file <- dir.listFiles()) {
    if (file.isDirectory) {
      recurseDirs(file)
    } else if (file.getName == "times.txt") {
      processTimeFile(file)
    }
  }
}

@main def processTimes(topDir: String): Unit = {
  recurseDirs(new File(topDir))
}

import breeze.linalg._
import breeze.stats.distributions.Beta

///------------
//python code
///------------

//exactCTR = 0.2
//num_iter = 400
//num_skip = 40
//alpha = 5
//beta = 5
//click = 0
//nonclick = 0
//xs = np.linspace(0, 1, 200)
//for i in range(num_iter):
//if np.random.random() < exactCTR:
//  click += 1
//else:
//nonclick += 1
//if i % num_skip == 0:
//  rv = scipy.stats.beta(alpha+click, beta+nonclick)
//ys = rv.pdf(xs)
//plt.plot(xs, ys, label="i = {0}".format(i))
//plt.xlabel("CTR")
//plt.legend()

///------------
//scala code
///------------
val exactCtr = 0.2
val numIter = 400
val numSkip = 40

randomDouble.apply()

import breeze.plot._

val f = Figure()
val p = f.subplot(0)

val xs = linspace(0, 1, 200)

(0 to numIter).foldLeft[(Int, Int)]((5, 5)) { case ((alpha, beta), i) =>

  val (newAlpha, newBeta) =
    if (randomDouble.apply() < exactCtr) (alpha + 1, beta)
    else (alpha, beta + 1)

  if(i % numSkip == 0) {
    val rv = new Beta(newAlpha, newBeta)
    p += plot(xs, xs.map(rv.pdf))
  }
  (newAlpha, newBeta)
}


p.xlabel = "x axis"
p.ylabel = "y axis"
f.saveas("lines.png")
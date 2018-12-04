package ml.em
import breeze.linalg._
import breeze.stats.distributions.MultivariateGaussian
import breeze.numerics._

object Em extends App {
  implicit class MultivariateGaussianOps(self: MultivariateGaussian) {
    def generate(n: Int) = DenseMatrix((1 to n).map(_ => self.draw()): _*)
  }
  implicit class DenseMatrixOps(self: DenseMatrix[Double]) {
    def matmul(that: DenseMatrix[Double]) =
      if (self.cols == that.rows) {
        val c = DenseMatrix.zeros[Double](self.rows, that.cols)
        for {
          i <- 0 until self.rows
          j <- 0 until that.cols
        } {
          var v = 0.0
          for (k <- 0 until self.cols) v += self(i, k) * that(k, j)
          c(i, j) = v
        }
        c
      } else {
        throw new IllegalArgumentException
      }
  }

  val data1 =
    MultivariateGaussian(mean = DenseVector(1000.0, 80.0),
                         covariance =
                           DenseMatrix((10000.0, 2500.0), (2500.0, 10000.0)))
      .generate(200)
  val data2 =
    MultivariateGaussian(mean = DenseVector(200.0, 100.0),
                         covariance =
                           DenseMatrix((12000.0, 1500.0), (1500.0, 12000.0)))
      .generate(100)
  val data3 =
    MultivariateGaussian(mean = DenseVector(500.0, -200.0),
                         covariance =
                           DenseMatrix((3000.0, -500.0), (-500.0, 2000.0)))
      .generate(100)
  val data4 =
    MultivariateGaussian(mean = DenseVector(100.0, 200.0),
                         covariance =
                           DenseMatrix((2000.0, 1500.0), (1500.0, 2000.0)))
      .generate(100)

  val x = DenseMatrix.vertcat(data1, data2, data3, data4)
  println(x)

  val N = x.rows
  println(N)
  val D = x.cols
  println(D)
  val K = 4

//  def f(mu, sigma):
//  def f_(x):
//  return ( (1 / (np.sqrt(2 * math.pi) ** D * np.sqrt(np.linalg.det(sigma)))) *
//    np.exp(-0.5 * np.linalg.solve(sigma, (x - mu)) @ (x - mu)).sum() )
//  return f_
  def f(mu: DenseVector[Double],
        sigma: DenseVector[Double]): DenseVector[Double] => Double = x => {
    val a = sigma.asDenseMatrix \ (x - mu).t
    1 / (pow(sqrt(2 * Math.PI), D) * sqrt(det(sigma.asDenseMatrix))) *
      exp(-0.5 * sum(a matmul (x - mu).asDenseMatrix.t))
  }

//  def estep(x, mu, sigma, pi):
//  #pif = np.array([[pk * f(mk, sk)(xn) for xn in x] for pk, mk, sk in zip(pi, mu, sigma)])
//  pif = np.array([pk * np.apply_along_axis(f(mk, sk), 1, x) for pk, mk, sk in zip(pi, mu, sigma)])
//  gamma = pif / np.sum(pif, axis=0)
//  return gamma
  def estep(x: DenseMatrix[Double],
            mu: DenseMatrix[Double],
            sigma: DenseMatrix[Double],
            pi: DenseMatrix[Double]) = {
    DenseMatrix.zipMap_d
    val pif: DenseMatrix[Double] = {
      for (i <- 0 until K) yield {
        val pk = pi(::, i)
        val mk = mu(i, ::).t
        val sk = sigma(i, ::).t
        val ff = f(mk, sk)
        ff(x(1, ::).t)
      }
      ???
    }

//      zipValues(pi, mu, sigma) { (pk, mk, sk) =>
//      val ff = f(mk, sk)
//      DenseMatrix.vertcat(
//        x(*, ::)
//          .foldLeft(List.empty[DenseMatrix[Double]]) {
//            case (acc, xn) =>
//              pk * ff(xn) :: acc
//          }
//          .reverse: _*)
//    }
    pif / sum(pif, Axis._0)
  }

  val minx = min(x(::, *)).t
  val maxx = max(x(::, *)).t
  val rangex = abs(minx) + abs(maxx)
  val eqdev = rangex / (K + 1.0)

  val mu =
    DenseMatrix((0 until K).map(k => eqdev * (k + 1.0) + minx): _*)

  val sigma =
    DenseMatrix((0 until K).map(k => eqdev * (k + 1.0) + minx): _*)

  val pi = (DenseVector.ones[Double](K) * (1.0 / K)).asDenseMatrix

  val gamma = estep(x, mu, sigma, pi)

}

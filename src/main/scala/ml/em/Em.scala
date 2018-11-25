package ml.em
import breeze.linalg._
import breeze.stats.distributions.MultivariateGaussian
import breeze.numerics._

object Em extends App {
  implicit class MultivariateGaussianOps(self: MultivariateGaussian) {
    def generate(n: Int) = DenseMatrix((1 to n).map(_ => self.draw()) : _*)
  }
  implicit class DenseMatrixOps(self: DenseMatrix[Double]) {
    def matmul(that: DenseMatrix[Double]) =
      if (self.cols == that.rows) {
        val c = DenseMatrix.zeros[Double](self.rows, that.cols)
        for {
          i <- 0 until self.rows
          j <- 0 until that.cols
        }{
          var v = 0.0
          for (k <- 0 until self.cols) v += self(i, k) * that(k, j)
          c(i, j) = v
        }
        c
      } else {
        throw new IllegalArgumentException
      }
  }

  val data1 = MultivariateGaussian(mean = DenseVector(1000.0, 80.0), covariance = DenseMatrix((10000.0, 2500.0),(2500.0, 10000.0))).generate(200)
  val data2 = MultivariateGaussian(mean = DenseVector(200.0, 100.0), covariance = DenseMatrix((12000.0, 1500.0),(1500.0, 12000.0))).generate(100)
  val data3 = MultivariateGaussian(mean = DenseVector(500.0,-200.0), covariance = DenseMatrix(( 3000.0, -500.0),(-500.0,  2000.0))).generate(100)
  val data4 = MultivariateGaussian(mean = DenseVector(100.0, 200.0), covariance = DenseMatrix(( 2000.0, 1500.0),(1500.0,  2000.0))).generate(100)
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
  def f(mu: DenseMatrix[Double], sigma: DenseMatrix[Double]) = (x: DenseMatrix[Double]) => {
    1 / (pow(sqrt(2 * Math.PI), D) * sqrt(det(sigma))) *
      exp(-0.5 * sum((sigma \ (x - mu)) matmul (x - mu)))
  }


}
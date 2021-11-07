import java.util.Random;

public class ExtendedDistributions {
    Random rand;
    Random extraRand;
    double[] PoissonK;
    double[] exp_X;
    final int maxGrid = 2048;
    final double probArea = 1.0d / (double) maxGrid;

    public ExtendedDistributions() {
        rand = new Random();
        exp_X = new double[maxGrid + 1];
        PoissonK = new double[]{
                Math.log(1.0d), Math.log(1.0d),
                Math.log(2.0d), Math.log(6.0d), Math.log(24.0d),
                Math.log(120.0d), Math.log(720.0d), Math.log(5040.0d),
                Math.log(40320.0d), Math.log(362880.0d)};

        // Compute table for exponential distribution
        exp_X[0] = 0.0d;
        exp_X[maxGrid] = 20.0d;
        for (int i = 1; i < maxGrid; i++)
            exp_X[i] = -Math.log(Math.exp(-exp_X[i - 1]) - probArea);

    }

    // Create many gamma variates, similar to the R function.
    public double[] rgamma(int n, double shape, double scale) {
        double[] rgamma = new double[n];
        double v, X, U, j, k = 1;

        if (shape < 1) {
            U = rand.nextDouble();
            shape++;
            k = Math.pow(U, 1.0d / shape);
        }

        final double d = shape - 0.3333333333333333d;
        final double c = 1.0d / (3.0d * Math.sqrt(d));

        // Note: requires that x < Math.sqrt( 9 * a - 3 );
        int i = 0;
        for (; i < n; i++) {
            do {
                do {
                    X = rand.nextGaussian();
                    j = 1.0d + c * X;
                } while (j <= 0);

                v = j * j * j;
                U = rand.nextDouble();

            } while (U > 1.0d - 0.0331 * (X * X * X * X)
                    && Math.log(U) > 0.5d * (X * X) + d * (1 - v + Math.log(v)));
            rgamma[i] = k * d * v * scale;
        }
        return rgamma;
    }

    // Based on
    // A Simple Method for Generating Gamma Variables
    // By Marsaglia / Tsang
    // ACM Transactions on Mathematical Software Vol. 26, No.3
    // September 2000, Pages 363-372
    public double nextGamma(double shape, double scale) {
        double k = 1;
        double v, X, U, j;

        if (shape < 1) {
            U = rand.nextDouble();
            shape++;
            k = Math.pow(U, 1.0d / shape);
        }

        final double d = shape - 0.333333333333333333333d;
        final double c = 1.0d / (3.0d * Math.sqrt(d));

        // Note: requires that x < Math.sqrt( 9 * a - 3 );
        do {
            do {
                X = rand.nextGaussian();
                j = 1.0d + c * X;
            } while (j <= 0);

            v = j * j * j;
            U = rand.nextDouble();

        } while (U > 1.0d - 0.0331 * (X * X * X * X)
                && Math.log(U) > 0.5d * (X * X) + d * (1 - v + Math.log(v)));
        return k * d * v * scale;
    }

    public double nextChiSquare(double DoF) {
        return nextGamma(DoF / 2.0d, 2.0d);
    }

    // This version is probably more accurate, but far slower than nextExponential()
    public double nextSlowExponential() {
        return -Math.log(rand.nextDouble());
    }

    public double nextPareto(double shape, double scale) {
        return scale / Math.pow(rand.nextDouble(), 1 / shape);
    }

    public double nextExpGamma(double shape, double scale) {
        return Math.log(nextGamma(shape, scale));
    }

    public double nextBeta(double alpha, double beta) {
        final double X1 = nextGamma(alpha, 1.0d);
        final double X2 = nextGamma(beta, 1.0d);
        return X1 / (X1 + X2);
    }

    public double nextFisher(double v1, double v2) {
        return (nextChiSquare(v1) / nextChiSquare(v2));
    }

    public double nextWeibull(double lambda, double k) {
        return (lambda * Math.pow(-Math.log(1.0d - rand.nextDouble()), 1.0d / k));
    }

    public double nextGumbel(double location, double scale) {
        return (location - scale * Math.log(-Math.log(rand.nextDouble())));
    }

    public double nextExtremeValue(double a, double b) {
        return a - b * Math.log(-Math.log(1.0d - rand.nextDouble()));
    }

    public double nextGeometric(double p) {
        return (Math.log(rand.nextDouble()) / Math.log((1.0d - p)));
    }

    public double nextLogistic(double mu, double s) {
        return (mu - s * Math.log(nextExponential() / nextExponential()));
    }

    public double nextCauchy(double mu, double sigma) {
        return (sigma * nextStudentT(1) + mu);
    }

    public double nextPareto(double a) {
        return (Math.pow(rand.nextDouble(), -(1.0d / a)));
    }

    public double nextLogNormal(double mu, double sigma) {
        return (Math.exp((mu + rand.nextGaussian() * sigma)));
    }

    // Based on the method from the paper:
    // Polar generation of random variates with the t-distribution
    // By Ralph W. Bailey, Mathematics of Computation, Volume 62, Number 206
    // April 1994, Pages 779-781
    public double nextStudentT(double DoF) {
        double U, V, W, C, R;

        do {
            U = rand.nextDouble();
            V = extraRand.nextDouble();
            U = (U + U) - 1.0d;
            V = (V + V) - 1.0d;
            W = (U * U) + (V * V);
        } while (W > 1);

        C = U / Math.sqrt(W);
        R = Math.sqrt(DoF * (Math.pow(W, -2 / DoF) - 1.0d));
        return R * C;
    }

    // PTRD Implemented from
    // The transformed rejection method for generating Poisson random variables
    // By W. Hörmann, Insurance: Mathematics and Economics 12
    // 1993, Pages 39-45

    public int nextPoisson(int lambda) {
        if (lambda >= 10) {
            double smu, b, a, alphaInv, vr, V, U, us, k;
            smu = Math.sqrt((double) lambda);
            b = 0.931d + 2.53d * smu;
            a = -0.059d + 0.02483d * b;
            alphaInv = 1.1239d + 1.1328d / (b - 3.4d);
            vr = 0.9277d - 3.6224 / (b - 2.0d);
            for (; ; ) {
                V = rand.nextDouble();
                if (V <= 0.86d * vr) {
                    U = (V / vr) - 0.43d;
                    return (int) (((a + a) / (0.5d - Math.abs(U)) + b) * U + (double) lambda + 0.445d);
                }

                if (V >= vr)
                    U = rand.nextDouble() - 0.5d;
                else {
                    U = V / vr - 0.93d;
                    if (U >= 0.0d)
                        U = 0.5d - U;
                    else U = -0.5d - U;
                    V = rand.nextDouble() * vr;
                }

                us = 0.5 - Math.abs(U);
                if (us < 0.013d && V > us)
                    continue;
                k = (((a + a) / us + b) * U + (double) lambda + 0.445d);
                V = V * (alphaInv) / (a / (us * us) + b);
                if (k >= 10 && Math.log(V * smu) <= (k + 0.5d) * Math.log((double) lambda / k)
                        - lambda - Math.log(2.506628274631000241612) + k -
                        (0.08333333333333333333d - 1.0d / (360.0d * k * k)) / k)
                    return (int) k;

                if (0 <= k && k <= 9 && Math.log(V) <= k * Math.log((double) lambda) -
                        (double) lambda - PoissonK[(int) k])
                    return (int) k;
            }
        } else {
            final double expLambda = Math.exp(-(double) lambda);
            double U, V = 1.0d;
            int rpois = -1;
            do {
                U = rand.nextDouble();
                V *= U;
                rpois++;
            } while (V > expLambda);
            return rpois;
        }
    }

    public int nextNegativeBinomial(int r, double p) {
        return nextPoisson((int) nextGamma((double) r, (1.0d - p) / p));
    }

    // Exponential distribution based on table lookup
    //
    // The exponential distribution is divided into n parts
    // all of equal probability. A variable can be generated
    // with equal probability from each of these regions. The variable
    // is randomized further with a uniform distribution within each region.
    // The method can be sped up considerably by only using one random draw.
    public double nextExponential() {
        double U = rand.nextDouble();
        int index = 1 + rand.nextInt(maxGrid - 1);

        return exp_X[index] + (exp_X[index + 1] - exp_X[index]) * U;
    }
}

// EXP-006 的數值一致性對照工具。
//
// 直接使用 amber_v4.3.jar 內 shade 的 commons-math3 3.6.1，對一組 (n, p, k) 印出
// BinomialDistribution.cumulativeProbability 的**位元表示**，供 C++ 端逐位元比對。
//
// 用法：
//   javac -cp <amber_v4.3.jar> -d <outdir> CdfConformance.java
//   java  -cp <amber_v4.3.jar>:<outdir> CdfConformance > java_cdf.tsv
//
// 取樣範圍刻意涵蓋實際會用到的組合：深度 n 取自真實資料的分佈、
// p 取自 peak search 的格點（含其一半，因為 hetPeak 用 level/2），k 掃過 0..n。
import org.apache.commons.math3.distribution.BinomialDistribution;

public class CdfConformance
{
    // k 的取樣：實際用到的是 cumulativeProbability(2)，以及 k = min(alt, n-alt) <= n/2。
    // 全掃 0..n 會產生過多列，這裡取 0..min(n, 40) 全掃再加上 n/2 附近與 n 邊界。
    static int[] sampleK(int n)
    {
        java.util.TreeSet<Integer> ks = new java.util.TreeSet<>();
        for(int k = 0; k <= Math.min(n, 40); k++) ks.add(k);
        for(int d = -2; d <= 2; d++) { int k = n / 2 + d; if(k >= 0 && k <= n) ks.add(k); }
        for(int d = 0; d <= 2; d++) { int k = n - d; if(k >= 0) ks.add(k); }
        int[] result = new int[ks.size()];
        int i = 0;
        for(int k : ks) result[i++] = k;
        return result;
    }

    public static void main(String[] args)
    {
        System.out.println("n\tp\tk\tcdfBits\tcdf");

        double step = 0.001;
        double current = 0.005;
        java.util.List<Double> levels = new java.util.ArrayList<>();
        while(current <= 0.37 + 0.0001)
        {
            double v = Math.round(current * 1000.0) / 1000.0;
            levels.add(v);
            levels.add(v / 2);
            current += step;
            step *= 1.05;
        }

        int[] depths = {1, 2, 3, 5, 8, 13, 21, 25, 34, 50, 75, 100, 150, 200, 317, 500, 1000};

        for(int n : depths)
        {
            for(double p : levels)
            {
                for(int k : sampleK(n))
                {
                    double cdf = new BinomialDistribution(n, p).cumulativeProbability(k);
                    System.out.println(n + "\t" + Double.toHexString(p) + "\t" + k + "\t"
                            + Long.toHexString(Double.doubleToRawLongBits(cdf)) + "\t" + cdf);
                }
            }
        }
    }
}

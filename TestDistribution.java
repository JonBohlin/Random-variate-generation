public class TestDistribution {
    public static void main(String[] args) {
        int nrolls = 10;
        ExtendedDistributions e = new ExtendedDistributions();
        for(int i = 0; i < nrolls; i++)
            System.out.println(e.nextPoisson( 10 ));
    }
}

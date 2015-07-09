namespace Free.Ports.Proj4.Approximation
{
	// Approximation structures and procedures
	public class Tseries		// Chebyshev or Power series structure
	{
		public projUV a, b;		// power series range for evaluation
								// or Chebyshev argument shift/scaling
		public PW_COEF[] cu, cv;
		public int mu, mv;		// maximum cu and cv index (+1 for count)
		public int power;		// != 0 if power series, else Chebyshev
	}
}

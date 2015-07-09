namespace Free.Ports.Proj4.FactorsDeriv
{
	public struct FACTORS
	{
		public DERIVS der;
		public double h, k;			// meridinal, parallel scales
		public double omega, thetap;// angular distortion, theta prime
		public double conv;			// convergence
		public double s;			// areal scale factor
		public double a, b;			// max-min scale error
		public int code;			// info as to analytics, see following
	}
}

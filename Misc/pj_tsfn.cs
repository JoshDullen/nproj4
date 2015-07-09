using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// determine small t
		public static double pj_tsfn(double phi, double sinphi, double e)
		{
			sinphi*=e;
			return Math.Tan(0.5*(HALFPI-phi))/Math.Pow((1.0-sinphi)/(1.0+sinphi), 0.5*e);
		}
	}
}

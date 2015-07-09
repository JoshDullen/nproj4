using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// determine constant small m
		public static double pj_msfn(double sinphi, double cosphi, double es)
		{
			return cosphi/Math.Sqrt(1.0-es*sinphi*sinphi);
		}
	}
}

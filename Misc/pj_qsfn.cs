using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// determine small q
		public static double pj_qsfn(double sinphi, double e, double one_es)
		{
			if(e>=EPS7)
			{
				double con=e*sinphi;
				return one_es*(sinphi/(1.0-con*con)-(0.5/e)*Math.Log((1.0-con)/(1.0+con)));
			}

			return sinphi+sinphi;
		}
	}
}

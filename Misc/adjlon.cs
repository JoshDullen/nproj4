using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// reduce argument to range +/- PI
		public static double adjlon(double lon)
		{
			const double SPI=3.14159265359;
			const double ONEPI=PI;

			if(Math.Abs(lon)<=SPI) return lon;
			lon+=ONEPI;							// adjust to 0..2pi rad
			lon-=TWOPI*Math.Floor(lon/TWOPI);	// remove integral # of 'revolutions'
			lon-=ONEPI;							// adjust back to -pi..pi rad
			return lon;
		}
	}
}

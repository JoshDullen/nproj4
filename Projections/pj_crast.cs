using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_crast : PJ
	{
		public override string Name { get { return "crast"; } }
		public override string DescriptionName { get { return "Craster Parabolic (Putnins P4)"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double XM=0.97720502380583984317;
		const double RXM=1.02332670794648848847;
		const double YM=3.06998012383946546542;
		const double RYM=0.32573500793527994772;
		const double THIRD=0.333333333333333333;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			lp.phi*=THIRD;
			xy.x=XM*lp.lam*(2.0*Math.Cos(lp.phi+lp.phi)-1.0);
			xy.y=YM*Math.Sin(lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=3.0*Math.Asin(xy.y*RYM);
			lp.lam=xy.x*RXM/(2.0*Math.Cos((lp.phi+lp.phi)*THIRD)-1);

			return lp;
		}

		public override PJ Init()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}

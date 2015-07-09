using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_mill : PJ
	{
		public override string Name { get { return "mill"; } }
		public override string DescriptionName { get { return "Miller Cylindrical"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=lp.lam;
			xy.y=Math.Log(Math.Tan(Proj.FORTPI+lp.phi*0.4))*1.25;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.lam=xy.x;
			lp.phi=2.5*(Math.Atan(Math.Exp(0.8*xy.y))-Proj.FORTPI);

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

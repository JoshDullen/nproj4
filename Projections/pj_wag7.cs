using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_wag7 : PJ
	{
		public override string Name { get { return "wag7"; } }
		public override string DescriptionName { get { return "Wagner VII"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		// sphere
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.y=0.90630778703664996*Math.Sin(lp.phi);
			double theta=Math.Asin(xy.y);
			double ct=Math.Cos(theta);
			lp.lam/=3.0;
			double D=1/(Math.Sqrt(0.5*(1+ct*Math.Cos(lp.lam))));
			xy.x=2.66723*ct*Math.Sin(lp.lam);
			xy.y*=1.24104*D;
			xy.x*=D;

			return xy;
		}

		public override PJ Init()
		{
			fwd=s_forward;
			inv=null;
			es=0.0;

			return this;
		}
	}
}

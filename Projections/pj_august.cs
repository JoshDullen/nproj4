using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_august : PJ
	{
		public override string Name { get { return "august"; } }
		public override string DescriptionName { get { return "August Epicycloidal"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double t=Math.Tan(0.5*lp.phi);
			double c1=Math.Sqrt(1.0-t*t);
			lp.lam*=0.5;
			double c=1.0+c1*Math.Cos(lp.lam);
			double x1=Math.Sin(lp.lam)*c1/c;
			double y1=t/c;
			double x12=x1*x1;
			double y12=y1*y1;

			xy.x=x1*(3.0+x12-3.0*y12)*4/3;
			xy.y=y1*(3.0+3.0*x12-y12)*4/3;

			return xy;
		}

		public override PJ Init()
		{
			inv=null;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}

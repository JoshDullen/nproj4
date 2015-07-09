using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_denoy : PJ
	{
		public override string Name { get { return "denoy"; } }
		public override string DescriptionName { get { return "Denoyer Semi-Elliptical"; } }
		public override string DescriptionType { get { return "PCyl, Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double C0=0.95;
		const double C1=-0.08333333333333333333;
		const double C3=0.00166666666666666666;
		const double D1=0.9;
		const double D5=0.03;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.y=lp.phi;
			xy.x=lp.lam;
			lp.lam=Math.Abs(lp.lam);
			xy.x*=Math.Cos((C0+lp.lam*(C1+lp.lam*lp.lam*C3))*(lp.phi*(D1+D5*lp.phi*lp.phi*lp.phi*lp.phi)));

			return xy;
		}

		public override PJ Init()
		{
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}

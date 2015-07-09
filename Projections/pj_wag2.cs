using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_wag2 : PJ
	{
		public override string Name { get { return "wag2"; } }
		public override string DescriptionName { get { return "Wagner II"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double C_x=0.92483;
		const double C_y=1.38725;
		const double C_p1=0.88022;
		const double C_p2=0.88550;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			lp.phi=Proj.aasin(ctx, C_p1*Math.Sin(C_p2*lp.phi));
			xy.x=C_x*lp.lam*Math.Cos(lp.phi);
			xy.y=C_y*lp.phi;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/C_y;
			lp.lam=xy.x/(C_x*Math.Cos(lp.phi));
			lp.phi=Proj.aasin(ctx, Math.Sin(lp.phi)/C_p1)/C_p2;

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
